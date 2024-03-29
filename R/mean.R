# drift/mean function models
# guess some parameters and check the model parameter sanity
drift.init <- function(data,CTMM)
{
  AXES <- length(CTMM$axes)
  z <- get.telemetry(data,CTMM$axes)

  # weights from errors
  if(any(CTMM$error>0))
  {
    error <- get.error(data,CTMM,circle=TRUE,calibrate=TRUE)
    w <- clamp(1/error,0,1) # zero-error GPS approximation versus km-error ARGOS
  }
  else
  { w <- rep(1,length(data$t)) }
  # normalize weights
  # w <- w/sum(w)

  drift <- get(CTMM$mean)
  u <- drift(data$t,CTMM)
  # IID estimate for (error + 0 movement) || (movement + 0 error)
  CTMM$mu <- PDsolve(t(u) %*% (w * u)) %*% (t(u) %*% (w * z))

  # don't return variance if no variance
  if(!CTMM$range) { return(CTMM) }

  n <- length(data$t)
  z <- z - (u %*% CTMM$mu)

  V1 <- sum(w)
  V2 <- sum(w^2)

  z <- vapply(1:n,function(i){ w[i] * (z[i,] %o% z[i,]) },diag(1,AXES)) # [AXES,AXES,n]
  dim(z) <- c(AXES^2,n)
  z <- rowSums(z)
  dim(z) <- c(AXES,AXES)
  z <- z/(V1-V2/V1)
  CTMM$sigma <- z

  if(n==1) { CTMM$sigma <- diag(mean(CTMM$mu[1,]^2),nrow=AXES) } # porque no?
  else if(n==2) { CTMM$sigma <- diag(mean(diag(CTMM$sigma)),nrow=AXES) }

  CTMM$sigma <- covm(CTMM$sigma,isotropic=CTMM$isotropic,axes=CTMM$axes)

  return(CTMM)
}

# how do we rescale time
drift.scale <- function(CTMM,time) { CTMM }

# how do we shift the mean by a constant location
# this assumes the first mean term is always the stationary mean value and all others are deviations
stationary.shift <- function(mu,dmu)
{
  mu[1,] <- mu[1,] + dmu
  return(mu)
}

# for when all mean coefficients are centered on the origin and must be shifted
uniform.shift <- function(mu,dmu)
{
  mu <- t(t(mu) + dmu)
  return(mu)
}

# increase number of linear model parameters
drift.refine <- function(CTMM) { list() }

# name of mean function
drift.name <- function(CTMM) { NULL }

# place to put optional summary information
drift.summary <- function(CTMM,level,level.UD) { NULL }

prototype.drift <- function(t,CTMM) { rep(1,length(t)) }
attr(prototype.drift,"init") <- drift.init
attr(prototype.drift,"name") <- drift.name
attr(prototype.drift,"refine") <- drift.refine
attr(prototype.drift,"scale") <- drift.scale
attr(prototype.drift,"shift") <- stationary.shift
attr(prototype.drift,"summary") <- drift.summary
attr(prototype.drift,"velocity") <- function(t,CTMM) { rep(0,length(t)) }

new.drift <- methods::setClass("drift",
              representation("function",
                             energy="function",
                             init="function",
                             is.stationary="function",
                             name="function",
                             refine="function",
                             scale="function",
                             shift="function",
                             speed="function",
                             summary="function",
                             svf="function",
                             velocity="function"),
              prototype=prototype.drift)

########################################
# !!! TODO !!! REDO EVERYTHING LIKE THIS

drift.complexify <- function(CTMM)
{
  if(CTMM$mean=="stationary") { return(list()) }
  fn <- get(paste0(CTMM$mean,".drift.complexify"))
  fn(CTMM)
}

drift.simplify <- function(CTMM)
{
  if(CTMM$mean=="stationary") { return(list()) }
  fn <- get(paste0(CTMM$mean,".drift.simplify"))
  fn(CTMM)
}

#################################
# stationary mean/drift function
#################################
stationary.drift <- function(t,CTMM) { cbind( array(1,length(t)) ) }
stationary.velocity <- function(t,CTMM) { cbind( array(0,length(t)) ) }

stationary.is.stationary <- function(CTMM) { TRUE }

# svf of mean function
stationary.svf <- function(CTMM)
{
  # default
  EST <- function(t) { 0 }
  VAR <- function(t) { 0 }

  return(list(EST=EST,VAR=VAR))
}

# mean square speed function: point estimate and variance
stationary.speed <- function(CTMM) { list(EST=0,VAR=0) }

# <t(U)%*%(U)> and <t(U')%*%(U')>
stationary.energy <- function(CTMM) { list(UU=1,VV=0) }

# combine this all together for convenience
stationary <- new.drift(stationary.drift,
                        is.stationary=stationary.is.stationary,
                        energy=stationary.energy,
                        speed=stationary.speed,
                        svf=stationary.svf,
                        velocity=stationary.velocity)

############################
# mean zero process
############################

zero.drift <- function(t,CTMM) { cbind( array(0,length(t)) ) }

zero.energy <- function(CTMM) { list(UU=1,VV=0) }

zero <- new.drift(zero.drift,
                  is.stationary=stationary.is.stationary,
                  energy=zero.energy,
                  speed=stationary.speed,
                  svf=stationary.svf,
                  velocity=stationary.velocity)

############################
# Periodic drift function
############################
periodic.is.stationary <- function(CTMM) { !sum(CTMM$harmonic) }

# periodic mean/drift function
periodic.drift <- function(t,CTMM)
{
  harmonic <- CTMM$harmonic
  period <- CTMM$period

  # constant term
  U <- stationary.drift(t,CTMM)

  omega <- periodic.omega(CTMM)
  if(length(omega))
  {
    omega <- t %o% omega
    U <- cbind( U , cos(omega) , sin(omega) )
  }

  return(U)
}

# periodic velocity mean
periodic.velocity <- function(t,CTMM)
{
  harmonic <- CTMM$harmonic
  period <- CTMM$period

  # constant term
  U <- stationary.velocity(t,CTMM)

  omega <- periodic.omega(CTMM)
  if(length(omega))
  {
    theta <- omega %o% t # backwards for multiplication by first dimension
    U <- cbind( U , -t(omega*sin(theta)) , t(omega*cos(theta)) )
  }

  return(U)
}

# guess parameters
periodic.init <- function(data,CTMM)
{
  # default period of 1 day
  if(is.null(CTMM$period)) { CTMM$period <- 1 %#% "day" }

  # default harmonics is none
  if(is.null(CTMM$harmonic)) { CTMM$harmonic <- numeric(length(CTMM$period)) }

  if(is.null(CTMM$mu)) { CTMM <- drift.init(data,CTMM) }

  # !!! maybe run generic drift here

  return(CTMM)
}

# how do we rescale time
periodic.scale <- function(CTMM,time)
{
  CTMM$period <- CTMM$period / time

  return(CTMM)
}

# SVF of mean function
periodic.svf <- function(CTMM)
{
  if(all(CTMM$harmonic==0)) { return( stationary.svf(CTMM) ) }

  STUFF <- periodic.stuff(CTMM)
  omega <- STUFF$omega
  A <- STUFF$A
  COV <- STUFF$COV

  EST <- function(t) { sum( 1/4 * A^2 * (1-cos(omega*t) )) }

  VAR <- function(t)
  {
    grad <- 1/2 * A * (1-cos(omega*t))
    return(c(grad %*% COV %*% grad))
  }

  return(list(EST=EST,VAR=VAR))
}


# name of mean function
periodic.name <- function(CTMM)
{
  NAME <- CTMM$harmonic
  NAME <- paste(NAME,collapse=" ")
  NAME <- paste("harmonic",NAME)
  return(NAME)
}

# calculate deterministic mean square speed and its variance
periodic.speed <- function(CTMM)
{
  if(all(CTMM$harmonic==0)) { return(stationary.speed(CTMM)) }

  STUFF <- periodic.stuff(CTMM)
  omega <- STUFF$omega
  A <- STUFF$A
  COV <- STUFF$COV

  if(is.null(omega)) { return( stationary.speed(CTMM) ) }

  EST <- sum(omega^2*A^2)/2
  grad <- omega^2*A
  VAR <- c(grad %*% COV %*% grad)

  return(list(EST=EST,VAR=VAR))
}

# UU and VV terms
periodic.energy <- function(CTMM)
{
  if(all(CTMM$harmonic==0)) { return(stationary.energy(CTMM)) }
  omega <- periodic.omega(CTMM)
  K <- length(omega)
  if(is.null(omega)) { return( stationary.speed(CTMM) ) }

  # potential energy (modulo amplitudes)
  # <1>==1
  # <sin^2>==<cos^2>==1/2
  UU <- c(1,rep(1/2,2*K))
  UU <- diag(UU)

  # kinetic energy(modulo amplitudes)
  VV <- omega^2/2
  VV <- c(0,VV,VV)
  VV <- diag(VV)

  return(list(UU=UU,VV=VV))
}


# calculate rotational indices
periodic.summary <- function(CTMM,level,level.UD)
{
  alpha <- 1-level
  SUM <- NULL

  if(all(CTMM$harmonic==0)) { return(SUM) }

  STUFF <- periodic.stuff(CTMM)
  omega <- STUFF$omega
  A <- STUFF$A
  COV <- STUFF$COV

  # ROTATIONAL VARIANCE INDEX
  ROT <- sum(A^2)/2 # sine and cosine average 1/2
  RAN <- var.covm(CTMM$sigma)
  MLE <- ROT/(ROT+RAN)

  GRAD <- A
  COV.ROT <- c(GRAD %*% COV %*% GRAD)

  if(CTMM$isotropic) { PARS <- c('major','major') } else { PARS <- c('major','minor') }
  COV.RAN <- CTMM$COV[PARS,PARS]
  GRAD <- c(1,1)
  COV.RAN <- c(GRAD %*% COV.RAN %*% GRAD)

  GRAD <- c(RAN,-ROT)/(ROT+RAN)^2
  VAR <- sum(GRAD^2 * c(COV.ROT,COV.RAN))

  CI <- ci.tau(MLE,VAR,alpha=alpha,min=0,max=1)
  CI <- 100*sqrt(rbind(CI))
  rownames(CI) <- "rotation/deviation %"

  SUM <- rbind(SUM,CI)

  # ROTATIONAL SPEED INDEX
  # CIs are too big to be useful
  if(length(CTMM$tau)>1 && CTMM$tau[2])
  {
    CTMM <- get.taus(CTMM)

    ROT <- sum((omega*A)^2)/2 # sine and cosine average 1/2

    GRAD <- omega^2*A
    COV.ROT <- c(GRAD %*% COV %*% GRAD)

    STUFF <- rand.speed(CTMM)
    RAN <- STUFF$MLE
    COV.RAN <- STUFF$COV

    MLE <- ROT/(ROT+RAN)
    GRAD <- c(RAN,-ROT)/(ROT+RAN)^2
    VAR <- sum(GRAD^2 * c(COV.ROT,COV.RAN))

    CI <- ci.tau(MLE,VAR,alpha=alpha,min=0,max=1)
    CI <- 100*sqrt(rbind(CI))
    rownames(CI) <- "rotation/speed %"

    SUM <- rbind(SUM,CI)
  }

  return(SUM)
}

# increase number of harmonics in model
periodic.drift.complexify <- function(CTMM)
{
  period <- CTMM$period
  harmonic <- CTMM$harmonic

  GUESS <- list()

  # dt <- stats::median(diff(data$t))
  # P <- attr(CTMM$mean,"par")$P
  # # max harmonics for Nyquist frequency
  # kmax <- P/(2*dt)

  for(i in 1:length(period))
  {
    GUESS[[length(GUESS)+1]] <- CTMM
    GUESS[[length(GUESS)]]$harmonic[i] <- harmonic[i] + 1
  }

  return(GUESS)
}

# reduce number of harmonics in model
periodic.drift.simplify <- function(CTMM)
{
  period <- CTMM$period
  harmonic <- CTMM$harmonic

  GUESS <- list()

  for(i in 1:length(period))
  {
    if(harmonic[i]>0)
    {
      GUESS[[length(GUESS)+1]] <- CTMM
      GUESS[[length(GUESS)]]$harmonic[i] <- harmonic[i] - 1
    }
  }

  return(GUESS)
}

# list of non-zero frequencies of model
periodic.omega <- function(CTMM)
{
  harmonic <- CTMM$harmonic
  period <- CTMM$period

  omega <- NULL
  for(i in 1:length(harmonic))
  {
    if(harmonic[i] > 0)
    {
      w <- (2*pi/period[i]) * (1:harmonic[i])
      omega <- c( omega , w )
    }
  }

  return(omega)
}

# useful info from periodic mean function parameters
periodic.stuff <- function(CTMM)
{
  harmonic <- CTMM$harmonic
  period <- CTMM$period
  omega <- periodic.omega(CTMM)

  # amplitudes and covariance
  A <- CTMM$mu
  COV <- CTMM$COV.mu

  # total number of harmonics
  k <- length(omega)

  # default uncertainty of none
  if(is.null(COV)) { COV <- 0 }
  # default structure
  COV <- array(COV,c(2,2*k+1,2*k+1,2))

  # delete stationary coefficients
  A <- A[-1,,drop=FALSE]
  COV <- COV[,-1,,,drop=FALSE]
  COV <- COV[,,-1,,drop=FALSE]

  # structure omega like A
  omega <- c(omega,omega)
  omega <- cbind(omega,omega)
  # this is now (2k)*2

  # flatten block-vectors and block-matrices
  A <- array(A,(2*k)*2)
  omega <- array(omega,(2*k)*2)
  COV <- aperm(COV,c(2,1,3,4))
  COV <- array(COV,c((2*k)*2,(2*k)*2))

  return(list(A=A,COV=COV,omega=omega))
}

# combine this all together for convenience
periodic <- new.drift(periodic.drift,
                      is.stationary=periodic.is.stationary,
                      energy=periodic.energy,
                      init=periodic.init,
                      name=periodic.name,
                      refine=periodic.drift.complexify,
                      scale=periodic.scale,
                      speed=periodic.speed,
                      summary=periodic.summary,
                      svf=periodic.svf,
                      velocity=periodic.velocity)

##################################
# change-point mean functions
##################################
# change.point slot is a data.frame with columns: start, stop, state

change.point.is.stationary <- function(CTMM)
{
  CP <- CTMM$change.point.mu # change points
  if(is.null(CP)) { CP <- CTMM$change.point } # default
  nlevels(CP$state)==1
}


# periodic mean/drift function
change.point.drift <- function(t,CTMM,velocity=FALSE)
{
  CP <- CTMM$change.point.mu # change points
  if(is.null(CP)) { CP <- CTMM$change.point } # default

  COL <- 0
  NAMES <- NULL
  for(s in levels(CP$state))
  {
    M <- CTMM[[s]]
    N <- nrow(M$mu)
    COL <- COL + N
    NAMES <- c(NAMES,paste0(s,".",1:N) )
  }

  U <- array(0,c(length(t),COL))
  colnames(U) <- NAMES

  # get drift function for each contributing state
  i1 <- i2 <- 1
  for(j in 1:nrow(CP))
  {
    while(t[i2+1]<CP$stop[j]) { i2 <- i2+1 }
    SUB <- i1:i2

    s <- as.character(CP$state[j])
    M <- CTMM[[s]]

    NAME <- paste0(s,".",1:nrow(M$mu)) # mu column names

    D <- get(M$mean)

    if(!velocity)
    { U[SUB,NAME] <- D$drift(t[SUB],M) }
    else
    { U[SUB,NAME] <- D$velocity(t[SUB],M) }

    i1 <- i2 <- i2+1
  }

  return(U)
}

# periodic velocity mean
change.point.velocity <- function(t,CTMM)
{ change.point.drift(t=t,CTMM=CTMM,velocity=TRUE) }

# guess parameters
change.point.init <- function(data,CTMM)
{
  CP <- CTMM$change.point.mu # change points
  if(is.null(CP)) { CP <- CTMM$change.point } # default

  # get times for each state
  SUBS <- list()
  i1 <- i2 <- 1
  for(j in 1:nrow(CP))
  {
    while(t[i2+1]<CP$stop[j]) { i2 <- i2+1 }
    SUB <- i1:i2

    s <- as.character(CP$state[j])
    SUBS[[s]] <- c(SUBS[[2]],SUB)

    i1 <- i2 <- i2+1
  }

  for(s in levels(CP$state))
  {
    SUB <- data[SUBS[[s]],]
    M <- CTMM[[s]]
    CTMM[[s]] <- get(M$mean)$init(SUB,M)
  }

  return(CTMM)
}

# how do we rescale time
change.point.scale <- function(CTMM,time)
{
  CP <- CTMM$change.point.mu # change points
  if(is.null(CP)) { CP <- CTMM$change.point } # default

  for(s in levels(CP$state))
  {
    M <- CTMM[[s]]
    CTMM[[s]] <- get(M$mean)$scale(M,time)
  }

  return(CTMM)
}

# SVF of mean function
change.point.svf <- function(CTMM,speed=FALSE)
{
  CP <- CTMM$change.point.mu # change points
  if(is.null(CP)) { CP <- CTMM$change.point } # default

  # get mixture weights (for approximate SVF)
  W <- array(0,nlevels(CP$state))
  names(W) <- levels(CP$state)
  for(i in 1:nrow(CP))
  {
    s <- as.character(CP$state[i])
    W[s] <- W[s] + CP$stop[i]-CP$start[i]
  }
  W <- W/sum(W)

  FNS <- list()
  for(s in levels(CP$state))
  {
    M <- CTMM[[s]]

    if(!speed)
    { FNS[[s]] <- get(M$mean)$svf(M) }
    else
    { FNS[[s]] <- get(M$mean)$speed(M) }
  }

  EST <- function(t) { rowSums( sapply(levels(CP$state),function(s){W[s]*FNS[[s]]$EST(t)}) ) }
  VAR <- function(t) { rowSums( sapply(levels(CP$state),function(s){W[s]^2*FNS[[s]]$VAR(t)}) ) }

  return(list(EST=EST,VAR=VAR))
}


# name of mean function
change.point.name <- function(CTMM)
{
  CP <- CTMM$change.point.mu # change points
  if(is.null(CP)) { CP <- CTMM$change.point } # default

  NAME <- NULL
  for(s in levels(CP$state))
  {
    M <- CTMM[[s]]
    N <- get(M$mean)$name(M)
    if(length(N)) { NAME <- paste0(NAME,N,sep="-") }
  }

  return(NAME)
}

# calculate deterministic mean square speed and its variance
change.point.speed <- function(CTMM)
{ change.point.svf(CTMM,speed=TRUE) }

# UU and VV terms
change.point.energy <- function(CTMM)
{
  CP <- CTMM$change.point.mu # change points
  if(is.null(CP)) { CP <- CTMM$change.point } # default

  uu <- vv <- list()
  for(s in levels(CP$state))
  {
    M <- CTMM[[s]]
    STUFF <- get(M$mean)$energy(M)
    uu[[s]] <- STUFF$UU
    vv[[s]] <- STUFF$VV
  }

  # make a giant block-diagonal matrix
  DIM <- sapply(uu,nrow)
  DIM <- sum(DIM)
  UU <- VV <- array(0,c(DIM,DIM))

  # fill matrix
  offset <- 0
  for(i in 1:length(UU))
  {
    DIM <- 1:nrow(uu[[i]])
    DIM <- c(DIM,DIM)
    UU[offset+DIM] <- uu[[i]]
    VV[offset+DIM] <- vv[[i]]
    offset <- offset + nrow(uu[[i]])
  }

  return(list(UU=UU,VV=VV))
}


# calculate rotational indices
change.point.summary <- function(CTMM,level,level.UD)
{
  CP <- CTMM$change.point.mu # change points
  if(is.null(CP)) { CP <- CTMM$change.point } # default

  CI <- NULL
  for(s in levels(CP$state))
  {
    M <- CTMM[[s]]
    CI <- rbind(CI, get(M$mean)$summary(M,level,level.UD) )
  }

  return(CI)
}

# increase number of harmonics in model
change.point.drift.complexify <- function(CTMM,simplify=FALSE)
{
  CP <- CTMM$change.point.mu # change points
  if(is.null(CP)) { CP <- CTMM$change.point } # default

  GUESS <- list()
  for(s in levels(CP$state))
  {
    M <- CTMM[[s]]
    if(!simplify)
    { GUESS <- c(GUESS, get(M$mean)$complexify(M) ) }
    else
    { GUESS <- c(GUESS, get(M$mean)$simplify(M) ) }
  }

  if(!length(GUESS)) { GUESS <- NULL }
  return(GUESS)
}

# reduce number of harmonics in model
change.point.drift.simplify <- function(CTMM)
{ change.point.drift.complexify(CTMM,simplify=TRUE) }

# combine this all together for convenience
change.point <- new.drift(change.point.drift,
                          is.stationary=change.point.is.stationary,
                          energy=change.point.energy,
                          init=change.point.init,
                          name=change.point.name,
                          refine=change.point.drift.complexify,
                          scale=change.point.scale,
                          speed=change.point.speed,
                          summary=change.point.summary,
                          svf=change.point.svf,
                          velocity=change.point.velocity)


#################################
# continuous uniform spline mean functions
#################################
uspline.is.stationary <- function(CTMM) { !sum(CTMM$knot) }

uspline.drift <- function(t,CTMM)
{
  degree <- CTMM$degree
  knot <- CTMM$knot

  STUFF <- uspline.stuff(CTMM)
  tknot <- STUFF$tknot
  dt <- STUFF$dt

  # midpoint / no real grid / degenerate case
  if(knot==1)
  {
    U <- stationary.drift(t,CTMM)
    if(degree==1) { return( U ) }

    U <- cbind( U , (t-tknot) )
    return(U)
  }

  U <- array(0,c(length(t),knot,degree))

  # current starting knot
  k <- 1
  for(i in 1:length(t))
  {
    s <- t[i]
    # iterate knot number k until we are between the appropriate knots
    while(s>tknot[k+1] && k<knot) { k <- k+1 }

    # relative distance betwee knots
    s <- (s - tknot[k])/dt

    if(degree==1) # linear interpolant
    {
      U[i,k,1] <- 1-s
      U[i,k+1,1] <- s
    }
    else # cubic Hermite interpolant
    {
      U[i,k,1] <- (2*s^3 - 3*s^2 + 1)
      U[i,k,2] <- (s^3 - 2*s^2 + s)*dt
      U[i,k+1,1] <- (-2*s^3 + 3*s^2)
      U[i,k+1,2] <- (s^3 - s^2)*dt
    }
  }

  U <- array(U,c(length(t),knot*degree))

  return(U)
}

# initialize default parameters
uspline.init <- function(data,CTMM)
{
  # degree of continuity: 1,2 - location,velocity
  if(is.null(CTMM$degree)) { CTMM$degree <- 1 }
  # default number of knots
  if(is.null(CTMM$knot)) { CTMM$knot <- 1 }
  # default domain of spline grid
  if(is.null(CTMM$domain)) { CTMM$domain <- c(data$t[1],last(data$t)) }

  if(is.null(CTMM$mu)) { CTMM <- drift.init(data,CTMM) }

  return(CTMM)
}

# scale times
uspline.scale <- function(CTMM,time)
{
  knot <- CTMM$knot
  degree <- CTMM$degree

  if(CTMM$degree>1)
  {
    scale <- array(1,c(knot,degree))
    scale[,2] <- time
    scale <- array(scale,knot*degree)

    CTMM$mu[,1] <- CTMM$mu[,1] * scale
    CTMM$mu[,2] <- CTMM$mu[,2] * scale
  }

  return(CTMM)
}

# svf

# refine
uspline.refine <- function(CTMM)
{
  CTMM$knot <- CTMM$knot + 1
  return(list(CTMM))
}


# name of mean function
uspline.name <- function(CTMM)
{
  NAME <- paste("degree",CTMM$degree,"knot",CTMM$knot)
  return(NAME)
}

# ms speed of spline function
uspline.speed <- function(CTMM)
{
  degree <- CTMM$degree
  knot <- CTMM$knot

  STUFF <- uspline.stuff()
  dt <- STUFF$dt

  if(degree==1 && knot==1)
  {
    EST <- 0
    VAR <- 0
  }
  else if(degree==1)
  {
    EST <- (diff(CTMM$mu[,1])^2 + diff(CTMM$mu[,2])^2)/dt^2/(knot-1)
    VAR <- NA
  }
  else
  {
    EST <- NA
    VAR <- NA
  }

  return(list(EST=EST,VAR=VAR))
}

# summary

# calculate spline grid
uspline.stuff <- function(CTMM)
{
  knot <- CTMM$knot
  domain <- CTMM$domain

  if(knot>1)
  {
    # location of knots
    tknot <- seq(domain[1],domain[2],length.out=knot)
    # grid spacing
    dt <- diff(domain)/(knot-1)
  }
  else
  {
    tknot <- mean(domain)
    dt <- diff(domain)
  }

  return(list(tknot=tknot,dt=dt))
}

# convenience function
uspline <- new.drift(uspline.drift,
                     is.stationary=uspline.is.stationary,
                     init=uspline.init,
                     scale=uspline.scale,
                     shift=uniform.shift,
                     refine=uspline.refine,
                     name=uspline.name)

#################################
# piecewise-stationary mean/drift function
#################################
pwstationary.is.stationary <- function(CTMM) { !length(CTMM$breaks) }

pwstationary.drift <- function(t,CTMM)
{
  breaks <- CTMM$breaks
  id <- factor(breaks$id)
  start <- breaks$start
  stop <- breaks$stop

  IDS <- levels(id)
  u <- array(0,c(length(t),length(IDS)))
  colnames(u) <- IDS

  i <- 1
  for(r in 1:nrow(breaks))
  {
    # stationary mode
    while(t[i] <= stop[r] && i <= length(t))
    {
      u[i,id[r]] <- 1
      i <- i + 1
    }

    # linear transition mode
    if(r < nrow(breaks))
    {
      while(t[i] < start[r+1])
      {
        frac <- (t[i]-stop[r])/(start[r+1]-stop[r])
        u[i,id[r]] <- frac
        u[i,id[r+1]] <- 1 - frac
        i <- i + 1
      }
    }
  }

  return(u)
}

############ init handle posixct start & stop

#############
pwstationary.scale <- function(CTMM,time)
{
  CTMM$breaks[,c('start','stop')] <- CTMM$breaks[,c('start','stop')] / time

  return(CTMM)
}

# combine this all together for convenience
pwstationary <- new.drift(pwstationary.drift,
                          is.stationary=pwstationary.is.stationary,
                          shift=uniform.shift,
                          speed=stationary.speed,
                          svf=stationary.svf)
