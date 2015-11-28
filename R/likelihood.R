new.ctmm <- methods::setClass("ctmm", representation(info="list"), contains="list")
new.covm <- methods::setClass("covm", representation(par="numeric",isotropic="logical"), contains="matrix")
# try putting contains as unnamed first element of representation?

#######################################
# convenience wrapper for new.ctmm
ctmm <- function(tau=NULL,isotropic=FALSE,error=FALSE,...)
{
  tau.names <- c("position","velocity","acceleration")
  dim.names <- c("x","y")
  List <- list(...)
  
  info <- List$info
  List$info <- NULL
  if(is.null(info)) { info=list() }
  
  List$error <- error
  
  # put covariance into universal format
  if(!is.null(List$sigma)) { List$sigma <- covm(List$sigma,isotropic=isotropic) }
  List$isotropic <- isotropic
  
  # label spatial elements
  if(!is.null(List$mu)) { List$mu <- as.numeric(List$mu) ; names(List$mu) <- dim.names }
  if(!is.null(List$COV.mu)) { dimnames(List$COV.mu) <- list(dim.names,dim.names) }
  
  # label tau elements
  K <- length(tau)
  if(length(tau)>0)
  {
    tau <- sort(tau,decreasing=TRUE)
    names(tau) <- tau.names[1:K]
    tau <- tau[tau>0]
  }
  List$tau <- tau
  # label tau covariance
  tau.names <- c("area",names(tau))
  if(!is.null(List$COV.tau)) { dimnames(List$COV.tau) <- list(tau.names,tau.names) }
  
  result <- new.ctmm(List,info=info)
  
  return(result)
}

# 2D covariance matrix universal format
covm <- function(pars,isotropic=FALSE)
{
  if(is.null(pars))
  { return(NULL) }
  else if(class(pars)=="covm")
  { return(pars) }
  else if(length(pars)==1)
  {
    pars <- c(pars,0,0)
    sigma <- diag(pars[1],2)
  }
  else if(length(pars)==3)
  { sigma <- sigma.construct(pars) }
  else if(length(pars)==4)
  {
    sigma <- pars
    pars <- sigma.destruct(sigma)
  }
  
  # isotropic error check
  if(isotropic)
  {
    pars <- c(mean(diag(sigma)),0,0)
    sigma <- diag(pars[1],2)
  }
  
  name <- c("x","y")
  dimnames(sigma) <- list(name,name)
  
  names(pars) <- c("area","eccentricity","angle")
  
  new.covm(sigma,par=pars,isotropic=isotropic)
}

# construct covariance matrix from 1-3 parameters
sigma.construct <- function(pars)
{
  GM <- pars[1]
  if(length(pars)==1)
  {
    e <- 0
    theta <- 0
  }
  else
  {
    e <- pars[2]
    theta <- pars[3]
  }
  
  u <- c(cos(theta),sin(theta))
  v <- c(-sin(theta),cos(theta))
  e <- exp(e/2)
  
  sigma <- GM * ( (u%o%u)*e + (v%o%v)/e )
  
  return(sigma)
} 

# reduce covariance matrix to 1-3 parameters
sigma.destruct <- function(sigma)
{
  stuff <- eigen(sigma)
  
  e <- stuff$values
  GM <- sqrt(prod(e))
  e <- log(e[1]/e[2])
  
  if(e==0)
  { theta <- 0 }
  else
  {
    theta <- stuff$vectors[,1]
    theta <- atan(theta[2]/theta[1])
  } 
  
  return(c(GM,e,theta))
}

# blank generic function
pars <- function(...) { return(NA) }

# return the canonical parameters of a covariance matrix
pars.covm <- function(COVM)
{
  if(COVM@isotropic)
  { return(COVM@par[1]) }
  else
  { return(COVM@par) }
}

# returns the canonical parameters of a tau vector
pars.tauv <- function(tau)
{
  if(length(tau)==0)
  { return(tau) }
  else if(tau[1] < Inf)
  { return(tau) }
  else if(length(tau)==1)
  { return(NULL) }
  else
  { return(tau[-1]) }
}

############################
# coarce infinite parameters into finite parameters appropriate for numerics
###########################
ctmm.prepare <- function(data,CTMM)
{
  K <- length(CTMM$tau)  # dimension of hidden state per spatial dimension
  
  range <- TRUE
  
  if(K>0)
  { 
    # numerical limit for rangeless processes
    if(CTMM$tau[1]==Inf)
    {
      range <- FALSE
      CTMM$tau[1] <- log(2^(52/2))*(last(data$t)-data$t[1])
      
      # diffusion -> variance
      if(!is.null(CTMM$sigma))
      { 
        T <- (CTMM$tau[1]-if(K>1){CTMM$tau[2]}else{0})
        CTMM$sigma <- CTMM$sigma*T
        CTMM$sigma@par[1] <- CTMM$sigma@par[1]*T
      }
    }
    
    # continuity reduction
    CTMM$tau = CTMM$tau[CTMM$tau!=0]
    K <- length(CTMM$tau)
  }
  # I am this lazy
  if(K==0) { K <- 1 ; CTMM$tau <- 0 }
  
  CTMM$range <- range
  return(CTMM)
}

ctmm.repair <- function(CTMM)
{
  if(!CTMM$range)
  {
    K <- length(CTMM$tau)
    
    # variance -> diffusion
    T <- (CTMM$tau[1]-if(K>1){CTMM$tau[2]}else{0})
    CTMM$sigma <- CTMM$sigma/T
    CTMM$sigma@par[1] <- CTMM$sigma@par[1]/T
    
    CTMM$tau[1] <- Inf
    
    # delete garbate estimates
    CTMM$mu <- NULL
    CTMM$COV.mu <- NULL
    CTMM$DOF.mu <- NULL
  }
  
  return(CTMM)
}

## prepare error array
error.prepare <- function(DATA,CTMM)
{
  n <- length(DATA$t)
  
  # model the error
  if(CTMM$error>0)
  {
    # is the data supplied with error estimates
    error <- DATA$error
    # if not, then use the modeled error
    if(is.null(error)) { error <- rep(as.numeric(CTMM$error),n) }
  }
  else
  { error <- rep(0,n) }
  
  return(error)
}

###############################
# Propagator/Green's function and Two-time correlation from Langevin equation for Kalman filter and simulations
langevin <- function(dt,CTMM)
{
  tau <- CTMM$tau
  sigma <- methods::getDataPart(CTMM$sigma)
  K <- length(tau)
  
  if(K==0)
  {
    Green <- array(0,c(1,1))
    Sigma <- array(1,c(1,1))
  }
  else if(K==1)
  {
    c0 <- exp(-dt/tau)
    Green <- array(c0,c(1,1))
    Sigma <- array(1-c0^2,c(1,1))
  }
  else if(K==2)
  {
    # representation good when tau are close
    f <- 1/tau
    Omega2 <- 1/prod(tau)
    T <- sum(tau)
    
    if(dt==Inf) # make this numerically relative in future
    {
      Green <- rbind( c(0,0) , c(0,0) )
      Sigma <- rbind( c(1,0) , c(0,Omega2) )
    }
    else
    {
      if(tau[1]>tau[2])
      {
        Exp <- exp(-dt/tau)/diff(tau)
        c0 <- diff(Exp*tau)
        c1 <- -diff(Exp)
        c2 <- diff(Exp/tau)
      }
      else
      {
        nu <- (f[1]-f[2])/2
        f <- (f[1]+f[2])/2
        
        Exp <- exp(-f*dt)
        SinhE <- sinh(nu*dt)*Exp
        CoshE <- cosh(nu*dt)*Exp
        
        c0 <- CoshE + (f/nu)*SinhE
        c1 <- -(Omega2/nu)*SinhE
        c2 <- -Omega2*(CoshE - (f/nu)*SinhE)
      }
        
      Green <- rbind( c(c0,-c1/Omega2) , c(c1,-c2/Omega2) )
      Sigma <- -T*c1^2  #off-diagonal term
      Sigma <- rbind( c(1,0)-c(c0^2+c1^2/Omega2,Sigma) , c(0,Omega2)-c(Sigma,c1^2+c2^2/Omega2) )
    }
  }
  
  return(list(Green=Green, Sigma=sigma*Sigma))
}

#############################################################
# Internal Kalman filter/smoother for multiple derivatives, dimensions, trends
# Kalman filter/smoother for matrix P operator and multiple mean functions
# this is a lower level function
# more for generalizability/readability than speed at this point
kalman <- function(z,dt,CTMM,error=rep(0,nrow(z)),smooth=FALSE)
{
  n <- nrow(z)
  DATA <- 1:ncol(z)
  
  # stationary mean function
  u <- array(1,n)
  z <- cbind(z,u)
  VEC <- ncol(z)
  
  tau <- CTMM$tau
  K <- length(tau)  # dimension of hidden state per spatial dimension
  
  # observed dimensions
  OBS <- 1
  Id <- diag(OBS)
  
  # observable state projection operator (will eventially upgrade to use velocity telemetry)
  P <- array(0,c(K,OBS))
  P[1:OBS,] <- 1
  
  # forecast estimates for zero-mean, unit-variance filters
  zFor <- array(0,c(n,K,VEC))
  sFor <- array(0,c(n,K,K))
  
  # forcast residuals
  zRes <- array(0,c(n,OBS,VEC))
  sRes <- array(0,c(n,OBS,OBS))
  
  # concurrent/filtered estimates
  zCon <- array(0,c(n,K,VEC))
  sCon <- array(0,c(n,K,K))
  
  # initial state info
  Langevin <- langevin(dt=dt[1],CTMM=CTMM)
  Green <- Langevin$Green
  Sigma <- Langevin$Sigma
  sFor[1,,] <- Sigma
  
  for(i in 1:n)
  {
    # residual covariance
    sForP <- sFor[i,,] %*% P # why do I need this?
    sRes[i,,] <- ((t(P) %*% sForP) + error[i]*Id)
    
    # forcast residuals
    zRes[i,,] <- z[i,] - (t(P) %*% zFor[i,,])
    
    if(all(abs(sRes[i,,])<Inf)){ Gain <- sForP %*% solve(sRes[i,,]) }
    else { Gain <- sForP %*% (0*Id) } # solve() doesn't like this case...
    
    # concurrent estimates
    zCon[i,,] <- zFor[i,,] + (Gain %*% zRes[i,,])
    sCon[i,,] <- (sFor[i,,] - (Gain %*% t(sForP)))
    
    # update forcast estimates for next iteration
    if(i<n)
    {
      # does the time lag change values? Then update the propagators.
      if(dt[i+1] != dt[i])
      {
        Langevin <- langevin(dt=dt[i+1],CTMM=CTMM)
        Green <- Langevin$Green
        Sigma <- Langevin$Sigma
      }
      #update forcast estimates now
      zFor[i+1,,] <- Green %*% zCon[i,,]
      sFor[i+1,,] <- ((Green %*% sCon[i,,] %*% t(Green)) + Sigma)
    }
  }
  
  # Profile stationary mean
  mu <- sapply(1:VEC,function(i){sum(zRes[,,i]*zRes[,,VEC]/sRes)})
  W <- mu[VEC]
  mu <- mu[DATA]/W
  
  # returned profiled mean
  if(!smooth)
  {
    # Detrend mean
    zRes[,1,DATA] <- zRes[,1,DATA] - (zRes[,1,VEC] %o% mu)
    
    return(list(W=W,mu=mu,zRes=zRes,sRes=sRes))
  }
  
  # NOW RUN SMOOTHER
  # THIS NEEDS TO BE RE-WRITTEN FROM ERROR CHANGES
  
  # delete residuals
  rm(zRes,sRes)
  
  # Finish detrending the effect of a stationary mean
  MEAN <- zCon[,,VEC] ; dim(MEAN) <- c(n,K)
  zCon[,,DATA] <- zCon[,,DATA,drop=FALSE] - (MEAN %o% mu)
  MEAN <- zFor[,,VEC] ; dim(MEAN) <- c(n,K)
  zFor[,,DATA] <- zFor[,,DATA,drop=FALSE] - (MEAN %o% mu)
  # there has to be a better way to do this?
  # why does R drop dimensions so randomly?

  # delete u(t)
  zCon <- zCon[,,DATA,drop=FALSE]
  zFor <- zFor[,,DATA,drop=FALSE]
  # drop=FALSE must be here for BM/OU and I don't fully understand why
  
  # upgrade concurrent estimates to Kriged estimates  
  Green <- langevin(dt=dt[n],CTMM=CTMM)$Green
  for(i in (n-1):1)
  {
    L <- sCon[i,,] %*% t(Green) %*% solve(sFor[i+1,,])
    
    # overwrite concurrent estimate with smoothed estimate
    zCon[i,,] <- zCon[i,,] + L %*% (zCon[i+1,,]-zFor[i+1,,])
    sCon[i,,] <- sCon[i,,] + L %*% (sCon[i+1,,]-sFor[i+1,,]) %*% t(L)
    
    # next time's propagator if necessary
    if(dt[i] != dt[i+1])
    {
      Green <- langevin(dt=dt[i],CTMM=CTMM)$Green
    }
  }
  
  # restore stationary mean to locations only
  zCon[,1,] <- zCon[,1,] + (u %o% mu)
  
  zname <- c("position")
  if(K>1) { zname <- c(zname,"velocity") }
  dimnames(zCon) <- list(NULL,zname,c("x","y")[DATA])
  dimnames(sCon) <- list(NULL,zname,zname)
  
  # return smoothed states
  # this object is temporary
  state <- list(CTMM=CTMM,Z=zCon,S=sCon)
  class(state) <- "state"
  
  return(state)
}

####################################
# log likelihood function
####################################
ctmm.loglike <- function(data,CTMM=NULL,verbose=FALSE)
{
  n <- length(data$t)
  
  # prepare model for numerics
  CTMM <- ctmm.prepare(data,CTMM)
  
  tau <- CTMM$tau
  K <- length(tau)  # dimension of hidden state per spatial dimension
  
  sigma <- CTMM$sigma
  range <- CTMM$range
  
  ###DECONSTRUCT SIGMA AND TRANSFORM DATA
  
  if(is.null(CTMM$isotropic)) { isotropic <- FALSE }
  else { isotropic <- CTMM$isotropic }
  
  t <- data$t
  x <- data$x
  y <- data$y
  
  n <- length(t)
  
  # time lags
  dt <- c(Inf,diff(t))
  
  error <- error.prepare(data,CTMM)
  
  # do we have to diagonalize the variance?
  if(!CTMM$error)
  {
    CTMM$sigma <- 1
    KALMAN <- kalman(cbind(x,y),dt=dt,CTMM=CTMM,error=error)
    
    zRes <- KALMAN$zRes
    sRes <- KALMAN$sRes
    ML.sigma <- mean(zRes[,,1]*zRes[,,2]/sRes)
    ML.sigma <- cbind( c( mean(zRes[,,1]^2/sRes) , ML.sigma) , c(ML.sigma , mean(zRes[,,2]^2/sRes)) )
    
    # profile variance if sigma unspecified
    if(is.null(sigma)) { sigma <- covm(ML.sigma,isotropic=isotropic) }
    
    mu <- KALMAN$mu
    DOF.mu <- KALMAN$W
    COV.mu <- methods::getDataPart(sigma)/DOF.mu
    
    loglike <- -sum(log(sRes)) -(n/2)*log(det(2*pi*sigma)) - (n/2)*sum(diag(ML.sigma %*% solve(sigma)))
  }
  else if(CTMM$isotropic)
  {
    CTMM$sigma <- sigma@par[1]
    KALMAN <- kalman(cbind(x,y),dt=dt,CTMM=CTMM,error=error)
    
    mu <- KALMAN$mu
    COV.mu <- diag(1/KALMAN$W,2)
    DOF.mu <- KALMAN$W*CTMM$sigma
    
    zRes <- KALMAN$zRes
    sRes <- KALMAN$sRes
    ML.sigma <- mean(c(zRes[,,1]^2/sRes, zRes[,,2]^2/sRes))
    
    loglike <- -sum(log(sRes)) -(n)*log(2*pi) - (n)*ML.sigma
  }
  else
  {
    #diagonalize data and then run two 1D Kalman filters with separate means
    z <- rbind(x,y)
    GM <- sigma@par[1]
    e <- sigma@par[2]
    theta <- sigma@par[3]
    R <- rotate(-theta)
    z <- R %*% z
    
    # major axis likelihood
    CTMM$sigma <- GM * exp(+e/2)
    KALMAN <- kalman(cbind(z[1,]),dt=dt,CTMM=CTMM,error=error)

    zRes <- KALMAN$zRes
    sRes <- KALMAN$sRes
    ML.sigma <- mean(zRes[,,1]^2/sRes)

    mu <- KALMAN$mu
    COV.mu <- 1/KALMAN$W
    DOF.mu <- KALMAN$W*CTMM$sigma
    
    loglike <- -(1/2)*sum(log(sRes)) -(n/2)*log(2*pi) - (n/2)*ML.sigma
    
    # minor axis likelihood
    CTMM$sigma <- GM * exp(-e/2)
    KALMAN <- kalman(cbind(z[2,]),dt=dt,CTMM=CTMM,error=error)

    zRes <- KALMAN$zRes
    sRes <- KALMAN$sRes
    ML.sigma <- c(ML.sigma,mean(zRes[,,1]^2/sRes))

    mu <- c(mu,KALMAN$mu)
    COV.mu <- c(COV.mu,1/KALMAN$W)
    DOF.mu <- c(DOF.mu,KALMAN$W*CTMM$sigma)

    loglike <- loglike -(1/2)*sum(log(sRes)) -(n/2)*log(2*pi) - (n/2)*ML.sigma[2]
    
    # these quantities are diagonalized matrices
    COV.mu <- diag(COV.mu,2)
    DOF.mu <- diag(DOF.mu,2)
    
    # transform results back
    R <- rotate(+theta)
    mu <- R %*% mu
    COV.mu <- R %*% COV.mu %*% t(R)
    DOF.mu <- R %*% DOF.mu %*% t(R)
  }

  if(verbose)
  {
    # assign variables
    CTMM$sigma <- sigma
    CTMM <- ctmm.repair(CTMM)
    
    if(range)
    {
      CTMM$mu <- mu
      CTMM$COV.mu <- COV.mu
      CTMM$DOF.mu <- DOF.mu
    }
    
    CTMM$loglike <- loglike
    attr(CTMM,"info") <- attr(data,"info")

    return(CTMM)
  }
  else  { return(loglike) }
}


###########################################################
# FIT MODEL WITH LIKELIHOOD FUNCTION (convenience wrapper to optim)
ctmm.fit <- function(data,CTMM=ctmm(),control=list(),...)
{
  n <- length(data$t)
  tau <- CTMM$tau
  sigma <- CTMM$sigma
  CTMM$mu <- NULL # can always profile mu analytically

  isotropic <- CTMM$isotropic
  error <- CTMM$error
  
  K <- length(tau)
  
  # extra conditions for endlessly diffusing processes
  if(K>0 && tau[1]==Inf){ range <- FALSE }
  else { range <- TRUE }

  # parameter indices for non-profiled parameters that we have to numerically optimize
  if(error>0) 
  { SIGMA <- 1:(if(isotropic){1}else{3}) }
  else
  { SIGMA <- NULL } 

  TAU <- NULL
  calPARS <- function()
  {
    if(K+range==1) # UM, BM
    { TAU <<- NULL }
    else
    { TAU <<- length(SIGMA) + 1:(K-(1-range)) }
  }
  calPARS()

  # If error is an error estimate rather than TRUE, and if there is no error annotated, then we will fit error
  if(error>0 && !any(data$error>0))
  { ERROR <- length(SIGMA) + length(TAU) + 1 }
  else
  { ERROR <- NULL }
  
  # numerically fit parameters
  PARS <- length(SIGMA) + length(TAU) + length(ERROR)
  # degrees of freedom, including the mean, variance/covariance, tau, and error model
  k <- (if(range){2}else{0}) + (if(isotropic){1}else{3}) + length(TAU) + length(ERROR)
  
  # OPTIMIZATION GUESS (pars)
  # also construct reasonable parscale
  pars <- NULL
  parscale <- NULL
  calpars <- function()
  {
    pars <<- NULL
    parscale <<- NULL
    
    if(length(SIGMA)>0)
    {
      # need some initial guess...
      if(is.null(sigma)) { sigma <<- covm(stats::cov(cbind(data$x,data$y)),isotropic=isotropic) }
      pars <<- pars.covm(sigma)[SIGMA]
      parscale <<- c(pars[1],1,pi/4)[SIGMA]
    }
    
    if(length(TAU)>0)
    {
      pars <<- c(pars,pars.tauv(tau))
      parscale <<- c(parscale,pars[TAU])
    }
    
    if(length(ERROR)>0)
    {
      #assuming GPS for now
      pars <<- c(pars,error)
      parscale <<- c(parscale,error) 
    }
  }
  calpars()
  
  # OPTIMIZATION FUNCTION (fn)
  # optional argument lengths: TAU, TAU+1, TAU+SIGMA
  fn <- function(p)
  {
    # Fix sigma if provided, up to degree provided
    if(length(SIGMA)==0)
    { sigma <- NULL }
    else
    {
      # write over inherited sigma with any specified parameters
      sigma <- pars.covm(sigma)
      sigma[SIGMA] <- p[SIGMA]
      
      # enforce positivity
      sigma[SIGMA[-3]] <- abs(sigma[SIGMA[-3]])
      
      sigma <- covm(sigma,isotropic)
      
      # for some reason optim jumps to crazy big parameters
      if(sigma@par[1]*exp(abs(sigma@par[2])/2)==Inf) { return(Inf) }
    }
    
    # fix tau from par
    if(length(TAU)==0)
    { 
      if(range) { tau <- NULL }
      else { tau <- Inf }
    }
    else
    {
      tau <- p[TAU]
      
      # enforce positivity
      tau <- abs(tau)
      tau <- sort(tau,decreasing=TRUE)
      
      if(!range) { tau <- c(Inf,tau) }
    }
    
    # fix error from par
    if(length(ERROR)>0)
    {
      error <- p[ERROR]
      
      # enforce positivity
      error <- abs(error)
    }

    CTMM <- ctmm(tau=tau,sigma=sigma,isotropic=isotropic,error=error)
    return(-ctmm.loglike(data,CTMM))
  }
  
  # NOW OPTIMIZE
  if(PARS==0) # EXACT
  {
    # Bi-variate Gaussian || Brownian motion with zero error
    CTMM <- ctmm.loglike(data,CTMM=CTMM,verbose=TRUE)
    
    if(K==0){ DOF <- n ; CTMM$tau <- NULL }
    else { DOF <- n-1 }
    
    GM <- CTMM$sigma@par[1]
    CTMM$COV.tau <- rbind(if(isotropic) { GM^2/DOF } else { GM^2/DOF/2 })
  }
  else # all further cases require optimization
  {
    if(PARS==1) # Brent is the best choice here
    {
      RESULT <- NULL
      # direct attempt that can caputure zero boundary
      ATTEMPT <- stats::optim(par=pars,fn=fn,method="Brent",lower=0,upper=10*pars)
      RESULT <- rbind(RESULT,c(ATTEMPT$par,ATTEMPT$value))
      # log scale backup that can't capture zero boundary
      ATTEMPT <- stats::nlm(function(p){f=fn(pars*exp(p))},p=0,stepmax=log(10))
      RESULT <- rbind(RESULT,c(pars*exp(ATTEMPT$estimate),ATTEMPT$minimum))
      # choose the better estimate
      MIN <- sort(RESULT[,2],index.return=TRUE)$ix[1]
      pars <- RESULT[MIN,1]
    }
    else # Nelder-Mead is generally the safest and is default
    {
      control$parscale <- parscale
      pars <- stats::optim(par=pars,fn=fn,control=control,...)$par
    }
    
    # write best estimates over initial guess
    tau <- abs(pars[TAU])
    if(!range){ tau <- c(Inf,tau) }
    
    # save sigma if numerically optimized
    if(error>0)
    {
      sigma <- pars[SIGMA]
      # expand to all parameters
      sigma <- pars.covm(covm(sigma))
      # enforce positivity
      sigma[-3] <- abs(sigma[-3])
      sigma <- covm(sigma,isotropic=isotropic)
    }
    else
    { sigma <- NULL }
    
    # save error magnitude if modeled
    if(length(ERROR)>0)
    { error <- abs(pars[ERROR]) }
    
    # verbose ML information
    CTMM <- ctmm(tau=tau,sigma=sigma,isotropic=isotropic,error=error)
    CTMM <- ctmm.loglike(data,CTMM,verbose=TRUE)
    sigma <- CTMM$sigma
    
    # calculate area-tau uncertainty
    SIGMA <- 1
    calPARS()
    ERROR <- NULL
    calpars()
    COV.tau <- numDeriv::hessian(fn,pars)
    COV.tau <- PDsolve(COV.tau) # robust inverse
    CTMM$COV.tau <- COV.tau
  }
  
  CTMM$AICc <- (2*k-2*CTMM$loglike) + 2*k*(k+1)/(n-k-1)
  return(CTMM)
}


####################################
# Newton-Raphson iterate to a ctmm model
# to test if optim worked
newton.ctmm <- function(data,CTMM)
{
  tau <- CTMM$tau
  
  # wrapper function to differentiate
  fn <- function(par)
  { 
    # will need to update this for telemetry error
    return(-ctmm.loglike(data,CTMM=ctmm(tau=par)))
  }
  
  D <- numDeriv::grad(fn,tau)
  H <- numDeriv::hessian(fn,tau)
  
  tau <- tau - (H %*% D)

  return(ctmm(tau=tau,info=attr(data,"info")))
}


##############################################
# SIMULATE DATA over time array t
simulate.ctmm <- function(object,nsim=1,seed=NULL,t=c(0),...)
{
  if(!is.null(seed)){ set.seed(seed) }
  
  if(is.null(object$sigma))
  { object$sigma <- diag(1,2) }

  object <- ctmm.prepare(list(t=t),object)
  
  tau <- object$tau
  if(is.null(tau)) { tau = 0 }
  K <- length(tau)
  
  sigma <- object$sigma
  
  mu <- object$mu
  if(is.null(mu)) { mu <- c(0,0) }
  
  # recast as explicitly positive-definite matrix
  sigma <- Matrix::Matrix(sigma,sparse=FALSE,doDiag=FALSE)
  sigma <- as(sigma,"dppMatrix")
  # had to do that so sqrtm doesn't complain
  Lambda <- expm::sqrtm(sigma)
  
  K <- length(tau)
  
  n <- length(t)
  dt <- c(Inf,diff(t)) # time lags
  
  # where we will store the data
  x <- rep(0,times=n)
  y <- rep(0,times=n)
  
  # initial hidden state, for standardized process
  Hx <- rep(0,times=K)
  Hy <- rep(0,times=K)
  
  object$sigma <- 1
  for(i in 1:n)
  {
    # tabulate propagators if necessary
    if((i==1)||(dt[i]!=dt[i-1]))
    {
      Langevin <- langevin(dt=dt[i],CTMM=object)
      Green <- Langevin$Green
      Sigma <- Langevin$Sigma
    }
    
    # standardized process
    Hx <- MASS::mvrnorm(n=1,mu=as.vector(Green %*% Hx),Sigma=Sigma)
    Hy <- MASS::mvrnorm(n=1,mu=as.vector(Green %*% Hy),Sigma=Sigma)
    
    # un-standardize the process
    r <- as.vector( Lambda %*% c(Hx[1],Hy[1]) ) + mu
    
    x[i] <- r[1]
    y[i] <- r[2]
  }
  
  data <- data.frame(t=t, x=x, y=y)
  data <- new.telemetry(data,info=attr(object,"info"))
  
  return(data)
}
#methods::setMethod("simulate",signature(object="ctmm"), function(object,...) simulate.ctmm(object,...))


###################################################
# Calculate good CIs for other functions
confint.ctmm <- function(model, alpha=0.05)
{
  # z-values for low, ML, high estimates
  z <- stats::qnorm(1-alpha/2)
  
  tau <- model$tau
  tau <- tau[tau<Inf]
  K <- length(tau)
  
  COV <- model$COV.tau
  
  par <- NULL
  
  # timescale uncertainty: can hit 0 and Inf
  if(K>0){
    for(k in 1:K)
    {
      par <- rbind(par,c(0,tau[k],0))
      
      # tau normal for lower CI
      par[k,1] <- max(0,tau[k] - z*sqrt(COV[k+1,k+1]))
      
      # lower CI of f==1/tau normal for upper CI of tau
      par[k,3] <- 1/max(0,(1/tau[k] - z*sqrt(COV[k+1,k+1]/tau[k]^4)))
    }
  }
  
  # standard area uncertainty: chi-square
  GM.sigma <- sqrt(det(model$sigma))
  COV <- COV[1,1] 
  # effective DOF derived from ML curvature
  DOF <- 2*GM.sigma^2/COV
  
  par <- rbind(chisq.ci(GM.sigma,COV,alpha),par)
  
  return(par)
}

summary.ctmm <- function(object,level=0.95,level.UD=0.95,...)
{
  CLASS <- class(object)
  if(CLASS=="ctmm")
  { return(summary.ctmm.single(object,level=level,level.UD=level.UD)) }
  else if(CLASS=="list")
  { return(summary.ctmm.list(object,level=level,level.UD=level.UD)) }
}
  

######################################################
# Summarize results
summary.ctmm.single <- function(object, level=0.95, level.UD=0.95, ...)
{
  alpha <- 1-level
  alpha.UD <- 1-level
  
  # z-values for low, ML, high estimates
  z <- stats::qnorm(1-alpha/2)*c(-1,0,1)
  
  tau <- object$tau
  if(length(tau)>0 && tau[1]==Inf){ range <- FALSE }else{ range <- TRUE }
  tau <- tau[tau<Inf]
  K <- length(tau)
  
  AM.sigma <- mean(diag(object$sigma))
  GM.sigma <- sqrt(det(object$sigma))
  
  COV <- object$COV.tau
  
  # where to store unit information
  name <- rep("",K+1)
  scale <- rep(1,K+1)
  
  par <- confint.ctmm(object,alpha=alpha)
  
  # standard area to home-range area
  par[1,] <- -2*log(alpha.UD)*pi*par[1,]
  
  # pretty area units
  unit.list <- unit(par[1,2],"area")
  name[1] <- unit.list$name
  scale[1] <- unit.list$scale
  
  # pretty time units
  if(K>0)
  {
    for(r in 1:K)
    {
      unit.list <- unit(par[r+1,2],"time")
      name[r+1] <- unit.list$name
      scale[r+1] <- unit.list$scale
    }
  }
  
  # will be improved after telemetry errors incorporated
  if(K>1 || (!range && K>0))
  {
    # RMS velocity
    log.rms <- log(AM.sigma/prod(tau))/2
    VAR.log.rms <- 1/4 * (1/c(GM.sigma,-tau)) %*% COV %*% (1/c(GM.sigma,-tau))
    log.rms <- log.rms + z*sqrt(VAR.log.rms)
    rms <- exp(log.rms) # meters/second
    
    # pretty units
    unit.list <- unit(rms[2],"speed")
    name <- c(name,unit.list$name)
    scale <- c(scale,unit.list$scale)
    
    par <- rbind(par,rms)  
  }
  
  # Fix unit choice
  par <- par/scale
  
  par.names <- "area"
  
  tau.names <- c("position","velocity","acceleration")
  if(!range){ tau.names <- tau.names[-1] }
  if(K>0){ par.names <- c(par.names,paste("tau",tau.names[1:K])) }
  if(K>1 || (!range && K>0)){ par.names <- c(par.names,"speed") }
  
  rownames(par) <- paste(par.names," (",name,")",sep="")
  
  colnames(par) <- c("low","ML","high")

  if(!range) { par <- par[-1,] } # delete off "area"
  
  return(par)
}
#methods::setMethod("summary",signature(object="ctmm"), function(object,...) summary.ctmm(object,...))

summary.ctmm.list <- function(object, level=0.95, level.UD=0.95, ...)
{
  IC <- attr(object,"IC")
  ICS <- sapply(object,function(m){m[[IC]]})
  ICS <- ICS - ICS[[1]]
  ICS <- array(ICS,c(length(ICS),1))
  rownames(ICS) <- names(object)
  colnames(ICS) <- paste("d",IC,sep="")
  return(ICS)
}

ctmm.select <- function(data,CTMM,verbose=FALSE,IC="AICc",...)
{
  MODELS <- list()
  NAMES <- c("IID","OU","OUF")

  # step down in autocorrelation degree
  for(K in length(CTMM$tau):0) 
  {
    MOD <- CTMM

    # fix K
    if(K == 0) { MOD$tau <- NULL } else { MOD$tau <- MOD$tau[1:K] }

    # name model
    NAME <- paste(NAMES[K+1],"anisotropic")
    
    # fit model
    MODELS[[NAME]] <- ctmm.fit(data,MOD,...)
    
    # step down in isotropy
    if(CTMM$isotropic==FALSE)
    {
      MOD$isotropic <- TRUE
      MOD$sigma <- covm(MOD$sigma,isotropic=TRUE)
      NAME <- paste(NAMES[K+1],"isotropic")
      MODELS[[NAME]] <- ctmm.fit(data,MOD)
    }
  }

  # sort models by AICc
  ICS <- sapply(MODELS,function(m){m[[IC]]})
  IND <- sort(ICS,method="quick",index.return=TRUE)$ix
  MODELS <- MODELS[IND]
  attr(MODELS,"IC") <- IC
  
  # return everything
  if(verbose) { return(MODELS) }
  else { return(MODELS[[1]]) }
}


#################################################
# SLOW LIKELIHOOD FUNCTION
# not to be exposed to end users
# THIS NEEDS TO BE FIXED WITH ERRORS
ctmm.loglike.slow <- function(data,CTMM)
{  
  t <- data$t
  x <- data$x
  y <- data$y

  tau <- CTMM$tau
  sigma <- CTMM$sigma
  mu <- CTMM$mu
  isotropic <- CTMM$isotropic
  
  K <- length(tau)
  n <- length(t)
  
  # lag matrix
  C <- outer(t,t,"-")
  C <- abs(C)
  
  # now the correlation matrix
  if(K==1)
  { C <- exp(-C/tau) }
  else if(K==2)
  { C <- (tau[1]*exp(-C/tau[1])-tau[2]*exp(-C/tau[2]))/(tau[1]-tau[2]) }
  
  logdetC <- determinant(C,logarithm=TRUE)$modulus
  
  u <- rep(1,n) # stationary mean function
  
  # now the inverse correlation matrix times vectors
  C <- solve(C,cbind(x,y,u))
  
  Cx <- C[,1]
  Cy <- C[,2]
  w <- C[,3]

  W <- sum(w)
  
  if(is.null(mu)) { mu <- c(sum(Cx), sum(Cy))/W }
  
  # detrend mean
  x <- x - mu[1]
  y <- y - mu[2]

  Cx <- Cx - mu[1]*w
  Cy <- Cy - mu[2]*w
  
  if(is.null(sigma))
  {
    sigma <- x %*% Cy
    sigma <- rbind( c(x %*% Cx, sigma) , c( sigma, y %*% Cy) )/n
    sigma <- covm(sigma,isotropic)
  }

  if(length(sigma) == 1) { sigma <- sigma * diag(2)  }
  
  COV.mu <- methods::getDataPart(sigma)/W
  
  loglike <- -logdetC - (n/2)*log(det(2*pi*sigma)) - n

  CTMM <- ctmm(loglike=loglike,tau=tau,sigma=sigma,mu=mu,COV.mu=COV.mu,DOF.mu=W,info=attr(data,"info"))
    
  return(CTMM)
}