new.ctmm <- methods::setClass("ctmm", representation(info="list"), contains="list")

#######################################
# convenience wrapper for new.ctmm
ctmm <- function(tau=NULL,isotropic=FALSE,...)
{
  List <- list(...)
  
  info <- List$info
  List$info <- NULL
  if(is.null(info)) { info=list() }
  
  # variance scalar -> covariance matrix
  if(length(List$sigma)==1)
  {
    List$sigma <- diag(List$sigma,2)
    dimnames(List$sigma) <- list(c("x","y"),c("x","y"))
  }

  # label tau elements
  K <- length(tau)
  tau.names <- c("position","velocity","acceleration")
  if(length(tau)>0) { names(tau) <- tau.names[1:K] }
  List$tau <- tau
  
  List$isotropic <- isotropic
  
  result <- new.ctmm(List,info=info)
  
  return(result)
}

###############################
# Propagator/Green's function and Two-time correlation from Langevin equation for Kalman filter and simulations
langevin <- function(dt,tau)
{
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
  
  return(list(Green=Green, Sigma=Sigma))
}


########################################################
# Kalman filter/smoother for matrix P operator and multiple mean functions
# more for generalizability/readability than speed at this point
kalman <- function(data,model,smooth=FALSE)
{
  t <- data$t
  x <- data$x
  y <- data$y
  n <- length(t)
  
  tau <- model$tau
  K <- length(tau)  # dimension of hidden state per spatial dimension
  
  # time lags
  dt <- c(Inf,diff(t))
  
  # observable state projection operator (will eventially upgrade to use velocity telemetry)
  P <- array(0,c(K,1))
  P[1,] <- 1
  
  # stationary mean function
  u <- array(1,n)
  
  # Observations to run through filter
  z <- cbind(x,y,u)
  
  # forecast estimates for zero-mean, unit-variance filters
  zFor <- array(0,c(n,K,3))
  sFor <- array(0,c(n,K,K))
  
  # forcast residuals
  zRes <- array(0,c(n,1,3))
  sRes <- array(0,c(n,1,1))
  
  # concurrent estimates
  zCon <- array(0,c(n,K,3))
  sCon <- array(0,c(n,K,K))
  
  # initial state info
  Langevin <- langevin(dt=dt[1],tau=tau)
  Green <- Langevin$Green
  Sigma <- Langevin$Sigma
  sFor[1,,] <- Sigma
  
  for(i in 1:n)
  {
    sForP <- sFor[i,,] %*% P # why do I need this?
    sRes[i,,] <- Matrix::symmpart(t(P) %*% sForP) # zero telemetry error
    
    # forcast residuals
    zRes[i,,] <- z[i,] - (t(P) %*% zFor[i,,])
    
    Gain <- sForP %*% solve(sRes[i,,])
    
    # concurrent estimates
    zCon[i,,] <- zFor[i,,] + (Gain %*% zRes[i,,])
    sCon[i,,] <- Matrix::symmpart(sFor[i,,] - (Gain %*% t(sForP)))
    
    # update forcast estimates for next iteration
    if(i<n)
    {
      # does the time lag change values? Then update the propagators.
      if(dt[i+1] != dt[i])
      {
        Langevin <- langevin(dt=dt[i+1],tau=tau)
        Green <- Langevin$Green
        Sigma <- Langevin$Sigma
      }
      #update forcast estimates now
      zFor[i+1,,] <- Green %*% zCon[i,,]
      sFor[i+1,,] <- Matrix::symmpart(Green %*% sCon[i,,] %*% t(Green) + Sigma)
    }
  }
  
  # Profile stationary mean
  W <- sum(zRes[,,3]^2/sRes)
  mu <- c(sum(zRes[,,1]*zRes[,,3]/sRes),sum(zRes[,,2]*zRes[,,3]/sRes))/W
  
  # returned profiled mean
  if(!smooth)
  {
    # Detrend mean
    zRes[,1,1:2] <- zRes[,1,1:2] - (zRes[,1,3] %o% mu)

    return(list(W=W,mu=mu,zRes=zRes,sRes=sRes))
  }
  
  # NOW RUN SMOOTHER

  # delete residuals
  rm(zRes,sRes)
  
  # Finish detrending stationary mean
  zCon[,1,1:2] <- zCon[,1,1:2] - (zCon[,1,3] %o% mu)
  zFor[,1,1:2] <- zFor[,1,1:2] - (zFor[,1,3] %o% mu)

  # delete u(t)
  zCon <- zCon[,,1:2]
  zFor <- zFor[,,1:2]
  
  # upgrade concurrent estimates to Kriged estimates  
  Green <- langevin(dt=dt[n],tau=tau)$Green
  for(i in (n-1):1)
  {
    L <- sCon[i,,] %*% t(Green) %*% solve(sFor[i+1,,])

    zCon[i,,] <- zCon[i,,] + L %*% (zCon[i+1,,]-zFor[i+1,,])
    sCon[i,,] <- sCon[i,,] + L %*% (sCon[i+1,,]-sFor[i+1,,]) %*% t(L)
    
    # next time's propagator if necessary
    if(dt[i] != dt[i+1])
    {
      Green <- langevin(dt=dt[i],tau=tau)$Green
    }
  }
  
  # restore stationary mean to locations only
  zCon[,1,] <- zCon[,1,] + (u %o% mu)
  
  zname <- c("position")
  if(K>1) { zname <- c(zname,"velocity") }
  dimnames(zCon) <- list(NULL,zname,c("x","y"))
  dimnames(sCon) <- list(NULL,zname,zname)
  
  # return smoothed states
  state <- list(model=model,dt=stats::median(dt),t=t,H=zCon,C=sCon)
  class(state) <- "krige"
  
  return(state)
}


####################################
# log likelihood function
####################################
ctmm.loglike <- function(data,CTMM=NULL,verbose=FALSE)
{
  n <- length(data$t)
  
  tau <- CTMM$tau
  K <- length(tau)  # dimension of hidden state per spatial dimension
  
  sigma <- CTMM$sigma
  range <- TRUE
  
  if(K>0)
  { 
    # numerical limit for rangeless processes
    if(tau[1]==Inf)
    {
      range <- FALSE
      tau[1] <- log(2^(52/2))*(last(data$t)-data$t[1])
      # diffusion -> variance
      if(!is.null(sigma)){ sigma <- sigma*(tau[1]-if(K>1){tau[2]}else{0}) }
    }
    
    # continuity reduction
    tau = tau[tau!=0]
    K <- length(tau)
  }
  # I am this lazy
  if(K==0) { K <- 1 ; tau <- 0 }
  CTMM$tau <- tau
  
  # bounds to constrain optim
  if(min(tau)<0){return(-Inf)}
  if(is.unsorted(rev(tau))){return(-Inf)} # should pass to CPF instead of OUF((tau))
  if(!is.null(sigma) && (min(eigen(sigma)$values)<=0)){return(-Inf)}
  
  if(is.null(CTMM$isotropic)) { isotropic <- FALSE }
  else { isotropic <- CTMM$isotropic }
  
  # get residuals from Kalman filter
  Kalman <- kalman(data,CTMM,smooth=FALSE)
  W <- Kalman$W
  mu <- Kalman$mu
  zRes <- Kalman$zRes
  sRes <- Kalman$sRes
    
  # This expression becomes approximate with telemetry error
  ML.sigma <- sum(zRes[,,1]*zRes[,,2]/sRes) # off-diagonal first
  ML.sigma <- rbind( c(sum(zRes[,,1]^2/sRes),ML.sigma) , c(ML.sigma,sum(zRes[,,2]^2/sRes)) )/n
  
  # Profile variance
  if(is.null(sigma))
  {
    sigma <- ML.sigma 
    if(isotropic){ sigma <- diag(mean(diag(sigma)),2) }
  }
  
  loglike <- -sum(log(sRes)) -(n/2)*log(det(2*pi*sigma)) - (n/2)*sum(diag(ML.sigma %*% solve(sigma)))

  if(verbose)
  {
    COV.mu <- sigma/W
    
    if(!range)
    {
      # variance -> diffusion
      sigma <- sigma/(tau[1]-if(K>1){tau[2]}else{0})
      tau[1] <- Inf
    }
    
    CTMM <- list(loglike=loglike,isotropic=isotropic,tau=tau,sigma=sigma,mu=mu,COV.mu=COV.mu,DOF.mu=W)
    CTMM <- new.ctmm(CTMM,info=attr(data,"info"))
    
    return(CTMM)
  }
  else  { return(loglike) }
}

###########################################################
# FIT MODEL WITH LIKELIHOOD FUNCTION (convenience wrapper to optim)
ctmm.fit <- function(data,CTMM=NULL,...)
{
  n <- length(data$t)
  tau <- CTMM$tau
  sigma <- CTMM$sigma
  CTMM$mu <- NULL # can always profile mu analytically
  
  if(is.null(CTMM$isotropic)) { isotropic <- FALSE }
  else { isotropic <- CTMM$isotropic }
  
  K <- length(tau)
  COVstring <- c("GM.sigma")
    
  # extra conditions for endlessly diffusing processes
  if(K>0 && tau[1]==Inf){ range <- FALSE }
  else { range <- TRUE }

  # MLE
  if(K==0 || (K==1 && !range)) # Bi-variate Gaussian & Brownian motion
  {
    CTMM <- ctmm.loglike(data=data,CTMM=CTMM,verbose=TRUE)

    if(K==0){ DOF <- n ; CTMM$tau <- NULL }
    else { DOF <- n-1 }
    
    GM.sigma <- sqrt(det(CTMM$sigma))
    CTMM$COV.tau <- rbind(if(isotropic) { GM.sigma/DOF } else { GM.sigma/DOF/2 })
  }
  else # all further cases require optimization
  {
    # function to optimize tau
    if(range) # OU & OUF
    {
      fn <- function(par){ return(-ctmm.loglike(data=data,CTMM=ctmm(tau=abs(par),isotropic=isotropic))) }
    }
    else # IOU fixes tau[1]=Inf
    {
      fn <- function(par){ return(-ctmm.loglike(data=data,CTMM=ctmm(tau=c(Inf,abs(par)),isotropic=isotropic))) }
    }
    
    if(K>1 && range) # multidimensional optim: OUF
    { CTMM.optim <- stats::optim(par=tau,fn=fn,...) }
    else # 1D optim: OU & IOU
    { CTMM.optim <- stats::optim(par=tau[K],fn=fn,method="Brent",lower=0,upper=100*tau[K],...) }

    # write best estimate over initial guess
    tau <- abs(CTMM.optim$par)
    if(!range){ tau <- c(Inf,tau) }
    # get ML mu and sigma
    CTMM <- ctmm.loglike(data,CTMM=ctmm(tau=tau,isotropic=isotropic),verbose=TRUE)
    sigma <- CTMM$sigma
    
    # covariance area/scale
    GM.sigma <- sqrt(det(sigma))
    # SO(2) orientation matrix
    R <- sigma/GM.sigma
    # wrapper function to estimate area correlation with tau
    # par = c(sigma.GM,tau)
    fn <- function(par)
    {
      CTMM.par <- CTMM
      CTMM.par$sigma <- par[1]*R
      CTMM.par$tau <- abs(par[-1])
      if(!range){ CTMM.par$tau <- c(Inf,CTMM.par$tau) }
      return(-ctmm.loglike(data=data,CTMM=CTMM.par))
    }
    COV.tau <- solve(numDeriv::hessian(fn,c(GM.sigma,if(range){tau}else{tau[-1]})))
    tau.names <- c("position","velocity","acceleration")
    COVstring <- c(COVstring,tau.names[(2-range):K])
    
    CTMM$COV.tau <- COV.tau
    
    #if(min(eigen(COV.tau)$values)<=0){ warning("Maximum likelihood not found.",immediate.=TRUE) }
  }
  
  name <- c("x","y")
  dimnames(CTMM$sigma) <- list(name,name)
  dimnames(CTMM$COV.mu) <- list(name,name)
  names(CTMM$mu) <- name
  
  dimnames(CTMM$COV.tau) <- list(COVstring,COVstring)
  
  k <- K - if(!range){1}else{0} + 2 + if(isotropic){1}else{3}
  AICc <- (2*k-2*CTMM$loglike) + 2*k*(k+1)/(n-k-1)

  CTMM <- c(list(AICc=AICc),CTMM)
  
  CTMM <- new.ctmm(CTMM,info=attr(data,"info"))
  
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
    CTMM.par <- CTMM
    CTMM.par$tau <- par
    # will need to update this for telemetry error
    return(-ctmm.loglike(data=data,CTMM=ctmm(tau=par)))
  }
  
  D <- numDeriv::grad(fn,tau)
  H <- numDeriv::hessian(fn,tau)
  
  tau <- tau - (H %*% D)

  return(ctmm(tau=tau))
}


##############################################
# SIMULATE DATA over time array t
simulate.ctmm <- function(object,nsim=1,seed=NULL,t=c(0),...)
{
  if(!is.null(seed)){ set.seed(seed) }
  
  tau <- object$tau
  if(is.null(tau)) { tau = 0 }
  K <- length(tau)
  
  sigma <- object$sigma
  if(is.null(sigma))
  { sigma <- diag(2) }
  else if(tau[1]==Inf)
  {
    tau[1] <- log(2^(52/2))*(last(t)-t[1])
    # diffusion -> variance
    sigma <- sigma*(tau[1]-if(K>1){tau[2]}else{0})
  }
  
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
  
  for(i in 1:n)
  {
    # tabulate propagators if necessary
    if((i==1)||(dt[i]!=dt[i-1]))
    {
      Langevin <- langevin(dt=dt[i],tau=tau)
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

######################################################
# Summarize results
summary.ctmm <- function(object, alpha=0.05, alpha.HR=0.05, ...)
{
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
  par[1,] <- -2*log(alpha.HR)*pi*par[1,]
  
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


#################################################
# SLOW LIKELIHOOD FUNCTION
# not to be exposed to end users
ctmm.loglike.slow <- function(data,model)
{  
  t <- data$t
  x <- data$x
  y <- data$y

  tau <- model$tau
  sigma <- model$sigma
  mu <- model$mu
  
  if(is.null(model$isotropic)) { isotropic <- FALSE }
  else { isotropic <- model$isotropic }
  
  K <- length(tau)
  n <- length(t)
  
  # lag matrix
  C <- outer(t,t,"-")
  C <- abs(C)
  
  # now the correlation matrix
  if(K==1)
  {
    C <- exp(-C/tau)
  }
  else if(K==2)
  {
    C <- (tau[1]*exp(-C/tau[1])-tau[2]*exp(-C/tau[2]))/(tau[1]-tau[2])
  }
 
  logdetC <- determinant(C,logarithm=TRUE)$modulus
  
  u <- rep(1,n) # stationary mean function
  
  # now the inverse correlation matrix times vectors
  C <- solve(C,cbind(x,y,u))
  
  Cx <- C[,1]
  Cy <- C[,2]
  w <- C[,3]

  W <- sum(w)
  
  if(is.null(mu))
  { mu <- c(sum(Cx), sum(Cy))/W }
  
  # detrend mean
  x <- x - mu[1]
  y <- y - mu[2]

  Cx <- Cx - mu[1]*w
  Cy <- Cy - mu[2]*w
  
  if(is.null(sigma))
  {
    sigma <- rbind( c(x %*% Cx, x %*% Cy) , c( y %*% Cx, y %*% Cy) )/n
    
    if(isotropic)
    { sigma <- sum(diag(sigma))/2 }
  }

  if(length(sigma) == 1)
  { sigma <- sigma * diag(2)  }
  
  COV.mu <- sigma/W
  
  loglike <- -logdetC - (n/2)*log(det(2*pi*sigma)) - n

  model <- ctmm(loglike=loglike,tau=tau,sigma=sigma,mu=mu,COV.mu=COV.mu,DOF.mu=W,info=attr(data,"info"))
    
  return(model)
}