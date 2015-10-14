# variogram class
new.variogram <- methods::setClass("variogram",representation(info="list"),contains="data.frame")

# extend subset method
subset.variogram <- function(x,...)
{
  info <- attr(x,"info")
  x <- subset.data.frame(x,...)
  x < - droplevels(x)
  new.variogram(x,info=info)
}


# variogram funcion wrapper
variogram <- function(data,dt=NULL,fast=NULL,CI="Markov")
{
  if(length(dt)<2) { return(variogram.dt(data,dt=dt,fast=fast,CI=CI)) }
  
  # calculate a variograms at each dt
  vars <- lapply(dt, function(DT) { variogram.dt(data,dt=DT,fast=fast,CI=CI) } )
  
  # subset each variogram to relevant range of lags
  dt <- c(dt,Inf)
  lag <- vars[[1]]$lag
  vars[[1]] <- vars[[1]][lag<=dt[2],]
  for(i in 1:(length(dt)-1))
  {
    lag <- vars[[i]]$lag
    vars[[i]] <- vars[[i]][(dt[i]<lag)&(lag<=dt[i+1]),]
  }
  
  # coalate
  var <- vars[[1]]
  for(i in 2:(length(dt)-1)) { var <- rbind(var,vars[[i]]) }

  var <- new.variogram(var,info=attr(data,"info"))
    
  return(var)
} 
  

# wrapper for fast and slow variogram codes, for a specified dt
variogram.dt <- function(data,dt=NULL,fast=NULL,CI="Markov")
{
  # intelligently select algorithm
  if(is.null(fast))
  {
    if(length(data$t)<1000) { fast <- FALSE }
    else { fast <- TRUE }
  }
  
  if(fast)
  {
    return(variogram.fast(data=data,dt=dt,CI=CI))
  }
  else
  { 
    return(variogram.slow(data=data,dt=dt,CI=CI))
  }
}

############################
# best initial time for a uniform grid
grid.init <- function(t,dt=stats::median(diff(t)),W=rep(1,length(t)))
{
  cost <- function(t0)
  { 
    grid <- (t-t0)/dt
    return( sum(W*(grid-round(grid))^2) )
  }

  t0 <- stats::optimize(cost,t[1]+c(-1,1)*dt/2)$minimum
  
  return(t0)   
}

############################
# smear data across a uniform grid
gridder <- function(t,x,y,dt)
{
  n <- length(t)
  
  # time lags
  DT <- diff(t)
  
  # default time step
  if(is.null(dt))
  { dt <- stats::median(DT) }
  
  # gap weights to prevent oversampling with coarse dt
  W <- clamp(c(DT[1],DT)/dt) # left weights
  W <- W + clamp(c(DT,DT[n-1])/dt) # + right weights
  W <- W/2 # average left and right
  
  # choose best grid alignment
  t <- t - grid.init(t,dt,W)
  
  # fractional grid index -- starts at >=1
  index <- t/dt
  while(index[1]<1) { index <- index + 1 }
  while(index[1]>=2) { index <- index - 1 }
  
  # uniform lag grid
  n <- ceiling(max(index))
  lag <- seq(0,n-1)*dt
  
  # continuously distribute times over uniform grid
  W.grid <- rep(0,n)
  X.grid <- rep(0,n)
  Y.grid <- rep(0,n)
  for(i in 1:length(t))
  {
    j <- index[i]
    
    if(floor(j)==ceiling(j))
    { # trivial case
      J <- round(j)
      w <- W[i] # total weight
      W.grid[J] <- W.grid[J] + w
      X.grid[J] <- X.grid[J] + w*x[i]
      Y.grid[J] <- Y.grid[J] + w*y[i]
    }
    else
    { # distribution information between adjacent grids
      
      # left grid portion
      J <- floor(j)
      w <- W[i]*(1-(j-J))
      W.grid[J] <- W.grid[J] + w
      X.grid[J] <- X.grid[J] + w*x[i]
      Y.grid[J] <- Y.grid[J] + w*y[i]
      
      # right grid portion
      J <- ceiling(j)
      w <- W[i]*(1-(J-j))
      W.grid[J] <- W.grid[J] + w
      X.grid[J] <- X.grid[J] + w*x[i]
      Y.grid[J] <- Y.grid[J] + w*y[i]
    }
  }
  
  # normalize distributed information
  for(i in 1:n)
  {
    if(W.grid[i]>0)
    {
      X.grid[i] <- X.grid[i]/W.grid[i]
      Y.grid[i] <- Y.grid[i]/W.grid[i]
    }
  }
  # continuous weights eff up the FFT numerics so discretize weights
  W <- sum(W) # now total DOF
  W.grid <- sign(W.grid) # discrete weights
  
  return(list(w=W.grid,x=X.grid,y=Y.grid,lag=lag,dt=dt))
}


############################
# FFT VARIOGRAM
variogram.fast <- function(data,dt=NULL,CI="Markov")
{
  t <- data$t
  x <- data$x
  y <- data$y
  
  # smear the data over an evenly spaced time grid
  GRID <- gridder(t,x,y,dt)
  W.grid <- GRID$w
  X.grid <- GRID$x
  Y.grid <- GRID$y
  lag <- GRID$lag
  
  n <- length(lag)
  
  W.grid <- Conj(FFT(pad(W.grid,2*n)))
  XX.grid <- FFT(pad(X.grid^2,2*n))
  YY.grid <- FFT(pad(Y.grid^2,2*n))
  X.grid <- FFT(pad(X.grid,2*n))
  Y.grid <- FFT(pad(Y.grid,2*n))

  # pair number. one for x and y data
  DOF <- round(Re(2*IFFT(abs(W.grid)^2)[1:n]))
  # SVF un-normalized
  SVF <- Re(IFFT(Re(W.grid*(XX.grid+YY.grid))-(abs(X.grid)^2+abs(Y.grid)^2))[1:n])
  
  # delete missing lags
  SVF <- data.frame(SVF=SVF,DOF=DOF,lag=lag)
  SVF <- subset(SVF,DOF>0)
  lag <- SVF$lag
  DOF <- SVF$DOF
  SVF <- SVF$SVF
  
  # normalize SVF
  SVF <- SVF/DOF
  
  # only count non-overlapping lags... not perfect
  if(CI=="Markov")
  {
    dof <- 2*(last(t)-t[1])/lag
    dof[1] <- 2*length(t)
  
    for(i in 1:length(lag))
    {
      if(dof[i]<DOF[i]) {DOF[i] <- dof[i] }
    }
  }
  else if(CI=="IID") # fix initial and total DOF
  {
    DOF[1] <- 2*length(t)
    DOF[-1] <- DOF[-1]/sum(DOF[-1])*(length(t)^2-length(t))
  }
  
  result <- data.frame(SVF=SVF,DOF=DOF,lag=lag)
  result <- new.variogram(result, info=attr(data,"info"))
  
  return(result)
}

##################################
# LAG-WEIGHTED VARIOGRAM
variogram.slow <- function(data,dt=NULL,CI="Markov")
{
  t <- data$t
  x <- data$x
  y <- data$y

  n <- length(t)
  
  # time lags
  DT <- diff(t)
  DT.L <- c(DT[1],DT)
  DT.R <- c(DT,DT[n-1])
  
  # default time step
  if(is.null(dt))
  { dt <- stats::median(DT) }

  # where we will store stuff
  lag <- seq(0,ceiling((t[n]-t[1])/dt))*dt
  SVF <- rep(0,length(lag))
  DOF <- rep(0,length(lag))
  DOF2 <- rep(0,length(lag))
  
  pb <- utils::txtProgressBar(style=3)
  for(i in 1:n)
  { 
    for(j in i:n)
    {
      tau <- t[j] - t[i]
      var <- ((x[j]-x[i])^2 + (y[j]-y[i])^2)/4
      
      # gap weight
      if(tau==0) { w <- 1 }
      else { w <- (clamp(DT.L[j]/tau)+clamp(DT.R[j]/tau))*(clamp(DT.L[i]/tau)+clamp(DT.R[i]/tau)) }
      
      # fractional index
      k <- tau/dt + 1
      
      if(floor(k)==ceiling(k))
      { # even sampling
        # lag index
        K <- round(k)
        # total weight
        W <- w
        # accumulate
        SVF[K] <- SVF[K] + W*var
        DOF[K] <- DOF[K] + W
        DOF2[K] <- DOF2[K] + W^2
      }
      else
      { # account for drift by distributing semi-variance
        
        # left index
        K <- floor(k)
        # left weight
        W <- w*(1-(k-K))
        # accumulate left portion
        SVF[K] <- SVF[K] + W*var
        DOF[K] <- DOF[K] + W
        DOF2[K] <- DOF2[K] + W^2
        
        # right index
        K <- ceiling(k)
        # right weight
        W <- w*(1-(K-k))
        # accumulate right portion
        SVF[K] <- SVF[K] + W*var
        DOF[K] <- DOF[K] + W
        DOF2[K] <- DOF2[K] + W^2
      }
    }
    utils::setTxtProgressBar(pb,(i*(2*n-i))/(n^2))
  }
  
  # delete missing lags
  SVF <- data.frame(SVF=SVF,DOF=DOF,DOF2=DOF2,lag=lag)
  SVF <- subset(SVF,DOF>0)
  lag <- SVF$lag
  DOF <- SVF$DOF
  DOF2 <- SVF$DOF2
  SVF <- SVF$SVF
  
  # normalize SVF
  SVF <- SVF/DOF
  # effective DOF from weights, one for x and y
  DOF <- 2*DOF^2/DOF2
  
  # only count non-overlapping lags... still not perfect
  if(CI=="Markov")
  {
    dof <- 2*length(t)
    if(dof<DOF[1]) { DOF[1] <- dof  }
    
    for(i in 2:length(lag))
    { # large gaps are missing data
      dof <- 2*sum(DT[DT<=lag[i]])/lag[i]
      if(dof<DOF[i]) { DOF[i] <- dof }
      
      utils::setTxtProgressBar(pb,i/length(lag))
    }
  }
  else if(CI=="IID") # fix initial and total DOF
  {
    DOF[1] <- 2*length(t)
    DOF[-1] <- DOF[-1]/sum(DOF[-1])*(length(t)^2-length(t))
  }
  
  close(pb)
  
  result <- data.frame(SVF=SVF,DOF=DOF,lag=lag)
  result <- new.variogram(result, info=attr(data,"info"))
    
  return(result)
}


##########
plot.svf <- function(lag,model,alpha=0.05,col="red",type="l",...)
{
  # pull out relevant model parameters
  tau <- model$tau
  sigma <- model$sigma
  COV <- model$COV.tau
  
  if(length(tau)>0 && tau[1]==Inf){ range <- FALSE }else{ range <- TRUE }
  tau <- tau[tau>0]
  tau <- tau[tau<Inf]
  K <- length(tau)
  
  # trace variance
  GM.sigma <- sqrt(det(sigma))
  sigma <- mean(diag(sigma)) # now AM.sigma
  
  # parameter covariances
  # default to no error considered
  if(is.null(COV)) { COV <- diag(0,K+1) }
  
  # transform parameterization
  if(K>0 && range)
  {
    f <- 1/tau
    grad <- diag(c(1,-f^2),nrow=K+1)
    COV <- grad %*% COV %*% t(grad)
  }
    
  if(K==0 && range) # Bivariate Gaussian
  { 
    svf <- function(t){ if(t==0) {0} else {sigma} }
    grad.svf <- function(t){ svf(t)/GM.sigma }
  }
  else if(K==0) # Brownian motion
  {
    svf <- function(t){ sigma*t }
    grad.svf <- function(t){ svf(t)/GM.sigma }
  }
  else if(K==1 && range) # OU motion
  {
    svf <- function(t){ sigma*(1-exp(-t*f)) }
    grad.svf <- function(t){ c( svf(t)/GM.sigma , sigma*(-t*exp(-t*f)) ) }
  }
  else if(K==1) # IOU motion
  {
    svf <- function(t) { sigma*(t-tau*(1-exp(-t/tau))) }
    grad.svf <- function(t){ c( svf(t)/GM.sigma , -sigma*(1-(1+t/tau)*exp(-t/tau)) ) }
  }
  else if(K==2) # OUF motion
  { 
    svf <- function(t){ sigma*(1 - diff(tau*exp(-t/tau))/diff(tau)) } 
    grad.svf <- function(t)
    { c(
      svf(t)/GM.sigma ,
      -sigma*f[2]/diff(f)^2*((1+diff(f)*t)*exp(-f[1]*t)-exp(-f[2]*t)) ,
      -sigma*f[1]/diff(f)^2*((1-diff(f)*t)*exp(-f[2]*t)-exp(-f[1]*t))
    ) }
  }

  SVF <- Vectorize(function(t){svf(t)})
  
  # point estimate plot
  graphics::curve(SVF,from=0,to=lag,n=1000,add=TRUE,col=col,type=type,...)
  
  # confidence intervals if COV provided
  if(any(diag(COV)>0))
  {
    Lags <- seq(0,lag,lag/1000)
    
    # efffective DOF of SVF if its chi square
    DOF.svf <- function(t){ 2*svf(t)^2/( grad.svf(t) %*% COV %*% grad.svf(t) ) }
    
    for(j in 1:length(alpha))
    {
      svf.lower <- Vectorize(function(t){ svf(t) * CI.lower(DOF.svf(t),alpha[j]) })
      svf.upper <- Vectorize(function(t){ svf(t) * CI.upper(DOF.svf(t),alpha[j]) })
      
      graphics::polygon(c(Lags,rev(Lags)),c(svf.lower(Lags),rev(svf.upper(Lags))),col=translucent(col,alpha=0.1/length(alpha)),border=NA,...)
    }
  }
  
}

###########################################################
# PLOT VARIOGRAM
###########################################################
plot.variogram <- function(x, CTMM=NULL, alpha=0.05, fraction=0.5, col="black", col.CTMM="red", ...)
{  
  # number of variograms
  if(class(x)=="variogram" || class(x)=="data.frame") { x <- list(x) }
  n <- length(x)
  
  # maximum lag in data
  max.lag <- sapply(x, function(v){ last(v$lag) } )
  max.lag <- max(max.lag)
  # subset fraction of data
  max.lag <- fraction*max.lag
  
  # subset all data to fraction of total period
  x <- lapply(x, function(v) { subset.data.frame(v, lag <= max.lag) })

  # maximum CI on SVF
  max.SVF <- max(sapply(x, function(v){ max(v$SVF * CI.upper(v$DOF,min(alpha))) } ))
  # limit plot range to twice max SVF point estimate (otherwise hard to see)
  max.cap <- 2*max(sapply(x, function(v){ max(v$SVF) } ))
  if(max.SVF>max.cap) { max.SVF <- max.cap }
  
  # choose SVF units
  SVF.scale <- unit(max.SVF,"area")
  SVF.name <- SVF.scale$name
  SVF.scale <- SVF.scale$scale
  
  # choose lag units
  lag.scale <- unit(max.lag,"time",2)
  lag.name <- lag.scale$name
  lag.scale <- lag.scale$scale
  
  xlab <- paste("Time-lag ", "(", lag.name, ")", sep="")
  ylab <- paste("Semi-variance ", "(", SVF.name, ")", sep="")
  
  # fix base plot layer
  plot(c(0,max.lag/lag.scale),c(0,max.SVF/SVF.scale), xlab=xlab, ylab=ylab, col=grDevices::rgb(1,1,1,0), ...)
  
  # color array for plots
  col <- array(col,n)
  
  for(i in 1:n)
  {
    lag <- x[[i]]$lag/lag.scale
    SVF <- x[[i]]$SVF/SVF.scale
    DOF <- x[[i]]$DOF
        
    # make sure plot looks nice and appropriate for data resolution
    type <- "l"
    if(length(lag) < 100) { type <- "p" }
    
    graphics::points(lag, SVF, type=type, col=col[[i]])
    
    for(j in 1:length(alpha))
    {
      SVF.lower <- SVF * CI.lower(DOF,alpha[j])
      SVF.upper <- SVF * CI.upper(DOF,alpha[j])
      
      graphics::polygon(c(lag,rev(lag)),c(SVF.lower,rev(SVF.upper)),col=translucent(col[[i]],alpha=0.1),border=NA)
    }
  }
  
  # NOW PLOT THE MODELS
  if(!is.null(CTMM))
  {
    if(class(CTMM)=="ctmm") { CTMM <- list(CTMM) }
    n <- length(CTMM) 
    
    # color array for plots
    col <- array(col.CTMM,n)
    type <- "l"
    
    for(i in 1:n)
    {
      # units conversion
      CTMM[[i]]$sigma <- CTMM[[i]]$sigma/SVF.scale
      if(length(CTMM[[i]]$tau)>0){ CTMM[[i]]$tau <- CTMM[[i]]$tau/lag.scale }
      
      K <- length(CTMM[[i]]$COV.tau[1,])
      scale <- SVF.scale
      if(K>1){ scale <- c(scale,rep(lag.scale,K-1)) }
      
      # variance -> diffusion adjustment
      if(length(CTMM[[i]]$tau)>0 && max(CTMM[[i]]$tau)==Inf)
      {
        CTMM[[i]]$sigma <- CTMM[[i]]$sigma*lag.scale
        scale[1] <- scale[1]/lag.scale
      }
      
      if(K>0)
      {
        scale <- diag(1/scale,length(scale))
        CTMM[[i]]$COV.tau <- scale %*% CTMM[[i]]$COV.tau %*% scale
      }
      
      plot.svf(max.lag/lag.scale,CTMM[[i]],alpha=alpha,type=type,col=col[[i]])
    }
  }
  
}
# PLOT.VARIOGRAM METHODS
#methods::setMethod("plot",signature(x="variogram",y="missing"), function(x,y,...) plot.variogram(x,...))
#methods::setMethod("plot",signature(x="variogram",y="variogram"), function(x,y,...) plot.variogram(list(x,y),...))
#methods::setMethod("plot",signature(x="variogram",y="ctmm"), function(x,y,...) plot.variogram(x,model=y,...))
#methods::setMethod("plot",signature(x="variogram"), function(x,...) plot.variogram(x,...))


#######################################
# plot a variogram with zoom slider
#######################################
zoom.variogram <- function(x, fraction=0.5, ...)
{
  # R CHECK CRAN BUG WORKAROUND
  z <- NULL
  
  # number of variograms
  n <- 1
  if(class(x)=="list") { n <- length(x) }
  else {x <- list(x) } # promote to list of one
  
  # maximum lag in data
  max.lag <- sapply(x, function(v){ last(v$lag) } )
  max.lag <- max(max.lag)
  
  min.lag <- sapply(x, function(v){ v$lag[2] } )
  min.lag <- min(min.lag)
  
  b <- 4
  min.step <- 10*min.lag/max.lag
  manipulate::manipulate( { plot.variogram(x, fraction=b^(z-1), ...) }, z=manipulate::slider(1+log(min.step,b),1,initial=1+log(fraction,b),label="zoom") )
}
#methods::setMethod("zoom",signature(x="variogram",y="missing"), function(x,y,...) zoom.variogram(x,...))
#methods::setMethod("zoom",signature(x="variogram",y="variogram"), function(x,y,...) zoom.variogram(list(x,y),...))
#methods::setMethod("zoom",signature(x="variogram",y="ctmm"), function(x,y,...) zoom.variogram(x,model=y,...))
#methods::setMethod("zoom",signature(x="variogram"), function(x,...) zoom.variogram(x,...))


####################################
# guess variogram model parameters #
####################################
variogram.guess <- function(variogram,range=TRUE)
{
  n <- length(variogram$lag)
  
  sigma <- mean(variogram$SVF[2:n])
  
  # peak diffusion rate estimate
  D <- max((variogram$SVF/variogram$lag)[2:n])
  
  # peak curvature estimate
  v2 <- 2*max((variogram$SVF/variogram$lag^2)[2:n])
  
  # position, velocity timescale estimate
  tau <- c(sigma/D,D/v2)
  
  if(!range) { sigma <- D ; tau[1] <- Inf }
  
  model <- ctmm(sigma=sigma,tau=tau,info=attr(variogram,"info"))
  return(model)
}


######################################################################
# visual fit of the variogram
######################################################################
variogram.fit <- function(variogram,range=TRUE,CTMM=NULL,name="variogram.fit.model",...)
{
  # R CHECK CRAN BUG WORKAROUND
  tau.v <- NULL
  store <- NULL ; rm(store)
  z <- NULL
  
  m <- 3 # slider length relative to point guestimate
  n <- length(variogram$lag)
  
  if(is.null(CTMM)) { CTMM <- variogram.guess(variogram,range=range) }
  
  sigma <- mean(diag(CTMM$sigma))
  tau <- CTMM$tau
  
  # slider maxima
  sigma.max <- m*sigma
  tau.max <- m*tau
  
  # parameters for logarithmic slider
  b <- 4
  min.step <- 10*variogram$lag[2]/variogram$lag[n]
  
  # manipulation controls
  manlist <- list(z = manipulate::slider(1+log(min.step,b),1,initial=1+log(0.5,b),label="zoom"),
                  sigma = manipulate::slider(0,sigma.max,initial=sigma,label="sigma variance (m^2)"),
                  tau.r = manipulate::slider(0,tau.max[1],initial=tau[1],label="tau position (s)"),
                  tau.v = manipulate::slider(0,tau.max[2],initial=tau[2],label="tau velocity (s)"),
                  store = manipulate::button(paste("save to",name))
  )
  
  if(!range)
  {
    manlist[[2]] <- manipulate::slider(0,sigma.max,initial=sigma,label="sigma diffusion (m^2/s)")
    manlist <- manlist[-3]
    tau.r <- Inf
  }
  
  manipulate::manipulate(
    {
      CTMM <- ctmm(sigma=sigma,tau=c(tau.r,tau.v),info=attr(variogram,"info"))
      if(store) { eval(parse(text=paste(name," <<- CTMM"))) } 
      plot.variogram(variogram,CTMM=CTMM,fraction=b^(z-1),...)
    },
    manlist
  )
}

# AVERAGE VARIOGRAMS
mean.variogram <- function(x,...)
{
  x <- c(x,list(...))
  IDS <- length(x)
  
  # assemble observed lag range
  lag.min <- rep(0,IDS) # this will be dt
  lag.max <- rep(0,IDS)
  for(id in 1:IDS)
  { 
    n <- length(x[[id]]$lag)
    lag.max[id] <- x[[id]]$lag[n]
    lag.min[id] <- x[[id]]$lag[2]
  }
  lag.max <- max(lag.max)
  lag.min <- min(lag.min)
  lag <- seq(0,ceiling(lag.max/lag.min))*lag.min
  
  # where we will store everything
  n <- length(lag)
  SVF <- rep(0,n)
  DOF <- rep(0,n)
  
  # accumulate semivariance
  for(id in 1:IDS)
  {
   for(i in 1:length(x[[id]]$lag))
   {
    # lag index
    j <- 1 + round(x[[id]]$lag[i]/lag.min)
    # number weighted accumulation
    DOF[j] <- DOF[j] + x[[id]]$DOF[i]
    SVF[j] <- SVF[j] + x[[id]]$DOF[i]*x[[id]]$SVF[i]
   }
  }
  
  # delete missing lags
  variogram <- data.frame(SVF=SVF,DOF=DOF,lag=lag)
  variogram <- subset(variogram,DOF>0)

  # normalize SVF
  variogram$SVF <- variogram$SVF / variogram$DOF
  
  # drop unused levels
  variogram <- droplevels(variogram)
  
  variogram <- new.variogram(variogram, info=mean.info(x))

  return(variogram)
}
#methods::setMethod("mean",signature(x="variogram"), function(x,...) mean.variogram(x,...))


# consolodate info attributes from multiple datasets
mean.info <- function(x)
{
  # mean identity
  identity <- sapply(x , function(v) { attr(v,"info")$identity } )
  identity <- unique(identity) # why did I do this?
  identity <- paste(identity,collapse=" ")
  
  info=attr(x[[1]],"info")
  info$identity <- identity
  
  return(info)
  }
