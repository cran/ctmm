################################
# Return hidden state estimates
################################
smoother <- function(DATA,CTMM)
{
  CTMM <- ctmm.prepare(DATA,CTMM)
  
  x <- DATA$x
  y <- DATA$y
  t <- DATA$t
  
  dt <- c(Inf,diff(t))
  n <- length(t)
  
  sigma <- CTMM$sigma
  K <- length(CTMM$tau)
  
  if(!CTMM$error)
  {
    CTMM$sigma <- 1
    KALMAN <- kalman(cbind(x,y),dt=dt,CTMM=CTMM,error=DATA$error,smooth=TRUE)
    
    R <- KALMAN$Z[,"position",]
    if(K>1) { V <- KALMAN$Z[,"velocity",] }
    COV <- KALMAN$S[,"position","position"] %o% sigma
  }
  else if(CTMM$isotropic)
  {
    CTMM$sigma <- sigma@par[1]
    KALMAN <- kalman(cbind(x,y),dt=dt,CTMM=CTMM,error=DATA$error,smooth=TRUE)
    
    R <- KALMAN$Z[,"position",]
    if(K>1) { V <- KALMAN$Z[,"velocity",] }
    COV <- KALMAN$S[,"position","position"] %o% diag(1,2)
  }
  else
  {
    #diagonalize data and then run two 1D Kalman filters with separate means
    R <- cbind(x,y)
    V <- array(0,dim=c(n,2))
    COV <- array(0,dim=c(n,2,2))
    
    GM <- sigma@par[1]
    e <- sigma@par[2]
    theta <- sigma@par[3]
    
    M <- rotate(-theta)
    R <- t(M %*% t(R))
    
    # major axis likelihood
    CTMM$sigma <- GM * exp(+e/2)
    KALMAN <- kalman(cbind(R[,1]),dt=dt,CTMM=CTMM,error=DATA$error,smooth=TRUE)
    
    R[,1] <- KALMAN$Z[,"position",1]
    if(K>1) { V[,1] <- KALMAN$Z[,"velocity",1] }
    COV[,1,1] <- KALMAN$S[,"position","position"]

    # minor axis likelihood
    CTMM$sigma <- GM * exp(-e/2)
    KALMAN <- kalman(cbind(R[,2]),dt=dt,CTMM=CTMM,error=DATA$error,smooth=TRUE)
    
    R[,2] <- KALMAN$Z[,"position",1]
    if(K>1) { V[,2] <- KALMAN$Z[,"velocity",1] }
    COV[,2,2] <- KALMAN$S[,"position","position"]
    
    # transform results back
    M <- rotate(+theta)
    
    # there MUST be an easier way to do a simple inner product but %*% fails
    COV <- aperm(COV,perm=c(2,1,3))
    dim(COV) <- c(2,2*n)
    COV <- M %*% COV
    dim(COV) <- c(2*n,2)
    COV <- COV %*% t(M)
    dim(COV) <- c(2,n,2)
    COV <- aperm(COV,perm=c(2,1,3))
    
    R <- t(M %*% t(R))
    V <- t(M %*% t(V))
  }
  
  CTMM$sigma <- sigma
  CTMM <- ctmm.repair(CTMM)
  
  RETURN <- list(t=t,R=R,COV=COV)
  if(K>1) { RETURN$V <- V }
  
  return(RETURN)
}

#################################
# Kriged Kernel Density Estimate
# H is your additional smoothing bandwidth matrix (zero by default)
# resolution is the number of kriged locations per median step
# cor.min is roughly the correlation required between locations to bridge them
# dt.max is (alternatively) the maximum gap allowed between locations to bridge them
#################################
occurrence <- function(data,CTMM,H=diag(0,2),res.time=20,res.space=1000,grid=NULL,cor.min=0.5,dt.max=NULL)
{
  info <- attr(data,"info")
  CTMM <- ctmm.prepare(data,CTMM)
  
  # prepare data error
  data$error <- error.prepare(data,CTMM)
  
  # FIX THE TIME GRID TO AVOID TINY DT
  # target resolution
  t <- data$t
  DT <- diff(t)
  dt <- stats::median(DT)/res.time
  
  # maximum gap to bridge
  if(is.null(dt.max)) { dt.max <- -log(cor.min)*CTMM$tau[1] }
  
  # gaps to bridge
  BRIDGE <- where(DT<=dt.max)
  # this regularization is not perfectly regular, but holds up to sampling drift in caribou data
  t.grid <- c()  # full (even-ish) grid
  dt.grid <- c() # local (numeric) sampling resolution 
  t.new <- c()   # new times in this even grid
  for(i in BRIDGE)
  {
    n.sub <- round(DT[i]/dt)+1
    t.sub <- seq(from=t[i],to=t[i+1],length.out=n.sub)
    
    dt.sub <- DT[i]/(n.sub-1)
    dt.sub <- rep(dt.sub,n.sub)
    
    t.grid <- c(t.grid,t.sub)
    dt.grid <- c(dt.grid,dt.sub)

    t.new <- c(t.new,t.sub[c(-1,-n.sub)])
  }

  # half weight repeated endpoints in grid
  w.grid <- dt.grid
  REPEAT <- where(diff(t.grid)==0)
  w.grid[REPEAT] <- w.grid[REPEAT]/2
  w.grid[REPEAT+1] <- w.grid[REPEAT+1]/2
  
#   # toss out times in big gaps
#   TOSS <- DT > dt.max
#   TOSS <- (1:length(TOSS))[TOSS]
#   # start and stop indices to avoid big gaps
#   START <- c(1,TOSS+1)
#   STOP <- c(TOSS,length(t))
#   # fine grid skipping big gaps  
#   t.grid <- sapply(1:length(START),function(i){ seq(t[START[i]],t[STOP[i]],dt) })
#   t.grid <- unlist(t.grid)
# 
#   # times not sampled in data
#   t.new <- t.grid[! (t.grid %in% t)]

  # empty observation row for these times
  blank <- data[1,]
  blank$error <- Inf
  blank <- blank[rep(1,length(t.new)),]
  blank$t <- t.new
  
  # attach empty measurements to data
  data <- rbind(data,blank)
  
  # sort times
  data <- data[sort.list(data$t,na.last=NA,method="quick"),]
  # this is now our fake data set to feed into the kalman smoother
  
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # FIX EVERYTING BELOW THIS LINE
  # NEED TO EITHER RUN ONE FILTER OR TWO FILTERS WITH A ROTATION
  state <- smoother(data,CTMM)
  
  # evenly sampled subset: data points (bridge ends) may be counted twice and weighted half
  GRID <- c( where(data$t %in% t.grid) , where(data$t %in% t.grid[where(diff(t.grid)==0)]) )
  GRID <- sort.int(GRID,method="quick")

#  GRID <- data$t %in% t.grid
  
#  t <- state$t[GRID]
  R <- state$R[GRID,]
  COV <- state$COV[GRID,,]
  n <- length(R[,1])
  
  # continuous velocities will give us more information to use
  if(length(CTMM$tau)>1)
  { V <- state$V[GRID,] }
  else # null velocity data otherwise
  { V <- array(0,c(n,2)) }

  # fake data
  data <- data.frame(x=R[,1],y=R[,2])

  # some covariances are slightly negative due to roundoff error
  # so here I clamp the negative singular values to zero
  COV <- lapply(1:n,function(i){ PDpart(COV[i,,]) })
  COV <- unlist(COV)
  dim(COV) <- c(2,2,n)
  COV <- aperm(COV,perm=c(3,1,2))
  
  # uncertainties/bandwidths for this data
  h <- H # extra smoothing
  H <- array(0,c(n,2,2))
  for(i in 1:n)
  {
    # total covariance is bandwidth + krige uncertainty + kinetic numerical correction
    H[i,,] <- h + COV[i,,] + dt.grid[i]^2/12*(V[i,] %o% V[i,])
    # maybe publish the last part as a technical note
  }
  # there is probably a faster way to do that
  
  # using the same data format as AKDE, but with only the ML estimate (alpha=1)
  KDE <- kde(data,H=H,W=w.grid,res=res.space)
  KDE$H <- diag(0,2)
  KDE <- list(ML=KDE)
  KDE <- new.UD(KDE,info=info,level=0)
  return(KDE)
}

# here is some sample code for this function
# KDE <- ctmm:::kkde(DATA[2:100,],MODEL,H=diag(5^2,2),res.time=100)
# image(x=KDE$x/1000,y=KDE$y/1000,z=100*KDE$CDF,col=gray(0:10^5/10^5))
# contour(x=KDE$x/1000,y=KDE$y/1000,z=100*KDE$CDF,levels=c(95,99))