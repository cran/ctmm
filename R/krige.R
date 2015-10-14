################################
# Return hidden state estimates
# more appropriate to call this "smooth", but maybe more confusing
################################
filter <- function(data,model)
{
  return(kalman(data,model,smooth=TRUE))
}

##################################
# Plot the state estimate as if it were telemetry data
# This needs to be listified, ...
##################################
plot.state <- function(state,...)
{
  info <- state$model$info
  info$type <- "Krige"
  
  state <- data.frame(t=state$t,x=state$H[,1,1],y=state$H[,1,2],vx=state$H[,2,1],vy=state$H[,2,2])
  state <- new.telemetry(state,info=info)
  
  plot.telemetry(state,...)
}

#################################
# Kriged Kernel Density Estimate
# H is your additional smoothing bandwidth matrix (zero by default)
# resolution is the number of kriged locations per median step
# cor.min is roughly the correlation required between locations to bridge them
# dt.max is (alternatively) the maximum gap allowed between locations to bridge them
#################################
kkde <- function(data,model,H=diag(0,2),res.time=20,res.space=1000,cor.min=0.5,dt.max=NULL)
{
  model <- ctmm.prepare(data,model)
  
  # attach zero error to observations... VERY TEMPORARY
  data$error <- rep(0,dim(data)[1])
  
  # FIX THE TIME GRID TO AVOID SMALL DT
  
  # target resolution
  t <- data$t
  DT <- diff(t)
  dt <- stats::median(DT)/res.time
  
  # maximum gap to bridge
  if(is.null(dt.max)) { dt.max <- -log(cor.min)*model$tau[1] }
  
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
  
  state <- filter(data,model)
  
  # evenly sampled subset: data points (bridge ends) may be counted twice and weighted half
  GRID <- c( where(data$t %in% t.grid) , where(data$t %in% t.grid[where(diff(t.grid)==0)]) )
  GRID <- sort.int(GRID,method="quick")

#  GRID <- data$t %in% t.grid
  
#  t <- state$t[GRID]
  x <- state$H[GRID,"position","x"]
  y <- state$H[GRID,"position","y"]
  C <- state$C[GRID,"position","position"]
  
  # continuous velocities will give us more information to use
  if(length(model$tau)>1)
  {
    v.x <- state$H[GRID,"velocity","x"]
    v.y <- state$H[GRID,"velocity","y"]
  }
  else # empty velocity data otherwise
  {
    v.x <- rep(0,length(x))
    v.y <- rep(0,length(y))
  }

  # fake data
  data <- data.frame(x=x,y=y)

  # some variances are slightly negative due to roundoff error
  C <- clamp(C,min=0,max=Inf)
  
  n <- length(x)
  # uncertainties/bandwidths for this data
  h <- H # extra smoothing
  H <- array(0,c(n,2,2))
  for(i in 1:n)
  {
    v <- c(v.x[i],v.y[i])
    # total covariance is bandwidth + krige uncertainty + kinetic numerical correction
    H[i,,] <- h + model$sigma*C[i] + dt.grid[i]^2/12*(v %o% v)
    # maybe publish the last part as a technical note
  }
  
  # need to fix and understand this alpha problem
  kde(data,H=H,W=w.grid,res=res.space)
}

# here is some sample code for this function
# KDE <- ctmm:::kkde(DATA[2:100,],MODEL,H=diag(5^2,2),res.time=100)
# image(x=KDE$x/1000,y=KDE$y/1000,z=100*KDE$CDF,col=gray(0:10^5/10^5))
# contour(x=KDE$x/1000,y=KDE$y/1000,z=100*KDE$CDF,levels=c(95,99))


################################
# TOTALLY UNFINISHED
# Krige at one time or over list of times
krige <- function(Krige, t)
{
  # sampled times
  T <- Krige$t
  N <- length(T)
  dt <- Krige$dt

  # Kriged time indices (bisection lookup)
  n <- length(t)
  I <- array(1,n)
  for(i in 1:n) 
  {
    I.min <- I[i]
    I.max <- N
    
    # iterate two methods:
    # biased/informed bisection that is exact for evenly sampled data
    # naive bisection that always converges
    while(I.max-I.min>1)
    {
      # VANILLA BISECTION
      dI <- floor((I.max-I.min)/2)
      I.test <- I.min + dI
      if(T[I.test]<=t) { I.min <- I.test }
      
      I.test <- I.max - dI
      if(T[I.test]<=t) { I.min <- I.test }
      
      # BIASED/INFORMED BISECTION  
      # How many steps we need to push forward for evenly sampled data
      I.test <- I.min + floor((t[i]-T[I.min])/dt)
      if(T[I.test]<=t) { I.min <- I.test }
      
      # How many steps we need to push backwards for evenly sampled data
      I.test <- I.max - floor((T[I.max]-t[i])/dt)
      if(T[I.test]>=t) { I.max <- I.test }
    }

    # Update max-min
    I[i] <- I.min
    if(i<n) { I[i+1] <- I[i] }
  }
  
  
  
}


