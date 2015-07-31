# akde object generator
# list of kde objects with info slots
new.akde <- methods::setClass("akde", representation(info="list",alpha="numeric"), contains="list")


# Slow lag counter
tab.lag.DOF <- function(data,fast=NULL,dt=NULL)
{
  t <- data$t
  # intelligently select algorithm
  n <- length(t)
  if(is.null(fast))
  {
    if(n<100) { fast <- FALSE }
    else { fast <- TRUE }
  }
  
  # calculate lag,n(lag) vectors
  if(fast)
  {
    lag.DOF <- variogram.fast(data,dt=dt,CI="IID") # count all lag pairs
    lag.DOF$SVF <- NULL
    
    # lag==0 should not be doubled
    lag.DOF$DOF[1] <- lag.DOF$DOF[1]/2 
  }
  else
  {
    lag <- outer(t,t,FUN="-")
    lag <- abs(lag) # may have to change far in the future
    lag <- c(lag) # collapse to 1D array
    
    # can we think of a faster sorter?
    r <- lag # to be the remainder
    lag <- NULL
    DOF <- NULL
    while(length(r)>0)
    {
      n <- length(r)
      lag <- c(lag,r[1])
      r <- r[r!=r[1]]
      DOF <- c(DOF,n-length(r))
    }
    
    lag.DOF <- list(lag=lag,DOF=DOF)
  }
  
  return(lag.DOF)
}


##################################
# Bandwidth optimizer
#lag.DOF is an unsupported option for end users
akde.bandwidth <- function(data,CTMM,fast=NULL,dt=NULL)
{
  lag.DOF <- tab.lag.DOF(data,fast=fast,dt=dt)
  # Extract lag data
  DOF <- lag.DOF$DOF
  lag <- lag.DOF$lag
    
  tau <- CTMM$tau
  tau <- tau[tau>0]
  K <- length(tau)
  
  # standardized SVF
  if(K==0)
  {
    g <- Vectorize( function(t) { if(t==0) {0} else {1} } )
  }
  else if(K==1)
  {
    g <- Vectorize( function(t) { 1-exp(-t/tau) } )
  }
  else if(K==2)
  {
    g <- Vectorize( function(t) { 1-diff(tau*exp(-t/tau))/diff(tau) } )
  }
  
  n <- length(data$t)
  
  # Mean Integrated Square Error modulo a constant
  MISE <- function(h)
  {
    if(h<=0) {Inf}
    else { (1/n^2)*sum(DOF/(2*g(lag)+2*h^2)) - 2/(2+h^2) + 1/2 }
  }
  
  h <- 1/n^(1/6) # User Silverman's rule of thumb to place lower bound
  h <- stats::optimize(f=MISE,interval=c(h/2,2))$minimum
  
  H <- h^2*CTMM$sigma
  
  rownames(H) <- c("x","y")
  colnames(H) <- c("x","y")
  
  return(H)
}

#######################################
# wrap the kde function for our telemetry data format and CIs.
akde <- function(data,CTMM,alpha=0.05,fast=NULL,dt=NULL,error=0.005)
{
  pb <- utils::txtProgressBar(min=-1,max=3,initial=-1,style=3)
  
  tau <- CTMM$tau
  K <- length(tau)
  
  sigma <- CTMM$sigma
  # covariance area/scale
  GM.sigma <- sqrt(det(sigma))
  # orientation matrix
  R <- sigma/GM.sigma

  # wrapper function to estimate bandwidth matrix
  # par = c(sigma.GM,tau)
  fn.H <- function(par)
  {
    CTMM.par <- CTMM
    CTMM.par$sigma <- par[1]*R
    if(K>0) { CTMM.par$tau <- par[-1] }
    H <- akde.bandwidth(data=data,CTMM=CTMM.par,fast=fast,dt=dt)
    return(H)
  }

  # wrapper function to estimate bandwidth area
  # par = c(sigma.GM,tau)
  fn.GM.H <- function(par)
  {
    H <- fn.H(par)
    GM.H <- sqrt(det(H))
    return(GM.H)
  }
  # How strongly optimal bandwidth varries with parameter estimates
  d.GM.H <- numDeriv::grad(fn.GM.H,c(GM.sigma,tau))
  
  # ML propagated curvature covariance
  COV <- d.GM.H %*% CTMM$COV.tau %*% d.GM.H
  # ML Bandwidth area
  GM.H <- fn.GM.H(c(GM.sigma,tau))
  # confidence intervals from chi^2
  GM.H <- chisq.ci(GM.H,COV,alpha)

  # data formatted for ks::kde
  x <- cbind(data$x,data$y)
  
  # object to store crap
  KDE <- list(low=0,ML=0,high=0)
  
  utils::setTxtProgressBar(pb,0)
  for(i in 1:3)
  {
    H <- GM.H[i]*R
    KDE[[i]] <- kde(data,H,alpha=error)
    utils::setTxtProgressBar(pb,i)
  }
  
  KDE <- new.akde(KDE,info=attr(data,"info"),alpha=alpha)
  
  close(pb)
  return(KDE)
}


##################################
# construct my own kde objects
# was using ks-package but it has some bugs
# alpha is the error goal in my total probability
kde <- function(data,H,alpha=0.005)
{
  x <- data$x
  y <- data$y
  
  min.x <- min(x)
  max.x <- max(x)
  
  min.y <- min(y)
  max.y <- max(y) 
  
  std.x <- sqrt(H[1,1])
  std.y <- sqrt(H[2,2])
  
  # how far to extend range from data as to ensure alpha significance in total probability
  z <- sqrt(-2*log(alpha))
  
  # kernel buffer sizes
  DX <- z*std.x
  DY <- z*std.y
  
  min.x <- min.x - DX
  max.x <- max.x + DX
  
  min.y <- min.y - DY
  max.y <- max.y + DY
  
  # how fine to make grid as to ensure alpha significance in total probability
  # ~0.4 calculated from asymptotic behavior of elliptic theta function from Riemmann sum approximation of Gaussian integral
  
  z <- (alpha/0.4)
  
  dx <- z*std.x
  dy <- z*std.y
  
  mu.x <- (min.x+max.x)/2
  mu.y <- (min.y+max.y)/2
  
  n.x <- ceiling((max.x-min.x)/2/dx)
  n.y <- ceiling((max.y-min.y)/2/dy)
  
  # grid locations
  X <- mu.x + ((-n.x):(n.x))*dx
  Y <- mu.y + ((-n.y):(n.y))*dy
  
  H.inv <- solve(H)

  # could maybe speed this up with compact subset of X,Y + translation
  Gauss <- function(x,y,r1,r2,c1,c2)
  {
    x <- x-X[r1:r2]
    y <- y-Y[c1:c2]
    # grid of -1/2 * standard variances without correlation
    R <- outer(-H.inv[1,1]/2*x^2,-H.inv[2,2]/2*y^2,"+")
    # now add correlation grid
    R <- R - H.inv[1,2]*(x %o% y)
    R <- exp(R)
    return(R)
  }
  
  pdf <- array(0,c(length(X),length(Y)))
  for(i in 1:length(data$x))
  {
    # calculate sub-grid indices
    r1 <- floor((x[i]-DX-X[1])/dx) + 1
    r2 <- ceiling((x[i]+DX-X[1])/dx) + 1

    c1 <- floor((y[i]-DY-Y[1])/dy) + 1
    c2 <- ceiling((y[i]+DY-Y[1])/dy) + 1
    
    pdf[r1:r2,c1:c2] <- pdf[r1:r2,c1:c2] + Gauss(x[i],y[i],r1,r2,c1,c2)
  }
  pdf <- pdf / (length(data$x) * 2*pi*sqrt(det(H)))

  result <- list(pdf=pdf,x=X,y=Y,H=H,dA=dx*dy)
  class(result) <- "kde"
  
  return(result)
}


##################
# find probability density of % contour
qkde <- function(kde,alpha=0.05,cells=FALSE)
{
  # here we make a lookup table to find the pdf @ a given cumulative probability from the max pdf down
  pdf <- c(kde$pdf)
  pdf <- sort(pdf,decreasing=TRUE,method="quick")
  cdf <- cumsum(pdf)*kde$dA
  
  # search sorted list
  n <- findInterval(1-alpha,cdf)
  p <- pdf[n]
  
  if(cells) { return(n) }
  else { return(p) }  
}


#######################
# summarize details of akde object
summary.akde <- function(object,alpha=0.05,...)
{
  area <- c(0,0,0)
  for(i in 1:3)
  {
    n <- qkde(object[[i]],alpha,cells=TRUE)
    area[i] <- n * object[[i]]$dA
  }
  
  unit.info <- unit(area,"area")
  name <- unit.info$name
  scale <- unit.info$scale
  
  area <- array(area/scale,c(1,3))
  colnames(area) <- c("low","ML","high")
  rownames(area) <- paste("area (",name,")",sep="")
  
  return(area)
}
#methods::setMethod("summary",signature(object="akde"), function(object,...) summary.akde(object,...))


################################
# create a raster of the ML akde
raster.akde <- function(AKDE,CI="ML")
{
  kde <- AKDE[[CI]]
  dx <- kde$x[2]-kde$x[1]
  dy <- kde$y[2]-kde$y[1]
  
  xmn <- kde$x[1]-dx/2
  xmx <- last(kde$x)+dx/2
  
  ymn <- kde$y[1]-dy/2
  ymx <- last(kde$y)+dy/2
  
  Raster <- raster::raster(t(kde$pdf[,dim(kde$pdf)[2]:1]),xmn=xmn,xmx=xmx,ymn=ymn,ymx=ymx,crs=attr(AKDE,"info")$projection)
  
  return(Raster)
}


################
# Is contour A inside contour B
inside <- function(A,B)
{
  result <- mode(sp::point.in.polygon(A$x,A$y,B$x,B$y))
  if(1<=result && result<=2) { return(1) } else { return(0) }
}


##############
SpatialPolygonsDataFrame.akde <- function(AKDE,alpha=0.05)
{
  ID <- paste(AKDE@info$identity," ",names(AKDE)," ",round(100*(1-alpha)),"%",sep="")

  polygons <- list()
  for(i in 1:length(AKDE))
  {
    kde <- AKDE[[i]]
    p <- qkde(kde,alpha=alpha)
    CL <- grDevices::contourLines(x=kde$x,y=kde$y,z=kde$pdf,levels=p)
    
    # create contour heirarchy matrix (half of it)
    H <- array(0,c(1,1)*length(CL))    
    for(row in 1:length(CL))
    {
      for(col in row:length(CL))
      {
        H[row,col] <- inside(CL[[row]],CL[[col]]) 
      }
    }
    
    # number of contours that this contour is inside
    I <- rowSums(H)
    
    # if I is odd, then you are a hole inside a positive area
    hole <- is.odd(I)
    
    # polygon
    polygons[[i]] <- list()
    for(j in 1:length(CL))
    {
      polygons[[i]][[j]] <- sp::Polygon( cbind( CL[[j]]$x , CL[[j]]$y ) , hole=hole[j] )
    }

    # polygonS
    polygons[[i]] <- sp::Polygons(polygons[[i]],ID=ID[i])
  }
  names(polygons) <- ID

    # spatial polygons
  polygons <- sp::SpatialPolygons(polygons, proj4string=sp::CRS(attr(AKDE,"info")$projection))

  # spatial polygons data frame  
  data <- data.frame(name=rev(ID))
  rownames(data) <- rev(ID)
  polygons <- sp::SpatialPolygonsDataFrame(polygons,data)
  
  return(polygons)
}


################
writeShapefile.akde <- function(AKDE, folder, file=AKDE@info$identity, alpha=0.05,  ...)
{
  SP <- SpatialPolygonsDataFrame.akde(AKDE,alpha=alpha)
  
  rgdal::writeOGR(SP, dsn=folder, layer=file, driver="ESRI Shapefile",...)
}

