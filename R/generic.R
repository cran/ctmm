# is.installed <- function(pkg) is.element(pkg, installed.packages()[,1]) 

# worst case FFT functions
FFT <- function(X) { stats::fft(X,inverse=FALSE) }
IFFT <- function(X) { stats::fft(X,inverse=TRUE)/length(X) }

# THIS SEEMS TO RUN FINE, BUT THEN CHECK FAILS
# .onLoad <- function(libname, pkgname)
# {
#   # Use good FFT library if available
#   if(is.installed("fftw"))
#   {
#     FFT <<- function(X) { fftw::FFT(as.numeric(X)) }
#     IFFT <<- function(X) { fftw::IFFT(as.numeric(X)) }
#   }
# }

zoom <- function(x,...) UseMethod("zoom") #S3 generic
#setGeneric("zoom",function(x,...) standardGeneric("zoom"),package="ctmm") #S4 generic

# forwarding function for list of a particular datatype
zoom.list <- function(x,...)
{
  CLASS <- class(x[[1]])
  utils::getS3method("zoom",CLASS)(x,...)
}
#methods::setMethod("zoom",signature(x="list"), function(x,...) zoom.list(x,...))


# forwarding function for list of a particular datatype
mean.list <- function(x,...)
{
  CLASS <- class(x[[1]])
  utils::getS3method("mean",CLASS)(x,...)
}
#methods::setMethod("mean",signature(x="list"), function(x,...) mean.list(x,...))


# forwarding function for list of a particular datatype
plot.list <- function(x,...)
{
  CLASS <- class(x[[1]])
  utils::getS3method("plot",CLASS)(x,...)
}
#methods::setMethod("plot",signature(x="list"), function(x,...) plot.list(x,...))


# forwarding function for list of a particular datatype
summary.list <- function(object,...)
{
  CLASS <- class(object[[1]])
  utils::getS3method("summary",CLASS)(object,...)
}


# parity tests
is.even <- Vectorize(function(x) {x %% 2 == 0})

is.odd <- Vectorize(function(x) {x %% 2 != 0})


# 2D rotation matrix
rotate <- function(theta)
{
  COS <- cos(theta)
  SIN <- sin(theta)
  R <- rbind( c(COS,-SIN), c(SIN,COS) )
  return(R)
}

# indices where a condition is met
where <- function(x)
{
  (1:length(x))[x]
}


# statistical mode
Mode <- function(x)
{
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
  mean(ux)
}

# adjoint of matrix
Adj <- function(M) { t(Conj(M)) }

# Hermitian part of matrix
He <- function(M) { (M + Adj(M))/2 }

# Positive definite part of matrix
PDpart <- function(M)
{ 
  # singular value decomposition
  M <- svd(M)
  M$d <- clamp(M$d,max=Inf) # toss out small negative values
  M$u <- (M$u + M$v)/2 # symmetrize
  
  M <- lapply(1:length(M$d),function(i){M$d[i]*(M$u[i,]%o%Conj(M$u[i,]))})
  M <- Reduce("+",M)
  
  return(M)
}

# Positive definite solver
PDsolve <- function(M)
{
  # symmetrize
  M <- He(M)
  
  # rescale
  W <- diag(M)
  W <- sqrt(W)
  W <- W %o% W
  
  # now a correlation matrix that is easy to invert
  M <- M/W
  M <- qr.solve(M)
  M <- M/W

  # symmetrize
  M <- He(M)

  return(M)
}

# COULD JUST MAKE A GENERAL PD-FUN OPERATION


# confidence interval functions
CI.upper <- Vectorize(function(k,Alpha){qchisq(Alpha/2,k,lower.tail=FALSE)/k})
CI.lower <- Vectorize(function(k,Alpha){qchisq(Alpha/2,k,lower.tail=TRUE)/k})


# calculate chi^2 confidence intervals from MLE and COV estimates
chisq.ci <- function(MLE,COV,alpha)
{
  DOF <- 2*MLE^2/COV
  CI <- MLE * c(CI.lower(DOF,alpha),1,CI.upper(DOF,alpha))
}


# last element of array
last <- function(vec) { vec[length(vec)] }


# CLAMP A NUMBER
clamp <- Vectorize(function(num,min=0,max=1) { if(num<min) {min} else if(num<max) {num} else {max} })


# PAD VECTOR
pad <- function(vec,size,padding=0,side="right")
{
  # this is now the pad length instead of total length
  size <- size - length(vec)
  
  if(side=="right"||side=="r")
  { return(c(vec,rep(padding,size))) }
  else if(side=="left"||side=="l")
  { return(c(rep(padding,size),vec)) }
}


# CHOOSE BEST UNITS FOR A LIST OF DATA
unit <- function(data,dimension,thresh=1)
{
  if(dimension=="length")
  {
    name.list <- c("meters","kilometers")
    scale.list <- c(1,1000)
  }
  else if(dimension=="area")
  {
    name.list <- c("square meters","hectares","square kilometers")
    scale.list <- c(1,100^2,1000^2) 
  }
  else if(dimension=="time")
  {
    name.list <- c("seconds","minutes","hours","days","years")
    scale.list <- c(1,60*c(1,60*c(1,24*c(1,365.24))))
  }
  else if(dimension=="speed")
  {
    name.list <- c("meters/day","kilometers/day")
    scale.list <- c(1,1000)/(60*60*24)
  }
  
  max.data <- max(abs(data))
  
  name <- name.list[1]
  scale <- scale.list[1]
  
  for(i in 2:length(name.list))
  {
    if(max.data > thresh*scale.list[i])
    {  
      name <- name.list[i]
      scale <- scale.list[i]
    }
  }
  
  return(list(scale=scale,name=name))
}


# convert units
setUnits <- function(arg1,arg2)
{
  return(arg1 %#% arg2) 
}


# convert units
`%#%` <- function(arg1,arg2)
{
  # convert to si units
  if(is.numeric(arg1))
  {
    num <- arg1
    name <- arg2
    pow <- +1
  }
  else # convert from si units
  {
    num <- arg2
    name <- arg1
    pow <- -1
  }
  
  alias <- list()
  scale <- c()
  
  add <- function(a,s)
  {
    n <- length(alias)
    alias[[n+1]] <<- a 
    scale[n+1] <<- s
  }
  
  # TIME
  add(c("s","s.","sec","sec.","second","seconds"),1)
  add(c("min","min.","minute","minutes"),60)
  add(c("h","h.","hr","hr.","hour","hours"),60^2)
  add(c("day","days"),24*60^2)
  add(c("week","weeks"),7*24*60^2)
  add(c("month","months"),365.24/12*7*24*60^2)
  add(c("yr","yr.","year","years"),365.24*7*24*60^2)
  
  # Distance conversions
  add(c("m","m.","meter","meters"),1)
  add(c("km","km.","kilometer","kilometers"),1000)
  
  # Area conversions
  add(c("m^2","m.^2","meter^2","meters^2","square meter","square meters","meter squared","meters squared"),1)
  add(c("ha","hectare","hectares"),100^2)
  add(c("km^2","km.^2","kilometer^2","kilometers^2","square kilometer","square kilometers","kilometer squared","kilometers squared"),1000^2)
  
  for(i in 1:length(alias))
  {
    if(name %in% alias[[i]]) { return(num*scale[i]^pow) }
  }
}