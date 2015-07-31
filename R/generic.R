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


# parity tests
is.even <- Vectorize(function(x) {x %% 2 == 0})

is.odd <- Vectorize(function(x) {x %% 2 != 0})


# statistical mode
Mode <- function(x)
{
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
  mean(ux)
}


# confidence interval functions
CI.upper <- Vectorize(function(k,Alpha){qchisq(Alpha/2,k,lower.tail=FALSE)/k})
CI.lower <- Vectorize(function(k,Alpha){qchisq(Alpha/2,k,lower.tail=TRUE)/k})


# calculate chi^2 confidence intervals from MLE and COV estimates
chisq.ci <- function(MLE,COV,alpha)
{
  DOF <- 2*MLE^2/COV
  CI <- MLE * c(CI.lower(DOF,alpha),1,CI.upper(DOF,alpha))
}


# add opacity to color
translucent <- Vectorize( function(col,alpha=0.5)
{
  col <- col2rgb(col)/255
  col <- grDevices::rgb(col[1],col[2],col[3],alpha)
  return(col)
} )


# last element of array
last <- function(vec) { vec[length(vec)] }


# CLAMP A NUMBER
clamp <- Vectorize(function(num,min.num=0,max.num=1) { if(num<min.num) {min.num} else if(num>max.num) {max.num} else {num} })


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