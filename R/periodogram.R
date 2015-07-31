# periodogram class
new.periodogram <- methods::setClass("periodogram",representation(info="list"),contains="data.frame")


# slow periodogram code
periodogram <- function(data,T=NULL,dt=NULL,resolution=1)
{
  t <- data$t
  x <- data$x
  y <- data$y
    
  # default sampling step
  if(is.null(dt)) { dt <- stats::median(diff(t)) }

  # Nyquist frequency
  F <- 1/(2*dt)
  
  # default sampling period
  if(is.null(T)) { T <- last(t)-t[1] }
  
  # frequency resolution
  df <- 1/(2*(T+dt)*resolution)
 
  # frequency grid
  f <- seq(df,F,df)
  
  # double angle matrix
  theta <- (4 * pi) * (f %o% t)

  # LSP lag shifts
  tau <- atan( rowSums(sin(theta)) / rowSums(cos(theta)) ) / (4*pi*f)
  
  # lagged angle matrix
  theta <- (2 * pi) * ((f %o% t) - (f * tau))

  # trigonometric matrices
  COS <- cos(theta)
  SIN <- sin(theta)
  
  LSP <- ((COS %*% x)^2 + (COS %*% y)^2)/rowSums(COS^2) + ((SIN %*% x)^2 + (SIN %*% y)^2)/rowSums(SIN^2)
  LSP <- LSP/4
  
  # sampling schedule periodogram
  SSP <- rowSums(COS)^2/rowSums(COS^2) + rowSums(SIN)^2/rowSums(SIN^2)
  SSP <- SSP/2
  
  # effective number of degrees of freedom for evenly sampled data
  # this is really just a placeholder for now
  DOF <- rep(1/resolution,length(f))
  
  result <- data.frame(LSP=LSP,SSP=SSP,DOF=DOF,f=f)
  result <- new.periodogram(result, info=attr(data,"info"))
  
  return(result)
}


# plot periodograms
plot.periodogram <- function(x,diagnostic=FALSE,col="black",transparency=0.25,...)
{
  # frequency in 1/days
  f <- x$f*24*60^2
  
  LSP <- x$LSP
  LSP <- log(LSP/max(LSP))
  
  SSP <- x$SSP
  SSP <- log(SSP/max(SSP))
  
  # human readable periods
  at <- c(365.24,365.24/4,29.53059,7,1)
  labels=c("year","season","month","week","day")
  for(d in 2:23)
  {
    at <- c(at,1/d)
    labels <- c(labels,paste("day/",d,sep=""))
  }
  at <- c(at,1/24)
  labels <- c(labels,"hour")
  
  plot(1/f,LSP,log="x",xaxt="n",xlab="Period",ylab="Log Spectral Density",col=translucent(col,alpha=((f[1]/f)^transparency)),...)
  if(diagnostic){ graphics::points(1/f,SSP,col=translucent("red",alpha=((f[1]/f)^transparency)),...) }
  graphics::axis(1,at=at,labels=labels)

}
#methods::setMethod("plot",signature(x="periodogram",y="missing"), function(x,y,...) plot.periodogram(x,...))
#methods::setMethod("plot",signature(x="periodogram"), function(x,...) plot.periodogram(x,...))
