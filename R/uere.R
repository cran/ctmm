uere <- function(data,diagnostic=FALSE)
{
  # promote to list of 
  if(class(data)=="telemetry" || class(data)=="data.frame") { data <- list(data)  }
  
  n <- 0
  mu <- list()
  W <- 0
  r <- 0
  # calculate mean and residuals
  for(i in 1:length(data))
  {
    if(is.null(data[[i]]$HDOP)) { stop("Failed to import GPS.HDOP column from Movebank file || missing HDOP column in telemetry object.") }
    
    w <- 1/data[[i]]$HDOP^2
    W[i] <- sum(w)
    n <- n + 2*length(w)
    mu[[i]] <- c( sum(w*data[[i]]$x) , sum(w*data[[i]]$x) )/W[i]
    r <- r + sum(w*(data[[i]]$x - mu[[i]][1])^2) + sum(w*(data[[i]]$y - mu[[i]][2])^2)
  }
  
  # total degrees of freedom
  n <- n - 2*length(data)
  
  # residuals and UERE
  UERE <- sqrt(r/n)
  
  if(diagnostic)
  {
    CTMM <- list()
    for(i in 1:length(data))
    {
      data[[i]]$error <- (UERE*data[[i]]$HDOP)^2/2
      sigma <- UERE^2/W[i]/2
      COV <- array(sigma^2 * 2/n,c(1,1))
      dimnames(COV) <- list("area","area")
      CTMM[[i]] <- ctmm(mu=mu[[i]],sigma=sigma,COV=COV,isotropic=TRUE)
    }
    
    plot(data,CTMM=CTMM,col.DF="black",col=grDevices::rainbow(length(data)))
    
    RETURN <- list()
    RETURN$UERE <- UERE
    RETURN$UERE95CI <- chisq.ci(UERE,DOF=n)
    
    return(RETURN)
  }
  else
  { return( UERE ) }
}