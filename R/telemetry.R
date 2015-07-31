#########################################
# telemetry class definition
#########################################
new.telemetry <- methods::setClass("telemetry", representation(info="list"), contains="data.frame")


#######################
# Import movebank into telemetry class
as.telemetry <- function(CSV,timezone="GMT",projection=NULL,...)
{
  type <- class(CSV)
  if (type=="data.frame") { CSV <- telemetry.data.frame(CSV,timezone=timezone,projection=projection) }
  else if (type=="character") { CSV <- telemetry.csv(CSV,timezone=timezone,projection=projection,...) }
  
  return(CSV)
}


# this assumes a MoveBank data.frame
telemetry.data.frame <- function(csv,timezone="GMT",projection=NULL)
{
  csv <- data.frame(id=csv$individual.local.identifier,timestamp=as.POSIXct(csv$timestamp,tz=timezone),longitude=csv$location.long,latitude=csv$location.lat)
  csv <- stats::na.omit(csv)
  csv$t <- as.numeric(csv$timestamp)
  
  xy <- cbind(csv$lon,csv$lat)
  colnames(xy) <- c("x","y")
  
  if(is.null(projection)) { projection <- suggest.projection(csv) }
  xy <- rgdal::project(xy,projection)
  
  csv$x <- xy[,1]
  csv$y <- xy[,2]
  
  id <- levels(csv$id)
  n <- length(id)
  
  telist <- list()  
  for(i in 1:n)
  {
    telist[[i]] <- csv[csv$id==id[i],]
    telist[[i]]$id <- NULL
    
    # clean through duplicates, etc..
    telist[[i]] <- telemetry.clean(telist[[i]])
        
    # combine data.frame with ancillary info
    info <- list(identity=id[i] ,timezone=timezone , projection=projection)
    telist[[i]] <- new.telemetry( telist[[i]] , info=info )
  }
  names(telist) <- id
  
  if (n>1) { return(telist) }
  else { return(telist[[1]]) }
}


#################
# clean up data
telemetry.clean <- function(data)
{  
  # sort in time
  data <- data[sort.list(data$t,na.last=NA,method="quick"),]
  
  # remove duplicate observations
  data <- unique(data)
  
  # exit with warning on duplicate times
  
}


# read in a MoveBank csv file
telemetry.csv <- function(file,timezone="GMT",projection=NULL,...)
{
  data <- utils::read.csv(file,...)
  data <- telemetry.data.frame(data,timezone=timezone,projection=projection)
  return(data)
}


########################################
# Suggest a good projection
########################################
suggest.projection <- function(data,datum="WGS84")
{
  # assume Movebank data.frame
  lon <- data$longitude
  lat <- data$latitude
  
  # as a first approximation use one-point equidistant at mean geolocation
  lon_0 <- mean(lon)
  lat_0 <- mean(lat)
  proj <- paste("+proj=aeqd +lon_0=",lon_0," +lat_0=",lat_0," +datum=",datum,sep="")
  xy <- rgdal::project(cbind(lon,lat),proj)
  
  # calculate and detrend mean
  mu <- c(mean(xy[,1]),mean(xy[,2]))
  xy <- xy - mu
  colnames(xy) <- c("x","y")
  
  # cross correlation
  cov <- mean(xy[,1]*xy[,2])
  # covariance matrix
  cov <- rbind( c( mean(xy[,1]^2) , cov ) , c( cov , mean(xy[,2]^2) ) )
  
  # figure out long axis (always first dim)
  R <- eigen(cov)$vectors
  # rotate data to long axis
  xy <- xy %*% R
  
  # bi-modal split of data
  xy1 <- xy[xy[,1]<0,]
  xy2 <- xy[xy[,1]>0,]
  
  # bi-modal modes
  mu1 <- c(mean(xy1[,1]),mean(xy1[,2]))
  mu2 <- c(mean(xy2[,1]),mean(xy2[,2]))
  
  # reverse rotation
  R <- solve(R)
  mu1 <- mu1 %*% R
  mu2 <- mu2 %*% R
  
  # re-trend mean
  mu1 <- mu1 + mu
  mu2 <- mu2 + mu
  
  # get long lat
  mu1 <- rgdal::project(mu1,proj,inv=TRUE)[1,]
  mu2 <- rgdal::project(mu2,proj,inv=TRUE)[1,]
  
  if(mu1[1] > mu2[1])
  {
    mu <- mu1
    mu1 <- mu2
    mu2 <- mu
  }
  
  proj <- paste("+proj=tpeqd +lon_1=",mu1[1]," +lat_1=",mu1[2]," +lon_2=",mu2[1]," +lat_2=",mu2[2]," +datum=",datum,sep="")
  
  return(proj)
}

################################
# ZOOM INTO TELEMETRY DATA
zoom.telemetry <- function(x,fraction=1,...)
{
  manipulate::manipulate(
  { plot.telemetry(x,fraction=fraction,...) },
  fraction = manipulate::slider(0, 1.0,initial=fraction)
  ) 
}
#methods::setMethod("zoom",signature(x="telemetry",y="missing"), function(x,y,...) zoom.telemetry(x,...))
#methods::setMethod("zoom",signature(x="telemetry",y="telemetry"), function(x,y,...) zoom.telemetry(list(x,y),...))
#methods::setMethod("zoom",signature(x="telemetry",y="ctmm"), function(x,y,...) zoom.telemetry(x,model=y,...))
#methods::setMethod("zoom",signature(x="telemetry",y="akde"), function(x,y,...) zoom.telemetry(x,akde=y,...))
#methods::setMethod("zoom",signature(x="telemetry"), function(x,...) zoom.telemetry(x,...))


#######################################
# PLOT TELEMETRY DATA
#######################################
plot.telemetry <- function(x,CTMM=NULL,AKDE=NULL,alpha.HR=0.05,alpha=0.05,CI=TRUE,PDF=TRUE,col="red",col.CI="black",col.PDF="blue",col.grid="grey",fraction=1,add=FALSE,xlim=NULL,ylim=NULL,...)
{
  # adjustments to CI-CIs to make them a bit diminished relatively
  trans <- c(0.5,1,0.5)
  lwd <- c(1,2,1)
  
  # listify everything for generality
  if(class(x)=="telemetry" || class(x)=="data.frame") { x <- list(x)  }
  if(!is.null(CTMM)) { if(class(CTMM)=="ctmm") { CTMM <- list(CTMM) } }
  if(!is.null(AKDE)) { if(class(AKDE)=="akde") { AKDE <- list(AKDE) } }
    
  dist.name <- "meters"
  dist.scale <- 1
  if(!add)
  {
    # bounding locations from data
    ext.x <- min(sapply(x, function(d){ min(d$x) } ))
    ext.x[2] <- max(sapply(x, function(d){ max(d$x) } ))
    
    ext.y <- min(sapply(x, function(d){ min(d$y) } ))
    ext.y[2] <- max(sapply(x, function(d){ max(d$y) } ))
    
    # bounding locations from AKDEs
    if(!is.null(AKDE))
    {
      ext.x[1] <- min(c(ext.x[1], sapply(AKDE, function(a) { sapply(a[1:3], function(k) { k$x[1] }) })))
      ext.x[2] <- max(c(ext.x[2], sapply(AKDE, function(a) { sapply(a[1:3], function(k) { last(k$x) }) })))
      
      ext.y[1] <- min(c(ext.y[1], sapply(AKDE, function(a) { sapply(a[1:3], function(k) { k$y[1] }) })))
      ext.y[2] <- max(c(ext.y[2], sapply(AKDE, function(a) { sapply(a[1:3], function(k) { last(k$y) }) })))
    }
    
    # bounding locations from Gaussian CTMM
    if(!is.null(CTMM))
    {
      z <- sqrt(-2*log(alpha.HR))
      
      for(i in 1:length(CTMM))
      {
        # proportionality constants for outer CIs
        sigma <- CTMM[[i]]$sigma

        # capture outer contour if present
        const <- 1
        if(!is.null(CTMM[[i]]$COV.tau))
        {
          K <- length(CTMM[[i]]$tau)
          const <- confint.ctmm(CTMM[[i]],alpha)[1,3]/sqrt(det(sigma))
        }

        buff <- z*sqrt(const*diag(sigma))
        
        ext.x[1] <- min(ext.x[1], CTMM[[i]]$mu[1] - buff[1])
        ext.x[2] <- max(ext.x[2], CTMM[[i]]$mu[1] + buff[1])
        
        ext.y[1] <- min(ext.y[1], CTMM[[i]]$mu[2] - buff[2])
        ext.y[2] <- max(ext.y[2], CTMM[[i]]$mu[2] + buff[2])
      }
    }
        
    # bounding box
    mu <- c(mean(ext.x),mean(ext.y))
    buff <- c(diff(ext.x),diff(ext.y))/2
    
    # now zoom in/out to some fraction of the grid
    buff <- fraction*buff
    
    ext.x <- mu[1] + buff[1]*c(-1,1)
    ext.y <- mu[2] + buff[2]*c(-1,1)
    
    # try to obey xlim/ylim if provided
    if(!is.null(xlim) || !is.null(ylim))
    {
      max.diff <- max(diff(xlim),diff(ylim))*c(-1,1)/2
      
      if(is.null(ylim))
      { ylim <- mu[2] + max.diff }
      else if(is.null(xlim))
      { xlim <- mu[1] + max.diff }

      ext.x <- xlim
      ext.y <- ylim
    }
    
    # Get best unit scale
    dist <- unit(abs(c(ext.x,ext.y)),"length")
    dist.name <- dist$name
    dist.scale <- dist$scale
    
    xlab <- paste("x ", "(", dist.name, ")", sep="")
    ylab <- paste("y ", "(", dist.name, ")", sep="")
    
    mu <- mu/dist.scale
    
    ext.x <- ext.x/dist.scale
    ext.y <- ext.y/dist.scale
    
    # empty base layer plot
    plot(ext.x,ext.y, xlab=xlab, ylab=ylab, col=grDevices::rgb(1,1,1,0), asp=1, ...)
  }
  
  #########################
  # PLOT GAUSSIAN CONTOURS AND DENSITY
  if(!is.null(CTMM))
  {
    # number of CTMM objects
    
    # contours colour
    col.CI <- array(col.CI,length(CTMM))
    col.PDF <- array(col.PDF,length(CTMM))
    
    for(i in 1:length(CTMM))
    {
      tau <- CTMM[[i]]$tau
      K <- length(tau)
      
      # scale coordinates
      CTMM[[i]]$mu <- CTMM[[i]]$mu/dist.scale
      CTMM[[i]]$sigma <- CTMM[[i]]$sigma/dist.scale^2
      
      # plot ML density function
      if(PDF)
      {
        # look how lazy I am
        PDF <- kde(list(x=CTMM[[i]]$mu[1],y=CTMM[[i]]$mu[2]),H=CTMM[[i]]$sigma)
        plot.pdf(PDF,col=col.PDF[[i]],...)
      }
      
      # plot CIs
      if(CI)
      {
        # plot ML estimate, regular style
        plot.ctmm(CTMM[[i]],alpha.HR,col=col.CI[[i]],lwd=lwd[2],...)
        
        # plot CIs dashed if present
        if(!is.null(CTMM[[i]]$COV.tau))
        {
          CTMM[[i]]$COV.tau <- CTMM[[i]]$COV.tau/dist.scale^4 # don't care about tau, just sigma
          
          # proportionality constants for outer CIs
          const <- confint.ctmm(CTMM[[i]],alpha)[1,c(1,3)]/sqrt(det(CTMM[[i]]$sigma))
          sigma <- CTMM[[i]]$sigma
          
          for(j in 1:2)
          {
            CTMM[[i]]$sigma <- const[j]*sigma
            plot.ctmm(CTMM[[i]],alpha.HR,col=translucent(col.CI[[i]],trans[1]),...)
          }
        }
      }
    
    }
        
  }

  ##########################################
  # PLOT akde CONTOURS... AND DENSITY
  if(!is.null(AKDE))
  {
    # number of akde objects
    
    # contours colour
    col.CI <- array(col.CI,length(AKDE))
    col.PDF <- array(col.PDF,length(AKDE))
    
    # UNIT CONVERSIONS
    for(i in 1:length(AKDE))
    {
      for(j in 1:3)
      {
        # unit conversion
        AKDE[[i]][[j]]$x <- AKDE[[i]][[j]]$x / dist.scale
        AKDE[[i]][[j]]$y <- AKDE[[i]][[j]]$y / dist.scale
        AKDE[[i]][[j]]$pdf <- AKDE[[i]][[j]]$pdf * dist.scale^2
        AKDE[[i]][[j]]$dA <- AKDE[[i]][[j]]$dA / dist.scale^2
        AKDE[[i]][[j]]$H <- AKDE[[i]][[j]]$H / dist.scale^2
      }
      
      # ML DENSITY PLOT
      if(PDF) { plot.pdf(AKDE[[i]]$ML,col=col.PDF[[i]],...) }
    }
    
    # CONTOURS
    if(CI)
    {
      for(i in 1:length(AKDE))
      {
        for(j in 1:3)
        {
          plot.kde(AKDE[[i]][[j]],alpha=alpha.HR,col=translucent(col.CI[[i]],trans[j]),lwd=lwd[j],...)
        }
      }
    }
          
    # RESOLUTION GRID
    if(!add)
    {
      dx <- sqrt(max(sapply( AKDE , function(a) { a$ML$H[1,1] } )))
      dy <- sqrt(max(sapply( AKDE , function(a) { a$ML$H[2,2] } )))
      
      # grid points per half
      gp <- ceiling(buff/c(dx,dy))
      col.grid <- translucent(col.grid,0.5)
      graphics::abline(v=(mu[1]+i*dx*(-gp[1]:gp[1])), col=col.grid)
      graphics::abline(h=(mu[2]+i*dy*(-gp[2]:gp[2])), col=col.grid)
    }
  }
  
  #########################
  # PLOT TELEMETRY DATA
  
  # color array for plots
  col <- array(col,length(x))
  
  # automagic the plot point size
  p <- sum(sapply(x, function(d) { length(d$t) } ))
  cex <- 1
  if(p>1000) { cex <- 1000/p }
  
  for(i in 1:length(x))
  {    
    graphics::points(x[[i]][,c("x","y")]/dist.scale, cex=cex, col=col[[i]],...)
  }
  
}
# SET METHODS FOR PLOT.TELEMETRY
#methods::setMethod("plot",signature(x="telemetry",y="missing"), function(x,y,...) plot.telemetry(x,...))
#methods::setMethod("plot",signature(x="telemetry",y="telemetry"), function(x,y,...) plot.telemetry(list(x,y),...))
#methods::setMethod("plot",signature(x="telemetry",y="ctmm"), function(x,y,...) plot.telemetry(x,model=y,...))
#methods::setMethod("plot",signature(x="telemetry",y="akde"), function(x,y,...) plot.telemetry(x,akde=y,...))
#methods::setMethod("plot",signature(x="telemetry"), function(x,...) plot.telemetry(x,...))


##################################
# plot PDF stored as KDE object
plot.pdf <- function(kde,col="blue",...)
{
  col <- translucent(col,alpha=(0:255)/255)
  zlim <- c(0,max(kde$pdf))
  graphics::image(kde$x,kde$y,kde$pdf,useRaster=TRUE,zlim=zlim,col=col,add=TRUE,...)
}


#############################
# Plot a KDE object's contours
plot.kde <- function(kde,alpha=0.05,col="black",...)
{
  # probability density of contour...
  p <- qkde(kde,alpha=alpha)

  graphics::contour(x=kde$x,y=kde$y,z=kde$pdf,levels=p,labels=round((1-alpha)*100),labelcex=1,col=col,add=TRUE,...)
}


##############################
# Plot Gaussian ctmm contours
plot.ctmm <- function(model,alpha=0.05,col="blue",...)
{
  mu <- model$mu
  sigma <- model$sigma
  
  Eigen <- eigen(sigma)
  std <- sqrt(Eigen$values)
  vec <- Eigen$vectors
  
  z <- sqrt(-2*log(alpha))
  
  num <- 100
  theta <- 2*pi*(0:num)/(num+1)
  Sin <- sin(theta)
  Cos <- cos(theta)
  
  x <- mu[1] + z*(Cos*std[1]*vec[1,1] + Sin*std[2]*vec[1,2])
  y <- mu[2] + z*(Cos*std[1]*vec[2,1] + Sin*std[2]*vec[2,2])
  
  graphics::xspline(x, y=y, shape=-1, open=FALSE, border=col, ...)
}


########################
# summarize telemetry data
summary.telemetry <- function(object,...)
{
  result <- attr(object,"info")
  
  dt <- stats::median(diff(object$t))
  units <- unit(dt,"time",thresh=1)
  result <- c(result,dt/units$scale)
  names(result)[length(result)] <- paste("sampling interval (",units$name,")",sep="")
  
  T <- last(object$t)-object$t[1]
  units <- unit(T,"time",thresh=1)
  result <- c(result,T/units$scale)
  names(result)[length(result)] <- paste("sampling period (",units$name,")",sep="")
  
  lon <- c(min(object$longitude),max(object$longitude))
  result <- c(result,list(lon=lon))
  names(result)[length(result)] <- paste("longitude range")
  
  lat <- c(min(object$latitude),max(object$latitude))
  result <- c(result,list(lat=lat))
  names(result)[length(result)] <- paste("latitude range")
  
  return(result)
}
#methods::setMethod("summary",signature(object="telemetry"), function(object,...) summary.telemetry(object,...))


#########################
# convert to spatialpoints object
SpatialPoints.telemetry <- function(data)
{
  if(class(data)=="telemetry" || class(data)=="data.frame")
  {
    return( sp::SpatialPoints( data[c("x","y")], proj4string=sp::CRS(attr(data,"info")$projection) ) )
  }
  else if(class(data)=="list")
  {
    SPL <- lapply( data, function(d) { sp::SpatialPoints( d[c("x","y")], proj4string=sp::CRS(attr(d,"info")$projection) ) } )
    return(SPL)
  }
}
#methods::setMethod("SpatialPoints",signature(coords="telemetry"), function(coords,...) SpatialPoints.telemetry(coords,...))


######################
# clean up of buffalo data
# this point is way off
#buffalo[[6]] <- ctmm:::new.telemetry(buffalo[[6]][-5720,],info=attr(buffalo[[6]],"info"))
# this time is duplicated and much less likely than the first
#buffalo[[5]] <- ctmm:::new.telemetry(buffalo[[5]][-869,],info=attr(buffalo[[5]],"info"))