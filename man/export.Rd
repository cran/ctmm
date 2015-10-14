\name{ctmm::export}
\alias{raster.akde}
\alias{SpatialPoints.telemetry}
\alias{SpatialPolygonsDataFrame.akde}
\alias{writeShapefile.akde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Export ctmm data formats
}
\description{Functions to export \code{ctmm} data formats into common \code{sp}, \code{raster}, and ESRI formats.}

\usage{
raster.akde(AKDE,CI="ML")

SpatialPoints.telemetry(data)

SpatialPolygonsDataFrame.akde(AKDE,alpha.HR=0.05)

writeShapefile.akde(AKDE,folder,file=AKDE@info$identity,alpha.HR=0.05,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{AKDE}{ \code{akde} kernel density object. }
  \item{CI}{Which confidence interval of the bandwidth uncertainty to use: the "ML" point estimate, "low" or "high".}
  \item{data}{ \code{telemetry} data object. }
  \item{alpha.HR}{ Significance level of the home-range area. I.e., the 50\% core home range is given by \code{alpha.HR=0.50}. }
  \item{folder}{Character name of folder for shapefile.}
  \item{file}{Character name of files for shapefile.}
  \item{...}{ Optional arguments passed to \code{\link[rgdal]{writeOGR}}.}
}

\details{\code{writeShapefile.akde} writes a shapefile to disk, with polygons corresponding to the low, ML, and high optimal-bandwidth AKDE home-range areas.}

\value{\code{raster.akde} returns a \code{raster} of the maximum likelihood (ML) probability density function.

\code{SpatialPoints.telemetry} returns a \code{spatialPoints} object for the \code{x}-\code{y} locations.

\code{SpatialPolygonsDataFrame.akde} returns a \code{SpatialPolygonsDataFrame} of the low, ML, and high optimal-bandwidth AKDE home-range areas, in the appropriate order for plotting.
}

%\references{}

\author{C. H. Fleming and K. Safi.}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{as.telemetry}}.  }

%\examples{
%# Load package and data
%library(ctmm)
%data(buffalo)

%# Extract movement data for a single animal
%billy <- buffalo[[1]]

%# Plot the data
%plot.telemetry(billy)
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
