\name{as.telemetry}
\alias{as.telemetry}
\alias{summary.telemetry}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Import, coerce, and summarize MoveBank data
}
\description{Functions to import MoveBank csv files or \code{data.frame} objects, coerce them into \code{telemetry} objects, and summarize them.}

\usage{
as.telemetry(CSV,timezone="GMT",projection=NULL,...)

\S3method{summary}{telemetry}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CSV}{The file name of the MoveBank csv file or a raw \code{data.frame} object from previously read MoveBank data.}
  \item{timezone}{Timezone argument for \code{\link{as.POSIXct}}.}
  \item{projection}{Optional PROJ.4 projection argument to be fed to \code{\link[rgdal]{project}}.}
  \item{...}{Optional arguments to be fed to \code{\link{read.csv}}.}
  \item{object}{\code{telemetry} object to summarize.}
}

\details{
If no projection argument is specified, a two-point equidistant projection is calculated that should be good for most range resident and migratory species.
Global migrations that are not relatively geodesic (locally straight) will probably suffer distortion.

The current MoveBank columns used are \code{individual.local.identifier}, \code{timestamp}, \code{location.long}, \code{location.lat}.
In the future, we will be supporting more GPS and ARGOS columns, including speed, heading, error classification, etc.
}

\value{\code{as.telemetry} returns a single \code{telemetry} object or list of \code{telemetry} objects if multiple animals are identified.}

%\references{}

\author{C. H. Fleming and K. Safi.}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{plot.telemetry}}, \code{\link{SpatialPoints.telemetry}}.  }

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
