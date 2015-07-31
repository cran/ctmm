\name{buffalo}
\alias{buffalo}
\docType{data}

\title{African buffalo GPS dataset from Kruger National Park, South Africa.}

\description{GPS data on six African buffalo, provided for demonstration purposes only! Contact data owner Paul Cross if you wish to do more with the data.}

\usage{data("buffalo")}

\format{ A list of 6 \code{telemetry} objects.}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\seealso{ \code{\link{as.telemetry}}, \code{\link{plot.telemetry}}. }

\examples{
# Load package and data
library(ctmm)
data(buffalo)

# Extract movement data for a single animal
cilla <- buffalo[[1]]

# Plot all sampled locations
plot(cilla)
}

\keyword{datasets}
