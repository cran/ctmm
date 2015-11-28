\name{akde}
\alias{akde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate an autocorrelated kernel density estimate
}
\description{
This function calculates autocorrelated kernel density estimates of different confidence levels from \code{telemetry} data and a continuous-time movement model.
}
\usage{
akde(data,CTMM,level=0.95,error=0.001,res=200,grid=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ 2D timeseries telemetry data represented as a \code{telemetry} object. }
  \item{CTMM}{A \code{ctmm} movement model from the output of \code{ctmm.fit}.}
  \item{level}{Confidence level for when propagating uncertainties in the \code{ctmm} parameter estimates.}
  \item{error}{Target error for kernel truncation.}
  \item{res}{Number of grid points along each axis.}
  \item{grid}{Optional grid specification with columns labeled \code{x} and \code{y}. Not yet supported.}
  \item{...}{Arguments passed to all instances of \code{\link{akde.bandwidth}}.}
}

%\details{}

\value{Returns a \code{UD} object: a list with component objects \code{low}, \code{ML}, \code{high}, and the corresponding confidence \code{level} for the optimal bandwidth area. Each component object contains the sampled grid line locations \code{x} and \code{y}, the probability density and cumulative distribution functions evaluated on the sampled grid locations \code{PDF} & \code{CDF}, the optimal bandwidth matrix \code{H}, and the area of each grid cell \code{dA}.
}

\references{
C. H. Fleming and W. F. Fagan and T. Mueller and K. A. Olson and P. Leimgruber and J. M. Calabrese (2015).
Rigorous home-range estimation with movement data: A new autocorrelated kernel-density estimator.
\href{http://www.esajournals.org/doi/abs/10.1890/14-2010.1}{Ecology, 96(5), 1182-1188}.

%D. Marcotte. (1996). Fast variogram computation with FFT. \emph{Computers and Geosciences} 22(10), 1175-1186.
}

\author{ C. H. Fleming. }

\note{ In the case of coarse grids, the value of \code{PDF} in a grid cell corresponds to the average probability density over the entire rectangular cell. }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde.bandwidth}}, \code{\link{raster.UD}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)
cilla <- buffalo[[1]]

# Fit a continuous-velocity model with tau ~ c(10 days, 1 hour)
# see help(variogram.fit)
GUESS <- ctmm(tau=c(10*24*60^2,60^2))
FIT <- ctmm.fit(cilla,GUESS)

# Compute akde object
UD <- akde(cilla,FIT)

# Plot data with AKDE contours
plot(cilla,UD=UD)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
