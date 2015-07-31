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
akde(data,CTMM,alpha=0.05,fast=NULL,dt=NULL,error=0.005)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ 2D timeseries telemetry data represented as a \code{telemetry} object or data.frame or list with components \code{t}, \code{x}, \code{y} all numeric values. }
  \item{CTMM}{A \code{ctmm} movement model from the output of \code{ctmm.fit}.}
  \item{alpha}{Significance level for when propagating uncertainties in the \code{ctmm} parameter estimates.}
  \item{fast, dt}{Arguments passed to all instances of \code{\link{akde.bandwidth}}.}
  \item{error}{Target error for discretization.}
}

\details{
The grid size and resolution is chosen to prevent probability errors larger than \code{error} on a worst case basis.

The grid size is chosen such that locations nearest the boundary will have at least \code{error} portion of their kernel's probability within an ellipse contained by the grid. For this the grid boundary is chosen to be at least a number of bandwidth standard deviations, \eqn{z}, from all sampled locations according to the 2D Gaussian quantile relation \eqn{z^2 =  -2 \log} \code{error}.

The choice of grid cell size is such that no single pixel contains a significant amount of probability by constraining the discretized probability density according to \eqn{ 2 \sum_{i=1}^\infty \frac{\exp{-\frac{dz^2}{2}}}{\sqrt{2\pi}} = 1 - } \code{error}. Limiting behavior of the Jacobi theta function then gives a suggested grid resolution of \eqn{dz = } \code{error} \eqn{/ 0.398942\cdots} bandwidth standard deviations.
}

\value{Returns an \code{akde} object: a list with component kernel-density objects \code{low}, \code{ML}, \code{high}, and the corresponding significance level \code{alpha} for the optimal bandwidth area. Each kernel-density object contains the sampled grid line locations \code{x} and \code{y}, the probability density function evaluated on the sampled grid locations \code{pdf}, the optimal bandwidth matrix \code{H}, and the area of each grid cell \code{dA}. }

\references{
C. H. Fleming and W. F. Fagan and T. Mueller and K. A. Olson and P. Leimgruber and J. M. Calabrese (2015).
Rigorous home-range estimation with movement data: A new autocorrelated kernel-density estimator.
\href{http://www.esajournals.org/doi/abs/10.1890/14-2010.1}{Ecology, 96(5), 1182-1188}.

%D. Marcotte. (1996). Fast variogram computation with FFT. \emph{Computers and Geosciences} 22(10), 1175-1186.
}

\author{ Chris H. Fleming. }
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde.bandwidth}}, \code{\link{raster.akde}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)
cilla <- buffalo[[1]]

# Fit a continuous-velocity model with tau ~ c(10 days, 1 hour)
# see help(variogram.fit)
m2 <- ctmm(tau=c(10*24*60^2,60^2))
M2 <- ctmm.fit(cilla,m2)

# Compute akde object
KD2 <- akde(cilla,M2)

# Plot data with AKDE contours
plot(cilla,AKDE=KD2)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
