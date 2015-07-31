\name{plot.telemetry}
\alias{plot.telemetry}
\alias{zoom.telemetry}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting methods for telemetry objects.
}
\description{
Produces simple plots of \code{telemetry} objects, possibly overlayed with a Gaussian \code{ctmm} movement model or an \code{akde} autocorrelated kernel density estimate.
}
\usage{
\S3method{plot}{telemetry}(x, CTMM=NULL, AKDE=NULL, alpha.HR=0.05, alpha=0.05, CI=TRUE, PDF=TRUE, col="red",
col.CI="black", col.PDF="blue", col.grid="grey", fraction=1, add=FALSE, xlim=NULL,
ylim=NULL, ...)

\S3method{zoom}{telemetry}(x,fraction=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ The 2D timeseries data represented as a \code{telemetry} object or list of such objects. All other arguments are optional.}
  \item{CTMM}{Optional Gaussian \code{ctmm} movement model from the output of \code{ctmm.fit} or list of such objects.}
  \item{AKDE}{Optional \code{akde} object from the output of \code{akde} or list of such objects.}
  \item{alpha.HR}{Significance level of \code{ctmm} model or \code{akde} density estimate contours to be displayed. I.e., \code{alpha.HR=0.50} yields the 50\% core home range within the rendered contours.}
  \item{alpha}{ Significant levels placed on the contour estimates themselves. I.e., the 50\% core home-range area can be estimated with 95\% confidence via \code{alpha=0.05}. }
  \item{CI}{Plot home-range contours.}
  \item{PDF}{Plot the maximum likelihood probability density function.}
  \item{col}{Color option for telemetry data. Can be an array.}
  \item{col.CI}{Color option for home-range contours. Can be an array.}
  \item{col.PDF}{Color option for the probability density function. Can be an array.}
  \item{col.grid}{Color option for the maximum likelihood \code{akde} bandwidth grid.}
  \item{fraction}{Fraction of the data, Gaussian \code{ctmm}, or \code{akde} range to plot, whichever is larger.}
  \item{add}{Setting to \code{TRUE} will disable the unit conversions and base layer plot, so that \code{plot.telemetry} can be overlayed atop other outputs more easily.}
  \item{xlim}{The \code{x} limits \code{c(x1, x2)} of the plot.}
  \item{ylim}{The \code{y} limits \code{c(y1, y2)} of the plot.}
  \item{\dots}{ Additional options passed to \code{plot}.}
}


\details{
Confidence intervals placed on the \code{ctmm} Gaussian home-range contour estimates only represent uncertainty in the area and not uncertainty in the mean location, eccentricity, or orientation angle.

For \code{akde} density estimates, the provided contours only represent the uncertainty in the optimal bandwidth area and their significance level is pre-calculated in the \code{akde} object itself, ignoring the \code{alpha} here.
With \code{akde} estimates, it is also important to note the scale of the bandwidth and, by default, grid cells are plotted with \code{akde} contours such that their length and width matches that of a bandwidth kernels' standard deviation in each direction. Therefore, this grid provides a visual approximation of the kernel-density estimate's ``resolution''.}

\value{Returns a plot of \eqn{x} vs. \eqn{y}, and, if specified, Gaussian \code{ctmm} or \code{akde} home range. \code{akde} plots also come with a standard resolution grid. \code{zoom} includes a zoom slider to manipulate \code{fraction}.}

\references{
C. H. Fleming, J. M. Calabrese, T. Mueller, K.A. Olson, P. Leimgruber, and W. F. Fagan. (2014). From fine-scale foraging to home ranges: A semi-variance approach to identifying movement modes across spatiotemporal scales. \href{http://www.jstor.org/discover/10.1086/675504}{The American Naturalist, 183(5), E154-E167.}

C. H. Fleming and J. M. Calabrese and T. Mueller and K. A. Olson and P. Leimgruber and W. F. Fagan (2014).
Non-Markovian maximum likelihood estimation of autocorrelated movement processes
\href{http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12176/abstract}{Methods in Ecology and Evolution, 5(5) 462-472}.

C. H. Fleming and W. F. Fagan and T. Mueller and K. A. Olson and P. Leimgruber and J. M. Calabrese (2015).
Rigorous home-range estimation with movement data: A new autocorrelated kernel-density estimator.
\href{http://www.esajournals.org/doi/abs/10.1890/14-2010.1}{Ecology, doi:10.1890/14-2010.1}.
}

\author{C. H. Fleming.}

\note{If \code{xlim} or \code{ylim} are provided, then the smaller or absent range will be expanded to ensure \code{asp=1}.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{ctmm.fit}}, \code{\link{plot}}, \code{\link{SpatialPoints.telemetry}}.  }

\examples{
# Load package and data
library(ctmm)
data(buffalo)

# Plot the data
plot(buffalo)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
