\name{ctmm.loglike}
\alias{ctmm.loglike}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the log-likelihood of a continuous-time movement model
}
\description{
This function evaluates an \eqn{O(n)} Kalman-filter implementation of the exact non-Markovian likelihood function of continuous-time movement models described in Fleming et al (2014) and Fleming et al (2015), given 2D animal tracking data.
}
\usage{
ctmm.loglike(data,CTMM,verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ \code{telemetry} data object for the 2D timeseries. }
  \item{CTMM}{ \code{ctmm} model object with the necessary component \code{tau} and optional component \code{isotropic}. }
  \item{verbose}{Return a \code{ctmm} model object, including the log-likelihood, instead of only the log-likelihood.}
}
\details{
See \code{\link{ctmm.fit}} for a complete description of the parameter inputs.
}

\value{
By default, only the log-likelihood value is returned, given the specified model parameters and data.
If \code{verbose=TRUE}, then a \code{ctmm} movement-model object list is returned with the following components:
\describe{
\item{\code{loglike}}{The log-likelihood.}
\item{\code{isotropic}}{A Boolean denoting whether or not the animal's covariance is circular or elliptical.}
\item{\code{tau}}{The argument \code{CTMM$tau}.}
\item{\code{sigma}}{The argument \code{CTMM$sigma} if specified, or the analytically profiled variance/covariance if not specified.}
\item{\code{mu}}{The argument \code{CTMM$mu} if specified, or the analytically profiled stationary mean vector if not.}
\item{\code{COV.mu}}{The covariance matrix of the vector \code{mu}, assuming that the point estimates \code{tau} and \code{sigma} are good.}
\item{\code{DOF.mu}}{The effective number of degrees of freedom in the estimate of \code{mu}, assuming that the point estimate \code{tau} is good.
This can be much smaller than \code{length(data$t)} if the data are autocorrelated.}
}}

\references{
C. H. Fleming, J. M. Calabrese, T. Mueller, K.A. Olson, P. Leimgruber, and W. F. Fagan. (2014). From fine-scale foraging to home ranges: A semi-variance approach to identifying movement modes across spatiotemporal scales. \href{http://www.jstor.org/discover/10.1086/675504}{The American Naturalist, 183(5), E154-E167.}

C. H. Fleming and Y. Subasi and J. M. Calabrese. (2015). A maximum-entropy description of animal movement. \href{http://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.032107}{Physical Review E, 91, 032107.}
}

\author{ C. H. Fleming. }
% \note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{ctmm.fit}} }

\examples{
#Load package and data
library(ctmm)
data(buffalo)

#Extract movement data for a single animal
cilla <- buffalo[[1]]

#Calculate likelihood of ~10 day range crossing time with discontinuous velocity
ctmm.loglike(cilla,ctmm(tau=10*24*60*60),verbose=TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
