\name{ctmm.fit}
\alias{ctmm.fit}
\alias{ctmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a continuous-time movement model
}
\description{
This function is a wrapper around \code{opim} and \code{ctmm.loglike} to maximize the likelihood function of continuous-time movement models described in Fleming et al (2014) and Fleming et al (2015), given 2D animal tracking data.
}
\usage{
ctmm.fit(data,CTMM=NULL,...)

ctmm(tau=NULL,isotropic=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ The 2D timeseries data represented as a \code{telemetry} object. }
  \item{CTMM}{ A \code{ctmm} movement-model object containing the initial guesses for \code{tau} and optional component \code{isotropic}.}
  \item{...}{Arguments passed to \code{\link{optim}}.}
  \item{tau}{Array of autocorrelation timescales explained below.}
  \item{isotropic}{A Boolean denoting whether or not the animal's covariance is circular or elliptical.}
}

\details{The initial \code{ctmm} parameter guess can be the output of \code{variogram.fit} or the function \code{ctmm(...)} with the argument \code{tau} explained below and optionally \code{isotropic=TRUE} for a distribution that is symmetric in \code{x} and \code{y}.

\code{tau} is an ordered array of autocorrelation timescales. 
If \code{length(tau)==0}, then an uncorrelated bi-variate Gaussian model is fit to the data.
If \code{length(tau)==1}, then an Ornstein-Uhlenbeck (OU) model (Brownian motion restricted to a finite home range) is fit the data, where \code{tau} is the position autocorrelation timescale. \code{tau=Inf} then yields Brownian motion (BM).
If \code{length(tau)==2}, then the OUF model (continuous-velocity motion restricted to a finite home range) is fit to the data, where \code{tau[1]} is again the position autocorrelation timescale and \code{tau[2]} is the velocity autocorrelation timescale. \code{tau[1]=Inf} then yields integrated Ornstein-Uhlenbeck (IOU) motion, which is a spatially unrestricted continuous-velocity process.
%\code{length(tau)>2} and central place foraging models are not yet implemented.

%To consider endlessly diffusing processes like Brownian motion and integrated OU motion, one can simply set the last element of \code{tau} to be much larger than the period of the data.

More models will be implemented in the future.

%In the absence of significant telemetry errors it is recommended that optional \code{model} components \code{sigma}, which can be a covariance matrix or variance scalar, and \code{mu}, which is the mean vector, are left absent so that they will be analytically profiled. (TELEMETRY ERROR NOT YET IMPLEMENTED)
}

\value{
\code{ctmm.fit} returns the maximum likelihood \code{ctmm} movement-model object with the following components:
\describe{
\item{\code{AICc}}{The corrected Akaike information criterion, under the assumption of a linear model. Not all model parameters here are linear, but the linear AICc is still implemented rather than AIC because AICc is known to behave well for analogous discrete-time models. Also, following the arguments of Burnham & Anderson (2002), AICc should be an improvement here over AIC because the model parameters are not extremely nonlinear, as would be the case with hinge features, etc..}
\item{\code{loglike}}{The log-likelihood.}
\item{\code{isotropic}}{A Boolean denoting whether or not the animal's covariance is circular or elliptical.}
\item{\code{tau}}{The maximum likelihood \code{tau} estimate.}
\item{\code{sigma}}{The maximum likelihood variance/covariance estimate. For the endlessly diffusing BM and IOU processes, this becomes the diffusion rate estimate.}
\item{\code{mu}}{The maximum likelihood stationary mean vector estimate.}
\item{\code{COV.mu}}{The covariance matrix of the vector \code{mu}, assuming that the point estimate \code{sigma} is good.}
\item{\code{DOF.mu}}{The effective number of degrees of freedom in the estimate of \code{mu}, assuming that the point estimate of \code{tau} is good. This can be much smaller than \code{length(data$t)} if the data are autocorrelated.}
\item{\code{COV.tau}}{Covariance of the vector \code{c(tau,GM.sigma)}, as derived from \code{hessian}, where \code{GM.sigma=sqrt(det(sigma))} is the geometric mean of the covariance matrix \code{sigma}. The animal's home-range area is proportional to \code{GM.sigma}. Typically, sequential components of \code{tau} are slightly correlated, while the last component of \code{tau} is extremely correlated to \code{sigma}, here summarized by its geometric mean.

\code{COV.tau} assumes that \code{tau>0}, which may not be the case if OU motion is more supported than OUF motion. 
}
}
}

\references{
C. H. Fleming, J. M. Calabrese, T. Mueller, K.A. Olson, P. Leimgruber, and W. F. Fagan. (2014). From fine-scale foraging to home ranges: A semi-variance approach to identifying movement modes across spatiotemporal scales. \href{http://www.jstor.org/discover/10.1086/675504}{The American Naturalist, 183(5), E154-E167.}

C. H. Fleming and Y. Subasi and J. M. Calabrese. (2015). A maximum-entropy description of animal movement. \href{http://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.032107}{Physical Review E, 91, 032107.}
}


\author{ Chris H. Fleming. }

\note{ Good starting guesses for the parameters can be obtained from \code{variogram.fit}.

The AICs of endlessly diffusing models like BM and IOU cannot be compared to the AICs of range resident models like bivariate Gaussian, OU, and OUF. Endlessly diffusing models have to be conditioned off of an initial state, which we derive in \code{ctmm} by taking the large range limit of a range-restricted process. I.e., BM is the limit OU(\code{Inf}) and IOU(\code{tau}) is the limit OUF(\code{Inf},\code{tau}). Using comparable likelihood functions would give up statistical efficiency and the objective prior.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{variogram.fit}}, \code{\link{summary.ctmm}}, \code{\link{ctmm.loglike}}, \code{\link{optim}}. }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)
cilla <- buffalo[[1]]

# Fit a continuous-velocity model with tau ~ c(10 days,1 hour)
# also see help(variogram.fit)
m2 <- ctmm(tau=c(10*24*60^2,60^2))
M2 <- ctmm.fit(cilla,m2)

# some human-readable information
summary(M2)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
