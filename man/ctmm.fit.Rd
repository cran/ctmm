\name{ctmm}
\alias{ctmm}
\alias{ctmm.fit}
\alias{ctmm.select}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specify, fit, and select continuous-time movement models
}
\description{These functions allow one to propose hypothethical movement models (with initial estimates), fit those models to the data, and select among those models via an information criteria.
The fitting functions wrap around \code{optim} and \code{ctmm.loglike} to maximize the likelihood function of continuous-time movement models described in Fleming et al (2014) and Fleming et al (2015), given 2D animal tracking data.
}
\usage{
ctmm(tau=NULL,isotropic=FALSE,range=TRUE,circle=FALSE,CPF=FALSE,error=FALSE,...)

ctmm.fit(data,CTMM=ctmm(),debias=TRUE,control=list(maxit=.Machine$integer.max),...)

ctmm.select(data,CTMM,verbose=FALSE,IC="AICc",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tau}{Array of autocorrelation timescales explained below.}
  \item{isotropic}{A Boolean denoting whether or not the animal's covariance is circular or elliptical.}
  \item{range}{A Boolean denoting whether or not the movement model has a finite range.}
  \item{circle}{The period it takes the animal to stochastically circle its mean location.}
  \item{CPF}{A Boolean denoting whether to consider an oscillatory CPF model.}
  \item{error}{A Boolean denoting whether or not to use annotated telemetry error estimates or an estimate of the error's variance if the data are not annotated with error estimates.}
  \item{data}{ The 2D timeseries data represented as a \code{telemetry} object. }
  \item{CTMM}{ A \code{ctmm} movement-model object containing the initial parameter guesses conforming to the basic structure of the model hypothesis. \code{ctmm.select} can accept a list of such objects.}
  \item{debias}{ A Boolean denoting whether or not to remove the lowest-order bias from the covariance \code{sigma} estimate. }
  \item{control}{An optional argument to be passed to \code{\link{optim}}, but with \code{parscale} overwritten with reasonable defaults.}
  \item{...}{Further arguments passed to \code{\link{optim}}.}
  \item{verbose}{Return every fit candidate model if \code{TRUE}, else return only the selected model.}
  \item{IC}{Information criteria used. Only \code{"AICc"} is currently supported.}
}

\details{The initial \code{ctmm} parameter guess can be the output of \code{ctmm.guess}, \code{variogram.fit} or the function \code{ctmm(...)} with the argument \code{tau} explained below and additonal model options described in \code{vignette("ctmm")}.

By default, \code{tau} is an ordered array of autocorrelation timescales. 
If \code{length(tau)==0}, then an IID bi-variate Gaussian model is fit to the data.
If \code{length(tau)==1}, then an Ornstein-Uhlenbeck (OU) model (Brownian motion restricted to a finite home range) is fit the data, where \code{tau} is the position autocorrelation timescale. \code{tau=Inf} then yields Brownian motion (BM).
If \code{length(tau)==2}, then the OUF model (continuous-velocity motion restricted to a finite home range) is fit to the data, where \code{tau[1]} is again the position autocorrelation timescale and \code{tau[2]} is the velocity autocorrelation timescale. \code{tau[1]=Inf} then yields integrated Ornstein-Uhlenbeck (IOU) motion, which is a spatially unrestricted continuous-velocity process.
%\code{length(tau)>2} and central place foraging models are not yet implemented.

If \code{CPF=TRUE}, then an oscillatory central place foraging model is assumed. In this case \code{tau} must be an length-2 array containaing the foraging period in \code{tau[1]} (e.g., 1 day \code{24*60^2} seconds) and the characteristic timescale overwhich similarity between foraging excursions persists in \code{tau[2]}.

\code{debias=TRUE} causes the maximum likelihood estimate (MLE) of the covariance \code{sigma} to be corrected by a prefactor of \code{n/(n-k)}, where \code{n} is the number of data points and \code{k} is the number of mean parameters. If \code{error=FALSE} and other autocorrelation parameters such as \code{tau} are exact, then this correction is exact. If \code{error=FALSE} and and the other autocorrelation parameters are simply their MLE values, then this correction removes the lowest-order bias in \code{sigma}. If \code{error=TRUE}, then this is an undercorrection that is better than doing nothing. Residual maximum likelihood estimation (REML) could be more exact for \code{error=TRUE}, buth there is a potential tradeoff between bias and variablity.
}

\value{
\code{ctmm.fit} returns the maximum likelihood \code{ctmm} movement-model object with all of the components of \code{CTMM} plust the following:
\describe{
\item{\code{AICc}}{The corrected Akaike information criterion, under the assumption of a linear model. Not all model parameters here are linear, but the linear AICc is still implemented rather than AIC because AICc is known to behave well for analogous discrete-time models. Also, following the arguments of Burnham & Anderson (2002), AICc should be an improvement over AIC here because the model parameters are not extremely nonlinear, etc..}
\item{\code{loglike}}{The log-likelihood.}
\item{\code{sigma}}{The maximum likelihood variance/covariance estimate (possibly debiased). For the endlessly diffusing BM and IOU processes, this is instead the diffusion rate estimate.}
\item{\code{mu}}{The maximum likelihood stationary mean vector estimate.}
\item{\code{COV.mu}}{The covariance matrix of the vector \code{mu}, assuming that the point estimate \code{sigma} is good.}
\item{\code{DOF.mu}}{The effective number of degrees of freedom in the estimate of \code{mu}, assuming that the point estimate of \code{tau} is good. This can be much smaller than \code{length(data$t)} if the data are autocorrelated.}
\item{\code{COV}}{Covariance of the parameter vector \code{c(area,tau,circle)}, as derived from \code{hessian}, where \code{area=sqrt(det(sigma))} is the geometric mean of the covariance matrix \code{sigma} or the standard area of the ellipse. The animal's home-range area is proportional to \code{area}. Typically, \code{sigma} is extremely correlated to \code{tau[1]}, and sequential components of \code{tau} are slightly correlated.
}
}
}

\references{
C. H. Fleming, J. M. Calabrese, T. Mueller, K.A. Olson, P. Leimgruber, and W. F. Fagan. (2014). From fine-scale foraging to home ranges: A semi-variance approach to identifying movement modes across spatiotemporal scales. \href{http://www.jstor.org/discover/10.1086/675504}{The American Naturalist, 183(5), E154-E167.}

C. H. Fleming and Y. Subasi and J. M. Calabrese. (2015). A maximum-entropy description of animal movement. \href{http://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.032107}{Physical Review E, 91, 032107.}
}


\author{ C. H. Fleming and G. Péron. }

\note{ Good starting guesses for the timescale and area parameters can be obtained from \code{variogram.fit}.
Good starting estimates for \code{circle} can be obtained from \code{ctmm.guess}.

Fitting without telemetry errors involves fewer parameters and is much faster than fitting with telemetry errors, as the variance/covariance can be profiled analytically without numerical optimization.
Therefore, fitting without telemetry errors, if possible, can provide improved initial estimates of \code{tau} and \code{sigma} for later fitting with telemetry errors.

The default \code{\link{optim}} method here is \code{"Nelder-Mead"} with the largest allowable number of maximum iterations. Annecdotally, the faster \code{"BFGS"} and \code{"L-BFGS-B"} methods often perform poorly on these types of problems.

The AICs of endlessly diffusing models like BM and IOU cannot be easily compared to the AICs of range resident models like bivariate Gaussian, OU, and OUF. Endlessly diffusing models have to be conditioned off of an initial state, which we derive in \code{ctmm} by taking the large range limit of a range-restricted process. I.e., BM is the limit OU(\code{Inf}) and IOU(\code{tau}) is the limit OUF(\code{Inf},\code{tau}). Using comparable likelihood functions would give up statistical efficiency and the objective prior.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{ctmm.loglike}}, \code{\link{ctmm.guess}}, \code{\link{optim}}, \code{\link{summary.ctmm}}, \code{\link{variogram.fit}}. }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)
cilla <- buffalo[[1]]

# Fit a continuous-velocity model with tau ~ c(10 days,1 hour)
# also see help(variogram.fit)
GUESS <- ctmm(tau=c(10*24*60^2,60^2))
FIT <- ctmm.fit(cilla,GUESS)

# some human-readable information
summary(FIT)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
