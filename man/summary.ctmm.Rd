\name{summary.ctmm}
\alias{summary.ctmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize a continuous-time movement model
}
\description{ This function returns a list of biologically interesting parameters in human readable format, as derived from a continuous-time movement model.}
\usage{
\S3method{summary}{ctmm}(object,alpha=0.05,alpha.HR=0.05,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{ctmm} movement-model object from the output of \code{ctmm.fit}. }
  \item{alpha}{ Significance level for parameter estimates. }
  \item{alpha.HR}{ Significance level for the Gaussian home-range area. }
  \item{...}{Unused options.}
}

%\details{}

\value{
A \code{matrix} with low, maximum likelihood, and high estimates for the following possible parameters:
\describe{
\item{\code{tau}}{The autocorrelation timescales.}
\item{\code{area}}{The Gaussian home-range area, where the point estimate has a significance level of \code{alpha.HR}. I.e., the core home range is where the animal is located 50\% of the time with \code{alpha.HR=0.50}.
This point estimate itself is subject to uncertainty, and is given confidence intervals derived from \code{alpha}.}
\item{\code{speed}}{The Gaussian root-mean-square (RMS) velocity, which is a convenient measure of average speed.}
}
}

%\references{}

\author{ C. H. Fleming. }

\note{ Confidence intervals on the autocorrelation timescales assume they are sufficiently greater than zero and less than infinity.  }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{ctmm.fit}}. }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)

# Extract movement data for a single animal
cilla <- buffalo[[1]]

# Find the best OU movement model
# also see help(variogram.fit)
m1 <- ctmm(tau=60*60*24*10)
M1 <- ctmm.fit(cilla,m1)

# Tell us something interpretable
summary(M1)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
