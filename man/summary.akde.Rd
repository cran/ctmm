\name{summary.akde}
\alias{summary.akde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize an autocorrelated kernel density estimate
}
\description{ This function returns a list of biologically interesting parameters in human readable format, as derived from an autocorrelated kernel density estimate.}
\usage{
\S3method{summary}{akde}(object,alpha.HR=0.05,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An \code{akde} autocorrelated kernel-density estimate from the output of \code{akde}. }
  \item{alpha.HR}{ Significance level for the home-range area. }
  \item{...}{Unused options.}
}

%\details{}

\value{
A \code{matrix} with low, maximum likelihood, and high estimates for the following parameters:
\describe{
\item{\code{area}}{The home-range area with fraction of inclusion \code{1-alpha.HR}. E.g., the 50\% core home range is estimated with \code{alpha.HR=0.50}. The confidence intervals reported here reflect the uncertainty in optimal smoothing, and are set by the \code{alpha} option of \code{akde()} when generating the \code{akde} object. In the case of IID data, where the uncertainty in optimal bandwidth is much smaller than the maximum likelihood optimal bandwidth (depicted by gridlines in \code{plot.telemetry()}), then the confidence intervals reported here are not particularly informative.}
}
}

\references{
C. H. Fleming and W. F. Fagan and T. Mueller and K. A. Olson and P. Leimgruber and J. M. Calabrese (2015).
Rigorous home-range estimation with movement data: A new autocorrelated kernel-density estimator.
\href{http://www.esajournals.org/doi/abs/10.1890/14-2010.1}{Ecology, 96(5), 1182-1188}.
}

\author{ C. H. Fleming. }

%\note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}. }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)

# Extract movement data for a single animal
cilla <- buffalo[[1]]

# Find the best OU movement model
# also see help(variogram.fit)
m1 <- ctmm(tau=6*24*60^2)
M1 <- ctmm.fit(cilla,m1)

# Estimate and summarize the AKDE
KD1 <- akde(cilla,M1)
summary(KD1)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
