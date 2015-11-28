\name{variogram.fit}
\alias{variogram.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visually fit a movement model to a variogram
}
\description{This function plots a \code{variogram} object overlayed with a continuous-time movement model guesstimated from the variogram's shape. Sliders are given to adjust the parameter guesstimates and the result can be saved to a global variable. The intention of this function is to facilitate good starting guesses for \code{ctmm.fit}.}

\usage{
variogram.fit(variogram,range=TRUE,CTMM=NULL,name="variogram.fit.model",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variogram}{ A \code{variogram} object from the output of \code{variogram}. }
  \item{range}{Whether or not processes with a finite home range are considered.}
  \item{CTMM}{Optional initial guesstimate of the model parameters, in \code{ctmm} object format. The slider range is triple this value.}
  \item{name}{Name of the global variable to store the guesstimate in.}
  \item{...}{Optional parameters passed to \code{plot.variogram}.}
}


\details{
When \code{range=TRUE}, \code{sigma} is the asymptote of the variogram and \code{tau} is an array of autocorrelation timescales. The position timescale is roughly the time lag it takes of the variogram to reach 63\% of its asymptote. The velocity autocorrelation timescale visually corresponds to width of the concave bowl shape at the beginning of the variogram. When \code{range=FALSE}, \code{sigma} is the asymptotic slope of the variogram and only the velocity timescale is finite.

By default, parameter values are estimated from the shape of the variogram. If this fails, the \code{CTMM} option can provide alternative initial guesstimates.
}

%\value{}

%\references{}

\author{ C. H. Fleming. }

%\note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{ctmm.fit}}, \code{\link{plot.variogram}}, \code{\link{variogram}}. }

\examples{
#Load package and data
library(ctmm)
data(buffalo)

#Extract movement data for a single animal
cilla <- buffalo[[1]]

#Calculate variogram
SVF <- variogram(cilla)

# generate a visual fit of the variogram (requires RStudio)
\donttest{variogram.fit(SVF)}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
