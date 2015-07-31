\name{periodogram}
\alias{periodogram}
\alias{plot.periodogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the Lomb-Scargle periodogram of animal-tracking data
}
\description{
This function calculates the exact, isotropic Lomb-Scargle periodogram (LSP, Scargle, 1982) from a telemetry object.
}
\usage{
periodogram(data,T=NULL,dt=NULL,resolution=1)

\S3method{plot}{periodogram}(x,diagnostic=FALSE,col="black",transparency=0.25,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ \code{telemetry} data object of the 2D timeseries data. }
  \item{T}{Sampling period for frequency resolution.}
  \item{dt}{Sampling interval for frequency cutoff.}
  \item{resolution}{Multiplier to inflate the frequency resolution.}
  \item{x}{Output object of \code{periodogram}.}
  \item{diagnostic}{Plot the sampling schedule's periodogram to check for spurious periodicities.}
  \item{col}{Color of periodogram.}
  \item{transparency}{Adds transparency to clustered data if greater than zero. Should be less than one.}
  \item{...}{Optional arguments fed to \code{\link{plot}}.}
}
\details{
If no \code{dt} is specified, the median sampling interval is used. This is typically a good assumption for most data, even when there are gaps and this choice corresponds to the discrete Fourier transform (DFT) periodogram for evenly-sampled data. If no \code{T} is specified, the full sampling period is used. 

The frequency grid interval is finally given by \code{1/(resolution*2*(T+dt))} and the frequency cutoff is given by \code{1/(2*dt)}, both in accordance with the DFT periodogram. Increasing \code{resolution} will make for a smooth periodogram, but sequential frequencies will be highly correlated.
}

\value{
Returns a periodogram object (class \code{periodogram}) which is a dataframe containing the frequency, \code{f}, the Lomb-Scargle periodogram at that frequency, \code{LSP}, and the approximate number of degrees of freedom associated with that value, \code{DOF}, with which its confidence intervals can be estimated.
}
\references{
J. D. Scargle. (1952). Studies in astronomical time-series analysis. II. Statistical aspects of spectral analysis of unevenly-sampled data. The Astrophysical Journal, 263, 835-853.
}

\author{ Chris H. Fleming. }

\note{ The LSP is totally inappropriate if you in any way alter the sampling rate within the dataset. Stick with variograms in that case. There is a \code{diagnostic} option in \code{plot.periodogram} that can check for spurious periodicities that result from an autocorrelated sampling schedule. This plot will not contain any periodicities if the LSP is appropriate.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{ \code{\link{plot.periodogram}}. }

\examples{
#Load package and data
library(ctmm)
data(buffalo)

#Extract movement data for a single animal
cilla <- buffalo[[1]]

#Calculate periodogram
per <- periodogram(cilla)

#Plot the periodogram
plot(per)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
