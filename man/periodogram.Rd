\name{periodogram}
\alias{periodogram}
\alias{plot.periodogram}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the Lomb-Scargle periodogram of animal-tracking data
}
\description{
This function calculates isotropic Lomb-Scargle periodogram (LSP, Scargle, 1982) from a telemetry object. One of two algorithms is used. The slow \eqn{O(n^2)} algorithm vectorizes the exact relations of Scargle (1982), while the fast \eqn{O(n \log n)} algorithm uses the FFT method described in Péron et al (2016). The latter method is only exact if the data are evenly scheduled, permitting gaps, and otherwise it produces some distortion at high frequencies.
}
\usage{
periodogram(data,CTMM=NULL,T=NULL,dt=NULL,res=1,fast=NULL)

\S3method{plot}{periodogram}(x,diagnostic=FALSE,col="black",transparency=0.25,grid=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ \code{telemetry} data object or list of such objects. }
  \item{CTMM}{ An optional \code{ctmm} model object for specifying the mean. }
  \item{T}{Sampling period for frequency resolution.}
  \item{dt}{Sampling interval for frequency cutoff.}
  \item{res}{Multiplier to inflate the frequency resolution.}
  \item{fast}{Use the exact algorithm if FALSE or the FFT algorithm if TRUE. The default behavior switches from slow to fast at length(data$t)>=10000.}
  \item{x}{Output object of \code{periodogram}.}
  \item{diagnostic}{Plot the sampling schedule's periodogram to check for spurious periodicities.}
  \item{col}{Color of periodogram.}
  \item{transparency}{Adds transparency to clustered data if greater than zero. Should be less than one.}
  \item{grid}{Whether or not to plot gridlines at common periodicities.}
  \item{...}{Optional arguments fed to \code{\link{plot}}.}
}
\details{
If no \code{dt} is specified, the median sampling interval is used. This is typically a good assumption for most data, even when there are gaps and this choice corresponds to the discrete Fourier transform (DFT) periodogram for evenly-sampled data. If no \code{T} is specified, the full sampling period is used. 

The frequency grid interval is finally given by \code{1/(2*(T+dt)*res)} and the frequency cutoff is given by \code{1/(2*dt)}, both in accordance with the DFT periodogram. Increasing \code{res} beyond \code{res=1} will make for a smooth periodogram, but sequential frequencies will be highly correlated.

If a \code{ctmm} argument is provided, the ML mean will be detrended from the data prior to calculating the periodogram. Otherwise, the sample mean will be detrended.

If a list of \code{telemetry} objects are fed into \code{periodogram}, then a mean \code{periodogram} object will be returned with the default \code{T} and \code{dt} arguments selected on a worst case basis according to the method described by Péron et al (2016).
}

\value{
Returns a periodogram object (class \code{periodogram}) which is a dataframe containing the frequency, \code{f} and the Lomb-Scargle periodogram at that frequency, \code{LSP}.
}
\references{
J. D. Scargle. (1952). Studies in astronomical time-series analysis. II. Statistical aspects of spectral analysis of unevenly-sampled data. The Astrophysical Journal, 263, 835-853.

G. Péron, C. H. Fleming, J. M. Calabrese. (2016).
}

\author{ C. H. Fleming and G. Péron }

\note{ The LSP is totally inappropriate if you in any way alter the sampling rate within the dataset. Stick with variograms in that case. There is a \code{diagnostic} option in \code{plot.periodogram} that can check for spurious periodicities that result from an autocorrelated sampling schedule. This plot will not contain any periodicities if the LSP is appropriate.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{ \code{\link{plot.periodogram}}. }

\examples{
#Load package and data
library(ctmm)
data(buffalo)

#Extract movement data for a single animal
cilla <- buffalo[[1]]

#Calculate periodogram (fast=TRUE for example)
LSP <- periodogram(cilla,fast=TRUE)

#Plot the periodogram
plot(LSP)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
