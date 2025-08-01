\name{ctmm-FAQ}
\alias{ctmm-FAQ}
\alias{ctmm-faq}
\encoding{UTF-8}

\description{Frequently asked questions for the \code{ctmm} package.}

\title{ctmm FAQ}

\details{
\bold{General recommendations}
\enumerate{
\item{Work through the vignettes \code{vignette("variogram")} and \code{vignette("akde")}. Also, see the \code{help} file for the method of interest, and its \code{example}.}
\item{Do not save workspaces between sessions. They will become corrupted over time. In RStudio, go to \code{Tools: Global Options: Workspace}, uncheck \code{Restore} and set \code{Save} to \code{Never}.}
\item{If RStudio is crashing frequently in Windows (or your display driver is crashing), try setting the rendering engine to \code{Software} under \code{Tools : Global Options : General : Advanced : Rendering Engine}.}
\item{Never edit or save your CSV in Microsoft Excel. The dates will be reformatted incorrectly and inconsistently.}
\item{If using Windows, make sure to have the suggested version of ``Rtools'' installed. If using MacOS, make sure to have ``Xcode'' installed. If using Ubuntu, make sure to have ``build-essential'' installed. Otherwise, you can sometimes run into problems when trying to update packages.}
\item{Upgrade \code{R} to the latest version and update all of your packages.}
\item{The development build can be installed via \code{remotes::install_github("ctmm-initiative/ctmm")}.}
\item{Sometimes installing from Github can silently fail to overwrite old files, requiring the package to be manually uninstalled, and then re-installed after restarting.}
\item{The \href{https://groups.google.com/g/ctmm-user}{\code{ctmm} user's group} is a good place to find and ask for help.}
\item{Bug reports and feature requests can be raised at the \href{https://github.com/ctmm-initiative/ctmm}{Github project page}.}
}

\bold{Help installing packages on Linux}

These are the packages I needed in Ubuntu:

\code{sudo apt install cmake ffmpeg fftw3 libfftw3-dev libgdal-dev libgeos-dev libgit2-dev libgmp-dev libgsl-dev libmpfr-dev libproj-dev libnode-dev libudunits2-dev r-base-core}

\bold{\code{as.telemetry} reports abnormal sampling intervals and speeds}

Make sure that you have the correct \code{timezone} and \code{timeformat} arguments specified. Also, see \code{\link{outlie}}.

\bold{rdb database corruption, "could not find function", "cannot coerce class", and other weird errors}

\code{R} might not have installed or loaded the package correctly---e.g., some files may have failed to overwrite previous versions---or the workspace/session might be corrupted. Uninstall \code{ctmm}, restart \code{R} without saving the workspace/session, and install \code{ctmm} again.

\bold{Infinite recursion and stack overflow errors}

\code{ctmm} has no recursive functions, so I am not exactly sure what causes this error, but it only occurs with certain versions of \code{R} on certain computer architectures. There are several solutions that have worked for people, including restarting \code{R} in a fresh session and updating their software. Alternatively:
\enumerate{
\item{Reboot your computer.}
\item{Increase the allowed number of nested expressions within \code{R} via \code{options(expressions=10000)} or some other large number.}
\item{Try a different computer.}
}

\bold{\code{plot} complains about the datatype or has weird errors}

Namespace collision sometimes occurs between \code{raster}, \code{sp}, \code{move}, and \code{ctmm}. Either restart \code{R} and only load the \code{ctmm} package, or run \code{ctmm::plot} instead of \code{plot}.

\bold{North is no longer up after importing data}

The default projection in \code{ctmm} does not preserve the direction of North, but better preserves distances for elongated distributions. See the \code{projection} argument in \code{\link{as.telemetry}} and the example in \code{\link[ctmm]{projection}}. The \code{\link{compass}} function is also useful for pointing north.

\bold{\code{projection} complains about the datatype and fails}

Namespace collision can occur between \code{raster} and \code{ctmm}. Either restart \code{R} and only load the \code{ctmm} package, or run \code{ctmm::projection} instead of \code{projection}.

\bold{\code{\link{ctmm.guess}} has no save button}

Maximize the plot window and/or increase your screen resolution.

\bold{\code{manipulate} panel does not popup in \code{\link{ctmm.guess}} or \code{\link{zoom,variogram-method}}}

Click the gear icon in the upper-left corner of the plot window.

\bold{Gear icon missing in \code{\link{ctmm.guess}} or \code{\link{zoom,variogram-method}}}

Recent versions of \code{manipulate} and/or RStudio seem to have some issues. Sometimes the gear icon does not render unless you re-run the function 2-5 times.

\bold{\code{manipulate::isAvailable} is not found}

You probably have an outdated copy of the \code{manipulate} package installed. Update \code{R} to the latest version and then update all of your packages. This seems to happen frequently with the MacOS release of \code{R}.
}


\author{C. H. Fleming}
