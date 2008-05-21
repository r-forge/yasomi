\name{plot.somtune}
\alias{plot.somtune}
\title{Plot SOM parameter tuning results}
\description{
Visualise the results of parameter tuning for the SOM
}
\usage{
\method{plot}{somtune}(x, relative = TRUE, ...)
}
\arguments{
  \item{x}{an object of class \code{somtune}}
  \item{relative}{specifies whether the best value of the criterion
    should be used as a baseline in the display (the default behaviour)
    or not.}
  \item{\dots}{additional parameters transmitted to the underlying
    \code{\link{barplot}} function}
}
\details{
This function produces a barplot with the quantisation error of all the
tested SOMs as well as an additional bar for an another error criterion
if the quantisation error was not used to choose the best SOM. 
}
%\value{
%}
%\references{}
\author{Fabrice Rossi}
\section{Warning}{This function works only when one parameter is tuned
  by \code{\link{som.tune}}. It won't work in particular if random
  initialisation is used several times together with a parameter tuning.}
\seealso{\code{\link{som.tune}} and \code{\link{som.tunecontrol}}}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a medium hexagonal grid
sg <- somgrid(xdim=15,ydim=15,topo="hex")

# choose a good SOM via Kaski and Lagus' error measure
tune.results <-
  som.tune(data,sg,som.tunecontrol(sg,criterion=error.kaskilagus))

# plot errors
plot(tune.results)
}
\keyword{hplot}
\keyword{cluster}
