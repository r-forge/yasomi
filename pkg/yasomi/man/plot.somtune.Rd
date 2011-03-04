\name{plot.somtune}
\alias{plot.somtune}
\title{Plot SOM parameter tuning results}
\description{
Visualise the results of parameter tuning for the SOM
}
\usage{
\method{plot}{somtune}(x, relative = TRUE,
sqrt.quant=!x$isquant, best.color="red", xlab, ylab, yaxs, legend.text, ...)
}
\arguments{
  \item{x}{an object of class \code{somtune}}
  \item{relative}{specifies whether the best (i.e., lowest) value of the
    criterion should be used as a baseline in the display (the default
    behaviour) or not}
  \item{sqrt.quant}{specifies whether the plot should display the
    quantisation error or its square root. When \code{\link{som.tune}}
    has been used with the quantisation error as the quality measure,
    the default is to show it directly, while is square root is used in
    other cases}
  \item{best.color}{gives the colour of the border of bar associated to
    the best configuration. Use \code{NA} to suppress this effect}
  \item{xlab}{a label for the x axis (see \code{\link{barplot}})}
  \item{ylab}{a label for the y axis (see \code{\link{barplot}})}
  \item{yaxs}{the style of axis interval calculation to be used for the
    y-axis (see \code{\link{par}})}
  \item{legend.text}{a vector of text used to construct a legend for the
    plot (see \code{\link{barplot}})}
  \item{\dots}{not used}
}
\details{
This function produces a barplot with the quantisation error of all the
tested SOMs as well as an additional bar for an another error criterion
if the quantisation error was not used to choose the best SOM. Most of
the parameters have reasonable default values computed from the
\code{somtune} object. The function tries to avoid superposition of the
legend with the bars. 
}
%\value{
%}
%\references{}
\author{Fabrice Rossi}
\section{Warning}{This function works only when one parameter is tuned
  by \code{\link{som.tune}}. It will not work in particular if random
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
