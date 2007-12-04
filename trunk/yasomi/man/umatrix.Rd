\name{umatrix}
\alias{umatrix}
\title{U-Matrix}
\description{
Plot a U-matrix for a Self-Organising Map.
}
\usage{
umatrix(som, ...)
}
\arguments{
  \item{som}{an object of class \code{"som"}}
  \item{\dots}{additional parameters transmitted to the low level plot function \code{\link{plot.sompdist}}}
}
\details{
  Many variations of Ultsch and Siemon's original U-matrix have been
  used in scientific publications. The \code{umatrix} function provides
  one of the simplest solution, based on a coloring of the glyphs used
  to represent the prior structure of the fitted SOM. With the default
  parameters, the colour of each cell encodes the mean distance between the
  corresponding prototype and its direct neighbours (using
  \code{\link{heat.colors}}). Details can be
  found in \code{\link{plot.sompdist}}. 
}
\value{
  \code{umatrix} silently returns the object of class \code{"sompdist"}
  build from the \code{som} parameters and used to construct the
  graphical representation. 
}
\references{
  Ultsch, A. and Siemon, H. P. (1990) Kohonen's self organizing feature
  maps for exploratory data analysis, in: Proceedings of International
  Neural Network Conference (INNC'90).
}
\author{Fabrice Rossi}
%\note{ ~~further notes~~ 
%}
\seealso{See \code{\link{plot.sompdist}} for details,
  \code{\link{distance.grid}} for possibly smoother plots.}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a medium hexagonal grid
sg <- somgrid(xdim=15,ydim=15,topo="hex")

# choose a good SOM via Kaski and Lagus' error measure
st <- som.tune(sg,data,som.tunecontrol(sg,,criterion=error.kaskilagus))
som <- st$best.som

# display the umatrix
umatrix(som)
}
\keyword{cluster}
\keyword{hplot}

