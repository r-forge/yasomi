\name{componentPlane}
\alias{componentPlane}
\title{Component Planes}
\description{
Plot a component plane for a Self-Organising Map.  The prior structure
is represented by a regular grid of hexagons or squares which are
colored according to the value of the chosen coordinate of the
associated prototype.  
}
\usage{
componentPlane(som, dim = 1, ...)
}
\arguments{
  \item{som}{an object of class \code{"som"}}
  \item{dim}{the coordinate to display}
  \item{\dots}{additional parameters transmitted to the low level plot
    function \code{\link{plot.somgrid}}} 
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\references{
  Vesanto, J. (2002) \emph{Data Exploration Process Based on the
    Self--Organizing Map}, PhD thesis, Helsinki University of Technology, Espoo (Finland), Acta Polytechnica Scandinavica, Mathematics and Computing Series No.
  115.
}
\author{Fabrice Rossi}
\seealso{\code{\link{plot.somgrid}}}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a medium hexagonal grid
sg <- somgrid(xdim=15,ydim=15,topo="hex")

# choose a good SOM
som <- som.tune(sg,data)$best.som

# display the component planes
spar <- par(mfrow=c(2,2))
for(i in 1:ncol(data)) {
    componentPlane(som,i)
}
par(spar)
}
\keyword{cluster}
\keyword{hplot}
