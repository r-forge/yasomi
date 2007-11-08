\name{somgrid}
\alias{somgrid}
\alias{print.somgrid}
\alias{summary.somgrid}
\title{Create a prior structure for a self-organising map}
\description{Create a prior structure (a.k.a., a grid) for a
  self-organising map.  In essence this is simple a list of points in two
  dimensions. 
}
\usage{
somgrid(xdim, ydim, topo = c("rectangular", "hexagonal"))
}
\arguments{
  \item{xdim}{the horizontal dimension of the grid}
  \item{ydim}{the vertical dimension of the grid}
  \item{topo}{the topology of the grid}
}
\value{
  An object of class \code{"somgrid"}, a list with components
  \item{pts}{a two-column matrix giving locations for the points of the
    prior structure}
  \item{dist}{a matrix containing the distances between all pairs of
    points in the grid}
  \item{diam}{the diameter of the grid, i.e., the largest distance
    between to points in the grid}
  \item{xdim,ydim,topo,basesize}{as in the arguments to \code{somgrid}}

}
\author{Fabrice Rossi}
\examples{
# display a rectangular grid
plot(somgrid(10,7))

# and a hexagonal one
plot(somgrid(9,11,topo="hexa"))
}
\keyword{cluster}

