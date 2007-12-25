\name{grid2lines}
\alias{grid2lines}
\title{Self-Organising Map grid}
\description{
Build a matrix that can be used to plot the direct neighbourhood relationship in the prior structure of a Self-Organising Map.
}
\usage{
grid2lines(som,prototypes)
}
\arguments{
  \item{som}{an object of class \code{"som"}}
  \item{prototypes}{the coordinates of the objects to link in the result
  matrix (if missing, the coordinates of the prototypes of the
  \code{som} are used).}
}
\details{
\code{grid2lines} returns a transformed version of the coordinate matrix
\code{prototypes} (or of the prototypes of a the fitted Self-Organising
Map \code{som}).  It can be used with,
e.g., \code{\link{lines}} to draw edges between prototypes that are
direct neighbours in the prior structure. 
}
\value{
a matrix with the same number of columns as the original dataset used by
the \code{som} object. 
}
\author{Fabrice Rossi}
\examples{
X <- cbind(rnorm(500),rnorm(500))

sg <- somgrid(xdim=7,ydim=7,topo="rect")

som <- batchsom(X,sg)

# display the data, the som's prototypes and the prior structure
plot(X,pch="+",col="red")
lines(grid2lines(som),type="b",pch=20)
}
\keyword{cluster}
\keyword{dplot}
