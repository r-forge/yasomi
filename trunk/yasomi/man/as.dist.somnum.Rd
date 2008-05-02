\name{as.dist.somnum}
\alias{as.dist.somnum}
\alias{as.matrix.somnum}
\title{Compute all the distances between prototypes in a fitted Self-Organising Map}
\description{
Computes all the pairwaise distances between the prototypes of a fitted
Self-Organising Map (SOM). 
}
\usage{
\method{as.dist}{somnum}(x, FUN = NULL)
\method{as.matrix}{somnum}(x, \dots)
}
\arguments{
  \item{x}{an object of class \code{"somnum"}}
  \item{\dots}{not used}
  \item{FUN}{not used}
}
\details{
Contrarily to the \code{\link{prototype.distances}} function, those
methods compute all the pairwise distances between the prototypes of the
SOM in the original space. This is useful in order to display those
prototypes via a nonlinear projection method such the Sammon mapping
(\code{\link{sammon}}). 
}
\value{
\code{as.matrix} returns a square matrix and \code{as.dist} an object
  of class \code{"dist"}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Fabrice Rossi}
%\note{
%}
\seealso{\code{\link{prototype.distances}} for a quite different
  approach based on path calculation in the prior structure.}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a small hexagonal grid
sg <- somgrid(xdim=8,ydim=8,topo="hex")

# fit the SOM
som <- batchsom(data,sg)

# display the prototypes via sammon mapping
library(MASS)
proto.sam <- sammon(as.dist(som))
plot(proto.sam$points)

# overlay the prior structure
lines(grid2lines(som,proto.sam$points))
}
\keyword{cluster}

