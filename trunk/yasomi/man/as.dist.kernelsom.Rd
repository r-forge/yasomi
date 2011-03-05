\name{as.dist.kernelsom}
\alias{as.dist.kernelsom}
\alias{as.matrix.kernelsom}
\title{Compute all the distances between prototypes in a fitted Self-Organising Map}
\description{
Computes all the pairwaise distances between the prototypes of a fitted
Kernel Self-Organising Map (SOM). 
}
\usage{
\method{as.dist}{kernelsom}(x, FUN = NULL)
\method{as.matrix}{kernelsom}(x, \dots)
}
\arguments{
  \item{x}{an object of class \code{"kernelsom"}}
  \item{\dots}{not used}
  \item{FUN}{not used}
}
\details{
Contrarily to the \code{\link{prototype.distances}} function, those
methods compute all the pairwise distances between the prototypes of the
SOM in the kernel space. This is useful in order to display those
prototypes via a nonlinear projection method such the Sammon mapping
(\code{\link{sammon}}).

The distances are computed in the kernel space. If the kernel used to
fit the SOM is not positive, negative values may be obtained. When such
values are generated, they are replaced by \code{NA}.
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
K <- as.kernelmatrix(tcrossprod(scale(iris[1:4])))

# a small hexagonal grid
sg <- somgrid(xdim=8,ydim=8,topo="hex")

# fit the SOM (random initialisation)
som <- batchsom(K,sg,init="random",method="cluster")

# display the prototypes via sammon mapping
library(MASS)
proto.sam <- sammon(as.dist(som))
plot(proto.sam$points)

# overlay the prior structure
lines(grid2lines(som,proto.sam$points))
}
\keyword{cluster}

