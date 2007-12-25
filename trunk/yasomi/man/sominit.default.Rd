\name{sominit.default}
\alias{sominit.default}
\title{Initialize the prototypes of a SOM with PCA}
\description{
Initialize the prototypes of a Self-Organising Map with Principal Component
Analysis.  The prototypes are regulary positioned (according to the prior
structure) in the subspace spanned by the two first principal components.
}
\usage{
\method{sominit}{default}(data, somgrid,...)
}
\arguments{
  \item{somgrid}{a \code{somgrid} object}
  \item{data}{a matrix or data frame of observations (should be scaled)}
  \item{...}{not used}
}
\value{
  A matrix containing appropriate initial prototypes
}
\author{Fabrice Rossi}
\seealso{\code{\link{sominit}} \code{\link{somgrid}}}
\examples{
X <- cbind(rnorm(500),rnorm(500))

sg <- somgrid(xdim=7,ydim=7,topo="rect")

proto <- sominit(X,sg)

plot(X,pch="+",col="red",xlim=range(X[,1],proto[,1]),
     ylim=range(X[,2],proto[,2]))
points(proto,pch=20)

}
\keyword{cluster}

