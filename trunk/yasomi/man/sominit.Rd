\name{sominit}
\alias{sominit}
\title{Initialize the prototypes of a SOM}
\description{
Generic initializing method for the prototypes of a Self-Organising Map.
}
\usage{
sominit(data, somgrid,...)
}
\arguments{
  \item{data}{the data to which the SOM will be fitted. This can be,
    e.g., a matrix or data frame of observations (which should be scaled), or a
    distance matrix}
  \item{somgrid}{a \code{somgrid} object}
  \item{...}{arguments to be passed to methods (such as the actual
    method used for initialisation)}
}
\value{
  A matrix containing appropriate initial prototypes. It should be
  compatible with the SOM prior structure (i.e., it should have as many
  rows as the size of the grid) and with the data. 
}
\author{Fabrice Rossi}
\seealso{\code{\link{somgrid}} and \code{\link{sominit.default}} for the
default PCA based initialisation for matrix data}
\examples{
X <- cbind(rnorm(500),rnorm(500))

sg <- somgrid(xdim=7,ydim=7,topo="rect")

proto <- sominit(X,sg)

plot(X,pch="+",col="red",xlim=range(X[,1],proto[,1]),
     ylim=range(X[,2],proto[,2]))
points(proto,pch=20)

}
\keyword{cluster}

