\name{sominit.pca.default}
\alias{sominit.pca.default}
\title{Initialise the prototypes of a SOM with PCA}
\description{
Initialise the prototypes of a Self-Organising Map with Principal Component
Analysis. The prototypes are regulary positioned (according to the prior
structure) in the subspace spanned by the two first principal
components. 
}
\usage{
\method{sominit.pca}{default}(data, somgrid, \dots)
}
\arguments{
  \item{data}{the data to which the SOM will be fitted, a matrix or data frame of observations (which should be scaled)}
  \item{somgrid}{a \code{somgrid} object}
  \item{\dots}{not used}
}
\value{
  A matrix containing appropriate initial prototypes
}
\author{Fabrice Rossi}
\seealso{\code{\link{somgrid}} for specifying the prior structure and
  \code{\link{sominit.random}} for random based initialisations.}
\examples{
X <- cbind(rnorm(500),rnorm(500))

sg <- somgrid(xdim=7,ydim=7,topo="rect")

proto <- sominit.pca(X,sg)

plot(X,pch="+",col="red",xlim=range(X[,1],proto[,1]),
     ylim=range(X[,2],proto[,2]))
points(proto,pch=20)

}
\keyword{cluster}
