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
\method{sominit.pca}{default}(data, somgrid, weights, with.princomp=FALSE, \dots)
}
\arguments{
  \item{data}{the data to which the SOM will be fitted, a matrix or data frame of observations (which should be scaled)}
  \item{somgrid}{a \code{somgrid} object}
  \item{weights}{optional weights for the data points}
  \item{with.princomp}{switch specifying whether the
    \code{\link{princomp}} should be used instead of the
    \code{\link{prcomp}} method for computing the principal components
    when no weights are given (see details)}
  \item{\dots}{not used}
}
\details{
When the optional \code{weights} are specified, the weighted covariance of the
data is computed via \code{\link{cov.wt}}. Then \code{\link{princomp}}
is used to find the two first principal components of the data. When
\code{weights} are missing, the PCA is conducted via
\code{\link{prcomp}}, expect is the function is told to use
\code{\link{princomp}} via the \code{with.princomp} parameter. As a
consequence, if
\code{with.princomp=FALSE}, the results of the function applied to
unweighted data points are likely to differ from the ones obtained on
the same data points with uniform weights. 
  }
\value{
  A list with the following components
  \item{prototype}{a matrix containing appropriate initial prototypes}
  \item{data.pca}{the results of the PCA conducted on the data via a
    call to \code{\link{prcomp}} or \code{\link{princomp}}}
}
\author{Fabrice Rossi}
\seealso{\code{\link{somgrid}} for specifying the prior structure and
  \code{\link{sominit.random}} for random based initialisations.}
\examples{
X <- cbind(rnorm(500),rnorm(500))

sg <- somgrid(xdim=7,ydim=7,topo="rect")

proto <- sominit.pca(X,sg)$prototypes

plot(X,pch="+",col="red",xlim=range(X[,1],proto[,1]),
     ylim=range(X[,2],proto[,2]))
points(proto,pch=20)

}
\keyword{cluster}
