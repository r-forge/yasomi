\name{as.kernelmatrix}
\alias{as.kernelmatrix}
\title{Turn a standard matrix into a Kernel matrix}
\description{
This function simply adds the \code{"kernelmatrix"} class to a standard
matrix for subsequent use for a Kernel Self-Organising Map. Optionally
performs some sanity check on the matrix to verify that it represents a
non negative kernel matrix.
}
\usage{
as.kernelmatrix(data, ...)
}
\arguments{
  \item{data}{an object of class \code{"matrix"}}
  \item{\dots}{not used}
}
\details{
The current version of the function only checks for symmetry of the
\code{data}. Future versions will provide optional non negativity tests.
}
\value{
  An object of class \code{"kernelmatrix"} which is in fact a
  matrix.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Fabrice Rossi}
%\note{ ~~further notes~~ 
%}
\seealso{See \code{\link{batchsom.kernelmatrix}} for the main use of
  Kernel matrices}
\examples{
## the linear kernel applied to the Iris dataset
data(iris)
# scaling and kernel computation
data <- as.kernelmatrix(tcrossprod(scale(iris[1:4])))
}
\keyword{cluster}

