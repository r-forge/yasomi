\name{batchkmeans.kernelmatrix}
\alias{batchkmeans.kernelmatrix}
\title{Kernel K-means Clustering}
\description{
Performs K-means clustering via an implicit mapping to a kernel induced feature space.
}
\usage{
\method{batchkmeans}{kernelmatrix}(data, ncenters, init = c("prototypes", "random", "cluster"),
prototypes, weights, max.iter, verbose = FALSE, keepdata = TRUE, ...)
}
\arguments{
  \item{data}{the data to cluster represented by a
    kernel matrix (an object of class \code{"kernelmatrix"}, see \code{\link{as.kernelmatrix}}).}
  \item{ncenters}{the number of clusters}
  \item{init}{the initialisation method (see \code{\link{batchkmeans}} for details)}
  \item{prototypes}{Initial values for the
    prototypes. It contains linear coefficients that describe the
    prototypes as linear combination of the initial data
    points mapped in the kernel induced feature space. It has therefore
    one row for each prototype (\code{ncenters}) and has many columns as data
    points. If missing, initial prototypes are chosen via
    the method specified by the \code{init} parameter (see details)}
  \item{weights}{optional weights for the data points}
  \item{max.iter}{maximal number of iterations of the algorithm}
  \item{verbose}{switch for tracing the clustering process}
  \item{keepdata}{if \code{TRUE}, the original data are returned as
    part of the result object}
  \item{\dots}{not used}
}
\details{
  This method implements Kernel K-means: the kernel is used to define an
  implicit mapping of the original data into a kernel induced feature
  space where a standard K-means is applied. The mapping is not
  explicitely carried out. 
}
\value{
  An object of class \code{"batchkmeans"} and of class
  \code{"kernelbatchkmeans"}, a list with components 
  including
  \item{prototypes}{a matrix containing the coordinates of the
    prototypes: each row of the matrix sums to one and contains the
    coefficients of a linear combination of the original observations in
  the kernel induced feature space}
  \item{classif}{a vector of integer indicating to which cluster each
    observation has been assigned}
  \item{errors}{a vector containing the evolution of the quantisation
    error during the fitting process}
  \item{data}{the original data if the function is called with
    \code{keepdata = TRUE}}
  \item{weights}{the weights of the data points if the function is called with
    \code{keepdata = TRUE} and if the \code{weights} is given}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Fabrice Rossi}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{See \code{\link{batchsom.kernelmatrix}} for Self-Oganising Map which
  provides both clustering and visualisation.}
\examples{
## inspired by the kmeans example in the stats package
require(graphics)

## 2 dimensional dataset
X <- cbind(rnorm(200,mean=2,sd=0.35),rnorm(200,mean=-1,sd=0.35))
Y <- cbind(runif(200,min=-1.5,max=-0.75),runif(200,min=0,max=0.5))
Z <- cbind(rnorm(200,sd=0.15),rnorm(200,sd=0.5))
M <- matrix(c(sin(pi/4),cos(pi/4),-cos(pi/4),sin(pi/4)),ncol=2)
U <- rbind(X,Y,Z\%*\%M+c(rep(0.25,200),rep(-0.5,200)))
U <- scale(U)

K <- as.kernelmatrix(tcrossprod(U))

km <- batchkmeans(K,3)

plot(U,col=km$classif,cex=0.5)

## map back the prototypes to the original space
points(km$prototypes\%*\%U,col=1:3,pch=20,cex=2)
}
\keyword{cluster}
