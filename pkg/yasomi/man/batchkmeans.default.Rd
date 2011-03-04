\name{batchkmeans.default}
\alias{batchkmeans.default}
\title{K-means Clustering for vector data}
\description{
Performs K-means clustering on vector data.
}
\usage{
\method{batchkmeans}{default}(data, ncenters, init = c("prototypes", "random", "cluster"),
prototypes, weights, max.iter, verbose = FALSE, keepdata = TRUE, ...)
}
\arguments{
  \item{data}{the data to cluster. This can be,
    a matrix or data frame of observations}
  \item{ncenters}{the number of clusters}
  \item{init}{the initialisation method ((see \code{\link{batchkmeans}} for details)}
  \item{prototypes}{Initial values for the
    prototypes (it must have the same number of columns as \code{data}
    and its row number must be equal to \code{ncenters}). If missing,
    initial prototypes are chosen via 
    the method specified by the \code{init} parameter (see details)}
  \item{weights}{optional weights for the data points}
  \item{max.iter}{maximal number of iterations of the algorithm}
  \item{verbose}{switch for tracing the clustering process}
  \item{keepdata}{if \code{TRUE}, the original data are returned as
    part of the result object}
  \item{\dots}{not used}
}
\details{
  This methods implements the standard (batch) K-means clustering
  algorithm (more precisely the Lloyd-Forgy algorithm). 
}
\value{
  An object of class \code{"batchkmeans"} and of class
  \code{"batchkmeansnum"}, a list with components
  including
  \item{prototypes}{a matrix containing the coordinates of the
    prototypes}
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
\seealso{See \code{\link{batchsom}} for Self-Oganising Map which
  provides both clustering and visualisation and \code{\link{kmeans}}
  for the version provided by the \code{stats} package.}
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

km <- batchkmeans(U,3)

plot(U,col=km$classif,cex=0.5)
points(km$prototypes,col=1:3,pch=20,cex=2)
}
\keyword{cluster}
