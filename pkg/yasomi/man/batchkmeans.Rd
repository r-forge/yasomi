\name{batchkmeans}
\alias{batchkmeans}
\title{Generic K-means Clustering}
\description{
Generic function to perform K-means clustering on some data.
}
\usage{
batchkmeans(data, ncenters, init = c("prototypes", "random", "cluster"),
prototypes, weights, max.iter, verbose = FALSE, keepdata = TRUE, ...)
}
\arguments{
  \item{data}{the data to cluster. Acceptable data
    type depend on the available methods, see details}
  \item{ncenters}{the number of clusters}
  \item{init}{the initialisation method (see details)}
  \item{prototypes}{Initial values for the
    prototypes (the exact representation of the prototypes depends on
    the data type). If missing, initial prototypes are chosen via
    the method specified by the \code{init} parameter (see details)}
  \item{weights}{optional weights for the data points}
  \item{max.iter}{maximal number of iterations of the algorithm}
  \item{verbose}{switch for tracing the clustering process}
  \item{keepdata}{if \code{TRUE}, the original data are returned as
    part of the result object}
  \item{\dots}{additional arguments to be passed to methods}
}
\details{
  In yasomi, the \code{batchkmeans} generic function is implemented by
  two methods which provide K-means for two distinct data representation:
\itemize{
  \item the default implementation \code{\link{batchkmeans.default}} is
  used when the dataset \code{data} is given by a matrix or a data
  frame: it provides a standard (batch) K-means 
  implementation;
  \item when the dataset is given as a kernel matrix (\code{data} is an
  object of class \code{"kernelmatrix"}, see
  \code{\link{as.kernelmatrix}}), the method
  \code{\link{batchkmeans.kernelmatrix}} implements the
  (batch) \emph{kernel} K-means algorithm. In this
  case, it is assumed that \code{data} contains all pairwise evaluation
  of a positive semi-definite kernel function and a batch K-means clustering is
  performed (implicitly) in the kernel induced feature space.
}
If the initial value of \code{prototypes} is not provided, it is 
obtained by one of the following method specified by the \code{init}
parameter:
  \describe{
    \item{\code{"prototypes"}}{the standard method proceeds by choosing
    randomly a subset of the data of the requested size (with repetition
    if the grid size is larger than the data size). If the
    \code{weights} parameter is given, the probability of choosing a data
  point is proportionnal to its weight.}
  \item{\code{"random"}}{the \code{"random"} method generate prototypes
    randomly and uniformly in the hypercube spanned by the data for
    standard Euclidean data. For dissimilarity data or for the Kernel
    data, the method generates prototypes via random convex combinations
    of the data points. In all cases, the optional
  \code{weights} are not taken into account by this method.}
  \item{\code{"cluster"}}{the clustering initialisation method build a
    random partition the data into balanced clusters and uses as initial
  prototypes the centre of mass of those clusters. The optional
  \code{weights} are not taken into account for balancing the clusters.}
    }
  }
\value{
  An object of class \code{"batchkmeans"}, a list with components
  including
  \item{prototypes}{a representation of the prototypes that depends on
    the actual method}
  \item{classif}{a vector of integer indicating to which cluster each
    observation has been assigned}
  \item{errors}{a vector containing the evolution of the quantisation
    error during the fitting process}
  \item{data}{the original data if the function is called with
    \code{keepdata = TRUE}}
  \item{weights}{the weights of the data points if the function is called with
    \code{keepdata = TRUE} and if the \code{weights} is given}
  The list will generally contain additional components specific to each
  implementation. The returned object will also generally have another
  class more specific than \code{"batchkmeans"}.   
}
%\references{ ~put references to the literature/web site here ~ }
\author{Fabrice Rossi}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{See \code{\link{batchsom}} for Self-Oganising Map which
  provides both clustering and visualisation.}
%\examples{
%}
\keyword{cluster}
