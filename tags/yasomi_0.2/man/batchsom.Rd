\name{batchsom}
\alias{batchsom}
\title{Generic Self-Organising Map fitting function}
\description{
\code{batchsom} is a generic function for fitting a Self-Organising Map
to some data, using a batch algorithm.
}
\usage{
batchsom(data, somgrid, init, prototypes, \dots)
}
\arguments{
  \item{data}{the data to which the SOM will be fitted. This can be,
    e.g., a matrix or data frame of observations (which should be scaled), or a
    distance matrix}
  \item{somgrid}{an object of class \code{"somgrid"} that specifies the prior
    structure of the Self-Organising Map: see \code{\link{somgrid}} }
  \item{init}{the initialisation method (see details)}
  \item{prototypes}{Initial values for the
    prototypes (the exact representation of the prototypes depends on
    the data type). If missing, initial prototypes are chosen via
    the method specified by the \code{init} parameter (see details)}
  \item{\dots}{additional arguments to be passed to methods}
}
\details{
If the initial value of \code{prototypes} is not provided, it is
obtained by a call to an initialisation method
specified via the \code{init} parameter (of character type). All
\code{batchsom} methods accept the \code{"random"} value for \code{init}
(see \code{\link{sominit.random}} for the corresponding
generic random initialisation method). all
\code{batchsom} methods accept also
the \code{"pca"} which covers ordered initialisation methods inspired by
principal component analysis in the standard vector case (see
\code{\link{sominit.pca}} for the corresponding generic method). The
default behaviour when neither \code{init} nor \code{prototypes} is
specified depends on the method. 
}
\value{
  An object of class \code{"som"}, a list with components including
  \item{somgrid}{as in the arguments to \code{batchsom}}
  \item{prototypes}{a representation of the prototypes that depends on
    the actual method}
  \item{classif}{a vector of integer indicating to which unit each
    observation has been assigned}
  \item{errors}{a vector containing the evolution of the quantisation
    error during the fitting process}
  \item{data}{the original data if the function is called with
    \code{keepdata = TRUE}}
  The list will generally contain additional components specific to each
  implementation. The returned object will also generally have another
  class more specific than \code{"som"}. 
}
%\references{ ~put references to the literature/web site here ~ }
\author{Fabrice Rossi}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{See \code{\link{sominit.pca}} and \code{\link{sominit.random}}
  for some control on the initial configuration of the prototypes, \code{\link{som.tune}} for the optimisation of some magic
  parameters (such as the radii), \code{\link{umatrix}} and
  \code{\link{distance.grid}} for visual analysis of the distances
  between the prototypes.}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a small hexagonal grid
sg <- somgrid(xdim=7,ydim=7,topo="hex")

# fit the SOM
som <- batchsom(data,sg)

# and display the prototypes
plot(som)
}
\keyword{cluster}
