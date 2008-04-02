\name{batchsom}
\alias{batchsom}
\alias{batchsom.default}
\title{Fit a Self-Organising Map to some data}
\description{
\code{batchsom} is used to fit a Self-Organising Map to some data, using a
batch algorithm. This is a generic method with a default implementation
for standard data.
}
\usage{
batchsom(data, somgrid, init, prototypes, \dots)
\method{batchsom}{default}(data, somgrid, init=c("pca","random"), prototypes,
         assignment = c("single", "heskes"),radii,nbRadii=30,
         maxiter = 75,kernel=c("gaussian","linear"),normalised,
         cut = 1e-07, verbose = FALSE, keepdata = TRUE, \dots)
}
\arguments{
  \item{data}{the data to which the SOM will be fitted. This can be,
    e.g., a matrix or data frame of observations (which should be scaled), or a
    distance matrix}
  \item{somgrid}{an object of class \code{"somgrid"} that specifies the prior
    structure of the Self-Organising Map: see \code{\link{somgrid}} }
  \item{prototypes}{a matrix of initial values for the
    prototypes. If missing, initial prototypes are chosen via
    the method specified by the \code{init} parameter (see details)}
  \item{init}{the initialisation method (see details)}
  \item{assignment}{the assignment method used to compute the best
    matching unit of an observation during training}
  \item{radii}{the radii of the neighbourhoods used during training. If
    missing, \code{nbRadii} radii are computed starting from a radius of 2/3
    of the diameter of the prior structure and ending with a radius small
    enough to prevent any influence of a node on its neighbours}
  \item{nbRadii}{the number of radii to generate if 'radii' is missing}
  \item{maxiter}{the maximum number of of iteration performed for each
    radius (the minimum is one)}
  \item{kernel}{the kernel to use to compute neighbouring iterations}
  \item{normalised}{switch for normalising the neighbouring
    interactions. Has no influence with the \code{"single"} assignment method}
  \item{cut}{minimal value below wich neighbouring interactions are not
    take into account}
  \item{verbose}{switch for tracing the fitting process}
  \item{keepdata}{if \code{TRUE}, the original data are returned as
    part of the result object}
  \item{\dots}{additional arguments to be passed to the initialisation method}
}
\details{
The default method is the standard batch algorithm for Self-Organising
Map. It applies to a matrix or a data frame of observations. In this
case, it is recommended to use the \code{"single"} assignment method
which seems to perform better in general. If the initial value of
\code{prototypes} is not provided, it is obtained by a call to
\code{\link{sominit.pca.default}}. 
}
\value{
  An object of class \code{"som"}, a list with components including
  \item{somgrid}{as in the arguments to \code{batchsom}}
  \item{prototypes}{a matrix containing the coordinates of the
    prototypes (also called the code book vectors)}
  \item{classif}{a vector of integer indicating to which unit each
    observation has been assigned}
  \item{errors}{a vector containing the evolution of the quantisation
    error during the fitting process}
  \item{data}{the original data if the function is called with
    \code{keepdata = TRUE}}
  In addition, for the default method, the object is also of class
  \code{"somnum"}. 
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
