\name{batchsom.default}
\alias{batchsom.default}
\title{Fit a Self-Organising Map to vector data}
\description{
\code{batchsom.default} is used to fit a Self-Organising Map to vector
data, using a batch algorithm. 
}
\usage{
\method{batchsom}{default}(data, somgrid, init=c("pca","random"), prototypes,
         assignment = c("single", "heskes"),radii,nbRadii=30,
         maxiter = 75,kernel=c("gaussian","linear"),normalised,
         cut = 1e-07, verbose = FALSE, keepdata = TRUE, \dots)
}
\arguments{
  \item{data}{the data to which the SOM will be fitted. This can be,
    a matrix or data frame of observations (which should be scaled)}
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
This methods implements the standard batch algorithm for Self-Organising
Map. It applies to a matrix or a data frame of observations. It is
recommended to use the \code{"single"} assignment method 
which seems to perform better in general. If the initial value of
\code{prototypes} is not provided, it is obtained by a call to
a function specified by the \code{init} parameter. If its value is
\code{"pca"} proprototypes are obtained by a call to 
\code{\link{sominit.pca.default}} (this is also the case when
\code{init} is not specified), while
\code{\link{sominit.random.default}} is called when \code{init} is
\code{"random"}. In both case, the additional parameters submitted to
the method are transmitted to the initialisation method. 
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
  The object is also of class
  \code{"somnum"}. 
}
%\references{ ~put references to the literature/web site here ~ }
\author{Fabrice Rossi}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a small hexagonal grid
sg <- somgrid(xdim=7,ydim=7,topo="hex")

# fit the SOM with random initial prototypes
som <- batchsom(data,sg,init="random")

# and display the prototypes
plot(som)
}
\keyword{cluster}
