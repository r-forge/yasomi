\name{batchsom}
\alias{batchsom}
\title{Generic Self-Organising Map fitting function}
\description{
Generic function for fitting a Self-Organising Map to some data, using a batch algorithm. 
}
\usage{
batchsom(data, somgrid, init=c("pca","random"), prototypes, weights,
         mode = c("continuous","stepwise"), min.radius, max.radius,
         steps, decrease = c("power", "linear"), max.iter,
         kernel = c("gaussian", "linear"), normalised,
         assignment = c("single", "heskes"), cut = 1e-07,
         verbose = FALSE, keepdata = TRUE, \dots)
}
\arguments{
  \item{data}{the data to which the SOM will be fitted. Acceptable data
    type depend on the available methods, see details}
  \item{somgrid}{an object of class \code{"somgrid"} that specifies the prior
    structure of the Self-Organising Map: see \code{\link{somgrid}} }
  \item{init}{the initialisation method (defaults to \code{"pca"}, see details)}
  \item{prototypes}{Initial values for the
    prototypes (the exact representation of the prototypes depends on
    the data type). If missing, initial prototypes are chosen via
    the method specified by the \code{init} parameter (see details)}
  \item{weights}{optional weights for the data points}
  \item{mode}{annealing mode:
    \describe{
      \item{\code{"continuous"} (default)}{this is the standard
	annealing strategy for SOM: the influence of neighbours changes
	at each epoch of the algorithm, from \code{max.radius} to
	\code{min.radius} in exactly \code{step} steps.}
      \item{\code{"stepwise"}}{in this strategy, the algorithm performs
	several epochs (a maximum of \code{max.iter}) for each of the
	\code{step} radii (from \code{max.radius} to
	\code{min.radius}). The algorithm changes the neighbours
	influence only when the classification remains stable from one
	epoch to another. The \code{max.iter} parameter provides a
	safeguard against cycling behaviours.}
    }
  }
  \item{min.radius}{the minimum neighbourhood influence radius. If
    missing, the value depends on the one of \code{kernel} but ensures
    in practice a local learning only (see details)}
  \item{max.radius}{the maximal neighbourhood influence radius. If missing
    two third of the prior structure diameter plus one}
  \item{steps}{the number of radii to use during annealing}
  \item{decrease}{the radii generating formula (\code{"power"} or
    \code{"linear"}), i.e., the way the \code{steps} radii are generated
    from the extremal values given by \code{min.radius} and
    \code{max.radius}}
  \item{max.iter}{maximal number of epochs for one radius in the
    \code{"stepwise"} annealing mode (defaults to 75)}
  \item{kernel}{the kernel used to transform distances in the prior
    structure into influence coefficients}
  \item{normalised}{switch for normalising the neighbouring
    interactions. Has no influence with the \code{"single"} assignment method}
  \item{assignment}{the assignment method used to compute the best
    matching unit (BMU) of an observation during training:
    \describe{
      \item{\code{"single"} (default)}{this is the standard BMU
	calculation approach in which the best unit for an observation
	is the one of the closest prototype of this observation }
      \item{\code{"heskes"}}{Tom Heskes' variant for the BMU in which a
	weighted fit of all the prototypes to an observation is used to
	compute the best unit. The rationale is that the BMU's prototype
	and its neighbouring units' prototypes must be close to the
	observation.}
    }
  }
  \item{cut}{minimal value below wich neighbouring interactions are not
    take into account}
  \item{verbose}{switch for tracing the fitting process}
  \item{keepdata}{if \code{TRUE}, the original data are returned as
    part of the result object}
  \item{\dots}{additional arguments to be passed to methods}
}
\details{
In yasomi, the \code{batchsom} generic function is implemented by three methods
which provide Self-Organising Map fitting for three distinct data
representation:
\itemize{
  \item the default implementation \code{\link{batchsom.default}} is
  used when the dataset \code{data} is given by a matrix or a data
  frame: it provides a standard batch SOM 
  implementation;
  \item for a dataset given by all the pairwise dissimilarities
  between the observations (\code{data} is an object of class
  \code{"dist"} from the \code{proxy} package), the method
  \code{\link{batchsom.dist}} implements 
  the \emph{relational} Self-Organising Map algorithm. It is
  based on prototypes represented by virtual linear combination of the
  original data points;
  \item when the dataset is given as a kernel matrix (\code{data} is an
  object of class \code{"kernelmatrix"}, see
  \code{\link{as.kernelmatrix}}), the method
  \code{\link{batchsom.kernelmatrix}} implements the
  batch \emph{kernel} Self-Organising Map algorithm. In this
  case, it is assumed that \code{data} contains all pairwise evaluation
  of a positive semi-definite kernel function and a batch SOM is
  performed (implicitly) in the kernel induced feature space.
  }
If the initial value of \code{prototypes} is not provided, it is 
obtained by a call to an initialisation method
specified via the \code{init} parameter (of character type). All
\code{batchsom} methods accept the \code{"random"} value for \code{init}
(see \code{\link{sominit.random}} for the corresponding
generic random initialisation method) as well as 
the \code{"pca"} value which covers ordered initialisation methods inspired by
principal component analysis in the standard vector case (see
\code{\link{sominit.pca}} for the corresponding generic method).
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
  \item{control}{a list containing all the parameters used to fit the SOM}
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
