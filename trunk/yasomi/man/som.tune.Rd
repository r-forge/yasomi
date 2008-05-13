\name{som.tune}
\alias{som.tune}
\alias{print.somtune}
\title{Parameter tuning for Self-Organising Maps}
\description{
This function tunes some parameters of a Self-Organising Map by
optimising a specified error measure. The prior structure is not
optimised by this function.
}
\usage{
som.tune(data, somgrid, control = som.tunecontrol(somgrid),
         verbose = FALSE, internalVerbose = FALSE)
\method{print}{somtune}(x,\dots)
}
\arguments{
  \item{data}{the data to which the SOM will be fitted. This can be,
    e.g., a matrix or data frame of observations (which should be scaled), or a
    distance matrix}
  \item{somgrid}{an object of class 'somgrid'}
  \item{control}{a control object which specifies the parameter to optimize}
  \item{verbose}{switch for tracing the tuning process}
  \item{internalVerbose}{switch for tracing each call to
    \code{\link{batchsom}}}
  \item{x}{an object of class \code{"somtune"}}
  \item{\dots}{not used}
}
\details{
The optimisation is done via a simple grid search that fit a
Self-Organising Map for each possible combination of the parameters to
tune. The prior structure itself is excluded from the optimisation as
most quantisation related error measures decrease with the size of
the grid.
}
\value{
  An object of class \code{"somtune"}, a list with components
  \item{best.som}{the best SOM according to the chosen error criterion}
  \item{errors}{the error for each configuration}
  \item{quantisation}{the quantisation error for each configuration}
  \item{isquant}{\code{TRUE} if the best SOM was chosen according to the
    quantisation error}
  \item{control}{the control object used by this call}
  \item{dimensions}{a list of strings with the names of the parameters
    that were varied by the function}
  \item{init,assignement,radii,annealing,kernel}{5 vectors containing
    the parameters used by each tested configuration}
  \item{best.index}{the index of the configuration used by the best SOM}
}
\author{Fabrice Rossi}
\note{this function is modelled after David Mayer's \link[e1071]{tune} function.
}
\seealso{\code{\link{som.tunecontrol}} and \code{\link{plot.somtune}}}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a medium hexagonal grid
sg <- somgrid(xdim=15,ydim=15,topo="hex")

# choose a good SOM via Kaski and Lagus' error measure
tune.results <-
  som.tune(data,sg,som.tunecontrol(sg,criterion=error.kaskilagus))

# print the results
print(tune.results)

# umatrix of the best SOM
umatrix(tune.results$best.som)
}
\keyword{cluster}

