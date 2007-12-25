\name{som.tune}
\alias{som.tune}
\title{Parameter tuning for Self-Organising Maps}
\description{
This function tunes some parameters of a Self-Organising Map by
optimising a specified error measure.  The prior structure is not
optimised by this function.
}
\usage{
som.tune(data, somgrid, control = som.tunecontrol(somgrid), verbose = FALSE)
}
\arguments{
  \item{data}{a matrix or a data frame of observations}
  \item{somgrid}{an object of class 'somgrid'}
  \item{control}{a control object which specifies the parameter to optimize}
  \item{verbose}{switch for tracing the tuning process}
}
\details{
The optimisation is done via a simple grid search that fit a
Self-Organising Map for each possible combination of the parameter to
tune. The prior structure itself is excluded from the optimisation as
most quantisation related error measures will decrease with the size of
the grid.
}
\value{
  A list with components
  \item{errors}{the error for each configuration}
  \item{best.som}{the best SOM according to the chosen error criterion}
}
\author{Fabrice Rossi}
\note{this function is modelled after David Mayer's \link[e1071]{tune} function.
}
\seealso{\code{\link{som.tunecontrol}}}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a medium hexagonal grid
sg <- somgrid(xdim=15,ydim=15,topo="hex")

# choose a good SOM via Kaski and Lagus' error measure
tune.results <-
  som.tune(data,sg,som.tunecontrol(sg,criterion=error.kaskilagus))

# plot errors
plot(tune.results$errors,type="h")
}
\keyword{cluster}

