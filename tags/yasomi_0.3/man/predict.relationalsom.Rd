\name{predict.relationalsom}
\alias{predict.relationalsom}
\title{Map new data to a fitted Relational Self-Organising Map}
\description{
Map new data to a fitted Self-Organising Map, i.e., compute the best
matching unit for each given observation. 
}
\usage{
\method{predict}{relationalsom}(object, newdata,
        with.secondwinner=FALSE,\dots)
}
\arguments{
  \item{object}{an object of class \code{"relationalsom"}}
  \item{newdata}{an object of class \code{"crossdist"}
    compatible with the fitted som (see details)}
  \item{with.secondwinner}{switch to specify whether the second best
    matching unit should be computed and returned}
  \item{\dots}{not used}
}
\details{
  The \code{newdata} object must be of class
  \code{"crossdist"} as returned by \code{\link[proxy]{dist}} and must
  contain the dissimilarities between the original data (used to fit the
  SOM) and the new data (for which the quantisation error is to be
  computed). 
}
\value{
  A list with components
  \item{classif}{a vector of integer indicating to which unit each
    observation has been assigned}
  \item{error}{the quantisation error of the observations by the
    prototypes of this Self-Organising Map (see
    \code{\link{error.quantisation.relationalsom}} for a definition of
    this value)}
  \item{rdist}{a matrix of squared dissimilarities between the new data
    and the prototypes, as obtained by the relational formula (with
    possible negative values preserved)}
  \item{winners}{if \code{with.secondwinner=TRUE} a two column matrix
    with a line per observation containing the index of the best
    matching unit and of the second best matching unit (in this order)}
}
\author{Fabrice Rossi}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# use only part of the data
train <- sample(1:nrow(data),100)
data.train <- data[train,]
d <- dist(data.train)

# a small hexagonal grid
sg <- somgrid(xdim=5,ydim=5,topo="hex")

# fit the SOM
som <- batchsom(d,sg)

# map remaining data
results <- predict(som,dist(data.train,data[-train,]))
print(paste("Learning quantisation error:",error.quantisation(som)))
print(paste("Test quantisation error:",results$error))
}
\keyword{cluster}
