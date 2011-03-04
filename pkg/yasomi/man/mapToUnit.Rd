\name{mapToUnit}
\alias{mapToUnit}
\title{Map additional data to a Self-Organising Map}
\description{
\code{mapToUnit} is used to map additional data to a fitted
Self-Organising Map.  The main use case is the mapping of variables that were not used during the fitting (e.g., factors).
}
\usage{
mapToUnit(som, values)
}
\arguments{
  \item{som}{an object of class \code{"som"}}
  \item{values}{a matrix or a data frame of observations to map on to
    the SOM's unit.  The number of observations must be equal to the number of observations used to fit the SOM}
}
\details{
The goal of this function is not to assign new observations to a fitted
SOM but rather to map additional variables (see
\code{\link{predict.somnum}} for new observations assignment).  The mapping is based only on the clustering result of the SOM: if observation number k was assigned to unit l at the end of the fitting process, then \code{values[k,]} is mapped to the l-th component of the resulting list. 
}
\value{
a list of length the number of units in the SOM. Each component contains the subset of \code{values} that is mapped to the corresponding unit. 
}
%\references{ ~put references to the literature/web site here ~ }
\author{Fabrice Rossi}
%\note{
%}
\seealso{\code{\link{predict.somnum}} for new observations assignment}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a small hexagonal grid
sg <- somgrid(xdim=7,ydim=7,topo="hex")

# fit the SOM
som <- batchsom(data,sg)

# map the class information and print the content of each cell
print(mapToUnit(som,iris[[5]]))
}
\keyword{cluster}

