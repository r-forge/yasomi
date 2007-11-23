\name{error.quantisation}
\alias{error.quantisation}
\title{Quantisation error for a Self-Organising Map}
\description{
Compute the quantisation error made by a fitted Self-Organising Map:
this is the mean of the squared Euclidean distance between each
observation and the prototype of its best matching unit.
}
\usage{
error.quantisation(som, data)
}
\arguments{
  \item{som}{an object of class \code{"som"}}
  \item{data}{a matrix or a data frame of observations}
}
%\details{
%}
\value{
  The quantisation error.
}
%\references{}
\author{Fabrice Rossi}
%\note{
%}
\seealso{\code{\link{error.kaskilagus}}, \code{\link{som.tunecontrol}}}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a small hexagonal grid
sg <- somgrid(xdim=7,ydim=7,topo="hex")

# fit the SOM
som <- batchsom(sg,data)

print(error.quantisation(som,data))

# a larger grid should in general give a lower quantisation error
sg <- somgrid(xdim=12,ydim=12,topo="hex")
som <- batchsom(sg,data)
print(error.quantisation(som,data))
}
\keyword{cluster}
