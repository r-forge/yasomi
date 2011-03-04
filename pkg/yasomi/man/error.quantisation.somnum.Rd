\name{error.quantisation.somnum}
\alias{error.quantisation.somnum}
\title{Quantisation error for a Self-Organising Map fitted on vector data}
\description{
Compute the quantisation error made by a fitted Self-Organising Map:
this is the mean of the Euclidean distance between each
observation and the prototype of its best matching unit.
}
\usage{
\method{error.quantisation}{somnum}(som, newdata, \dots)
}
\arguments{
  \item{som}{an object of class \code{"somnum"}}
  \item{newdata}{an optional matrix or a data frame of observations
    compatible with the fitted som}
  \item{\dots}{not used}
}
%\details{
%}
\value{
  If \code{newdata} is not given, the function returns the quantisation
  error made by the fitted som on the data used to fit it. When
  \code{newdata} is specified, the function returns the quantisation
  error of the
  fitted som on the corresponding data. 
}
%\references{}
\author{Fabrice Rossi}
%\note{
%}
\seealso{\code{\link{error.kaskilagus.somnum}}, \code{\link{som.tunecontrol}}}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a small hexagonal grid
sg <- somgrid(xdim=7,ydim=7,topo="hex")

# fit the SOM
som <- batchsom(data,sg)

print(error.quantisation(som))

# a larger grid should in general give a lower quantisation error
sg <- somgrid(xdim=12,ydim=12,topo="hex")
som <- batchsom(data,sg)
print(error.quantisation(som))
}
\keyword{cluster}
