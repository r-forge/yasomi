\name{umatrix}
\alias{umatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{U-Matrix}
\description{
Plot a U-matrix for a Self-Organising Map.
}
\usage{
umatrix(som, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{som}{an object of class \code{"som"}}
  \item{\dots}{additional parameters transmitted to the low level plot function \code{\link{plot.somgrid}}}
}
%\details{
%}
%\value{
%}
\references{Ultsch's papers}
\author{Fabrice Rossi}
%\note{ ~~further notes~~ 
%}
%\seealso{}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# a medium hexagonal grid
sg <- somgrid(xdim=15,ydim=15,topo="hex")

# choose a good SOM via Kaski and Lagus' error measure
som <-
som.tune(sg,data,som.tunecontrol(sg,,criterion=error.kaskilagus))$best.som

# display the umatrix
umatrix(som)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{hplot}

