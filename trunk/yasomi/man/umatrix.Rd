\name{umatrix}
\alias{umatrix}
\title{U-Matrix}
\description{
Plot a U-matrix for a Self-Organising Map.
}
\usage{
umatrix(som, ...)
}
\arguments{
  \item{som}{an object of class \code{"som"}}
  \item{\dots}{additional parameters transmitted to the low level plot function \code{\link{plot.somgrid}}}
}
%\details{
%}
%\value{
%}
\references{
  Ultsch, A. and Siemon, H. P. (1990) Kohonen's self organizing feature
  maps for exploratory data analysis, in: Proceedings of International
  Neural Network Conference (INNC'90).
}
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
\keyword{cluster}
\keyword{hplot}

