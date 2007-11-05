\name{grid2lines}
\alias{grid2lines}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Self-Organising Map grid}
\description{
Build a matrix that can be used to plot the direct neighbourhood relationship in the prior structure of a Self-Organising Map.
}
\usage{
grid2lines(som)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{som}{an object of class \code{"som"}}
}
\details{
\code{grid2lines} returns a transformed version of the coordinate matrix
of the prototypes of a fitted Self-Organising Map.  It can be used with, e.g., \code{\link{lines}} to draw edges between prototypes that are direct neighbours in the prior structure. 
}
\value{
a matrix with the same number of columns as the original dataset used by the \code{som} object. 
}
%\references{ ~put references to the literature/web site here ~ }
\author{Fabrice Rossi}
%\note{ ~~further notes~~ 
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
X <- cbind(rnorm(500),rnorm(500))

sg <- somgrid(xdim=7,ydim=7,topo="rect")

som <- batchsom(sg,X)

# display the data, the som's prototypes and the prior structure
plot(X,pch="+",col="red")
lines(grid2lines(som),type="b",pch=20)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{dplot}
