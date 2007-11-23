\name{prototype.distances}
\alias{prototype.distances}
\title{Compute distances between prototypes in a fitted Self-Organising Map}
\description{
Compute the Euclidean distances between the prototypes of a fitted
Self-Organising Map for neighbouring units in the prior structure.
}
\usage{
prototype.distances(som)
}
\arguments{
  \item{som}{an object of class \code{"som"}}
}
\details{
  The distance matrix has one row for each prototype and 6 or 8 columns,
  depending on the topology of the prior structure (6 for a hexagonal
  grid and 8 for a rectangular one). When a prototype has not the
  maximum number of neighbours because of its position in the grid,
  the corresponding entries are encoded with \code{NA}. 
}
\value{
  An object of class\code{"sompdist"}, a list with components
  \item{pdist}{a matrix containing the distances between prototypes. }
  \item{somgrid}{an object of class \code{"somgrid"}, the prior
    structure used for fitting the \code{som} object}
}
%\references{}
\author{Fabrice Rossi}
%\note{
%}
\seealso{\code{\link{error.kaskilagus}}}
\examples{
}
\keyword{cluster}

