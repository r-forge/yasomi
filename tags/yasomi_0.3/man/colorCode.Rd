\name{colorCode}
\alias{colorCode}
\title{Data coloring via a Self-Organising Map}
\description{
\code{colorcode} maps observations to integers in an organised way:
close observations should be assigned close integers (and vice versa).
The coding is done by fitting a one dimensional Self-Organising Map to
the data.  The main use of this method is to provide color based linking
between views. 
}
\usage{
colorCode(data, nbcolor)
}
\arguments{
  \item{data}{a matrix or a data frame of observations}
  \item{nbcolor}{the number of colors in the resulting mapping (i.e., the output will belong to [1,\code{nbcolor}]}
}
%\details{
%}
\references{
  Vesanto, J. (2002) \emph{Data Exploration Process Based on the
    Self--Organizing Map}, PhD thesis, Helsinki University of Technology, Espoo (Finland), Acta Polytechnica Scandinavica, Mathematics and Computing Series No.
  115.
}
\author{Fabrice Rossi}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# code calculation
code <- colorCode(data,50)

# display the result
pairs(data,bg=rainbow(50)[code],pch=21)

}
\keyword{cluster}
\keyword{dplot}
