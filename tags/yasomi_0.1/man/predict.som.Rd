\name{predict.som}
\alias{predict.som}
\title{Map new data to a fitted Self-Organising Map}
\description{
Map new data to a fitted Self-Organising Map, i.e., compute the best
matching unit for each given observation. 
}
\usage{
\method{predict}{som}(object, newdata, ...)
}
\arguments{
  \item{object}{an object of class \code{"som"}}
  \item{newdata}{a matrix or a data frame of data compatible with the data
    used to fit the \code{object} Self-Organising Map (that is with
    identical dimension)}
  \item{\dots}{not used}
}
%\details{
%}
\value{
  A list with components
  \item{classif}{a vector of integer indicating to which unit each
    observation has been assigned}
  \item{error}{the quantisation error of the observations by the
    prototypes of this Self-Organising Map}
}
\author{Fabrice Rossi}
\examples{
data(iris)
# scaling
data <- scale(iris[1:4])

# use only part of the data
train <- sample(1:nrow(data),100)
data.train <- data[train,]

# a small hexagonal grid
sg <- somgrid(xdim=5,ydim=5,topo="hex")

# fit the SOM
som <- batchsom(sg,data)

# map remaining data
results <- predict(som,data[-train,])
print(paste("Quantisation error:",results$error))
}
\keyword{cluster}

