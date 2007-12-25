error.kaskilagus <- function(som,data) {
    pre <- KaskiLagus(som,data)
    mean(pre$quant+pre$path)
}

error.quantisation <- function(som,data) {
    bmu(som$prototypes,data)$error
}

## numerical SOM
KaskiLagus.somnum <- function(som,data) {
    data <- as.matrix(data)
    winners <- secondBmu(som$prototypes,data)
    pdist <- as.matrix(prototype.distances(som))
    list(quant=sqrt(rowSums((data-som$prototypes[winners[,1],])^2)),path=pdist[winners])
}

## relational SOM
KaskiLagus.relationalsom <- function(som,data) {
    data <- as.matrix(data^2,diag=0)
    second <- relationalsecondbmu.R(som$prototypes,data)
    pdist <- as.matrix(prototype.distances(som))
    list(quant=second$error,path=pdist[second$winners])
}

