error.kaskilagus <- function(som,newdata) {
    if(missing(newdata)) {
        if(is.null(som$data)) {
            stop("cannot compute the quality measure without saved data or new data")
        }
        newdata <- som$data
    }
    pre <- KaskiLagus(som,newdata)
    mean(pre$quant+pre$path)
}

## works partially for relationSOM
error.quantisation <- function(som,newdata) {
    if(missing(newdata)) {
        som$errors[length(som$errors)]
    } else {
        newdata <- as.matrix(newdata)
        ## this is also verified in bmu
        if(ncol(newdata) != ncol(som$prototypes)) {
            stop("'newdata' and 'som$prototypes' have different dimensions")
        }
        bmu(som$prototypes,newdata)$error
    }
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

