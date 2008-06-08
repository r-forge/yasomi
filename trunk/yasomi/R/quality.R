## quantisation error

error.quantisation.som <- function(som,newdata,...) {
    ## basic case without newdata
    if(!missing(newdata)) {
        stop("cannot compute the quantisation for new data")
    }
    som$errors[length(som$errors)]
}

error.quantisation.somnum <- function(som,newdata,...) {
    if(missing(newdata)) {
        NextMethod()
    } else {
        newdata <- as.matrix(newdata)
        ## data compatibility is verified in bmu
        bmu(som$prototypes,newdata)$error
    }
}

error.quantisation.relationalsom <- function(som,newdata,...) {
    if(missing(newdata)) {
        NextMethod()
    } else {
        ## data compatibility is verified in predict
        predict(som,newdata)$error
    }
}

error.quantisation.kernelsom <- function(som,newdata,norms,...) {
    if(missing(newdata)) {
        NextMethod()
    } else {
        baseerror <- predict(som,newdata)$error
        if(missing(norms)) {
            warning("cannot compute a full quantisation error with newdata norms")
        } else {
            if(length(norms)!=ncol(newdata)) {
                stop("'newdata' and 'norms' have different dimensions")
            }
            baseerror <- baseerror + mean(norms)
        }
        baseerror
    }
}

## Kasksi and Lagus' error measure

error.kaskilagus.somnum <- function(som,newdata,...) {
    if(missing(newdata)) {
        if(is.null(som$data)) {
            stop("cannot compute the quality measure without saved data or new data")
        }
        pre <- KaskiLagus.somnum(som,som$data,som$weights)
    } else {
        pre <- KaskiLagus.somnum(som,newdata)
    }
    mean(pre$quant+pre$path)
}

KaskiLagus.somnum <- function(som,data,weights) {
    data <- as.matrix(data)
    winners <- secondBmu(som$prototypes,data)
    pdist <- as.matrix(prototype.distances(som))
    if(missing(weights)) {
        list(quant=sqrt(rowSums((data-som$prototypes[winners[,1],])^2)),path=pdist[winners])
    } else {
        list(quant=sqrt(weights*rowSums((data-som$prototypes[winners[,1],])^2))/(mean(weights)),path=pdist[winners])
    }
}

error.kaskilagus.relationalsom <- function(som,newdata,...) {
    if(missing(newdata)) {
        pre <- KaskiLagus.relationalsom(som)
    } else {
        pre <- KaskiLagus.relationalsom(som,newdata)
    }
    mean(pre$quant+pre$path)
}

KaskiLagus.relationalsom <- function(som,newdata) {
    if(missing(newdata)) {
        ## on the original data, we use Dalpha and nf
        second <- relationalsecondbmu.R(som$Dalpha,som$nf,som$weights)
    } else {
        ## on new data we use predict in extended mode
        second <- predict(som,newdata,with.secondwinner=TRUE)
    }
    pdist <- as.matrix(prototype.distances(som))
    list(quant=sqrt(second$error),path=pdist[second$winners])
}

error.kaskilagus.kernelsom <- function(som,newdata,norms,...) {
    if(missing(newdata)) {
        pre <- KaskiLagus.kernelsom(som)
    } else {
        pre <- KaskiLagus.kernelsom(som,newdata,norms)
    }
    mean(pre$quant+pre$path)
}

KaskiLagus.kernelsom <- function(som,newdata,norms) {
    if(missing(newdata)) {
        ## on the original data, we use saved information
        second <- kernelsomsecondbmu(som$prototypes,som$data,som$weights)
    } else {
        if(length(norms)!=ncol(newdata)) {
            stop("'newdata' and 'norms' have different dimensions")
        }
        ## on new data we use predict in extended mode
        second <- predict(som,newdata,with.secondwinner=TRUE)
        second$error <- second$error + mean(norms)
    }
    pdist <- as.matrix(prototype.distances(som))
    list(quant=sqrt(second$error),path=pdist[second$winners])
}

