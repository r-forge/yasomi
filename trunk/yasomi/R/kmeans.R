prototypes.random <- function(data,nb,method=c("prototypes","random","cluster"),
                              weights) {
    method <- match.arg(method)
    data <- as.matrix(data)
    nb.data <- nrow(data)
    if(method=="prototypes" || (method=="cluster" && nb>=nb.data)) {
        if(missing(weights) || is.null(weights)) {
            data[sample(nb.data,size=nb,replace=nb>nb.data),]
        } else {
            data[sample(nb.data,size=nb,replace=nb>nb.data,prob=weights),]
        }
    } else if(method=="random") {
        ## weights do not need to be taken into account here
        result <- matrix(0,ncol=ncol(data),nrow=nb)
        for(i in 1:ncol(data)) {
            therange <- range(data[,i])
            result[,i] <- runif(nb,min=therange[1],max=therange[2])
        }
        result
    } else {
        ## nb <nb.data
        result <- matrix(0,ncol=ncol(data),nrow=nb)
        if(missing(weights) || is.null(weights)) {
            clusters <- cut(sample(nb.data),nb,labels=FALSE,include.lowest=TRUE)
            for(i in 1:nb) {
                result[i,] <- colMeans(as.matrix(data[clusters==i,],ncol=ncol(data)))
            }
        } else {
            totalWeight <- sum(weights)
            subsets <- sample(nb.data)
            breaks <- seq(from=0,to=totalWeight,length.out=nb+1)
            sums <- cumsum(weights[subsets])
            for(i in 1:nb) {
                the.subset <- subsets[sums>breaks[i] & sums<=breaks[i+1]]
                result[i,] <- apply(sweep(matrix(data[the.subset,],ncol=ncol(data)),1,weights[the.subset],"*"),2,sum)/sum(weights[the.subset])
            }
        }
        result
    }
}

batchkmeans.default <- function(data,ncenters,
                                init=c("prototypes","random","cluster"),
                                prototypes,weights,max.iter=75,
                                verbose=FALSE,keepdata=TRUE,...) {
    if(missing(prototypes)) {
        prototypes <- prototypes.random(data,ncenters,init,weights)
    } else {
        if(ncol(prototypes)!=ncol(data)) {
            stop("'prototypes' and 'data' have different dimensions")
        }
        if(nrow(prototypes)!=ncenters) {
            stop("'prototypes' and 'ncenters' are not compatible")
        }
    }
    if(!missing(weights)) {
        if(length(weights)!=nrow(data)) {
            stop("'weights' and 'data' have different dimensions")
        }
    } else {
        weights <- rep(1,nrow(data))
    }
    classif <- rep(NA,nrow(data))
    errors <- rep(NA,max.iter)
    for(i in 1:max.iter) {
        bmus <- bmu(prototypes,data,weights)
        nclassif <- bmus$clusters
        noChange <- identical(classif,nclassif)
        classif <- nclassif
        errors[i] <- bmus$error
        if(verbose) {
            print(paste(i,errors[i]))
        }
        if(noChange) {
            break;
        }
        ## prototypes update (this is slow)
        for(i in 1:ncenters) {
            theclass <- classif==i
            prototypes[i,] <- apply(sweep(data[theclass,],1,weights[theclass],"*"),2,sum)/sum(weights[theclass])
        }
    }
    if(!noChange) {
        bmus <- bmu(prototypes,data,weights)
        classif <- bmus$clusters
        errors <- c(errors,bmus$errors)
        print(paste("warning: can't reach a stable configuration after",max.iter,"iterations"))
    } else {
        errors <- errors[!is.na(errors)]
    }
    res <- list(prototypes=prototypes,classif=classif,errors=errors)
    if(keepdata) {
        res$data  <- data
        res$weights <- weights
    }
    class(res) <- c("batchkmeansnum","batchkmeans")
    res

}

error.quantisation.batchkmeans <- function(som,newdata,...) {
    ## basic case without newdata
    if(!missing(newdata)) {
        stop("cannot compute the quantisation for new data")
    }
    som$errors[length(som$errors)]
}

convex.prototypes.random <- function(data,nb,
                                     method=c("prototypes","random","cluster"),
                                     weights) {
    method <- match.arg(method)
    nb.data <- nrow(data)
    if(method=="prototypes" || (method=="cluster" && nb>=nb.data)) {
        protos <- matrix(0,ncol=nb.data,nrow=nb)
        if(missing(weights) || is.null(weights)) {
            protos[cbind(1:nb,sample(nb.data,size=nb,replace=nb>nb.data))] <- 1
        } else {
            protos[cbind(1:nb,sample(nb.data,size=nb,replace=nb>nb.data,prob=weights))] <- 1
        }
    } else if(method=="random") {
        if(missing(weights) || is.null(weights)) {
            protos <- matrix(runif(nb.data*nb),ncol=nb.data,nrow=nb)
        } else {
            protos <- matrix(runif(nb.data*nb,max=rep(weights,each=nb)),ncol=nb.data,nrow=nb)
        }
        protos <- sweep(protos,1,rowSums(protos),"/")
    } else {
        ## nb <nb.data
        protos <- matrix(0,ncol=nb.data,nrow=nb)
        if(missing(weights) || is.null(weights)) {
            clusters <- cut(sample(nb.data),nb,labels=FALSE,include.lowest=TRUE)
            protos[cbind(clusters,1:nb.data)] <- 1
            protos <- sweep(protos,1,rowSums(protos),"/")
        } else {
            totalWeight <- sum(weights)
            subsets <- sample(nb.data)
            breaks <- seq(from=0,to=totalWeight,length.out=nb+1)
            sums <- cumsum(weights[subsets])
            for(i in 1:nb) {
                the.subset <- subsets[sums>breaks[i] & sums<=breaks[i+1]]
                protos[i,the.subset] <- weights[the.subset]
            }
            protos <- sweep(protos,1,rowSums(protos),"/")
        }
    }
    protos
}

fastKernelKmeansbmu <- function(cluster,nclust,K) {
###FIXME: could be faster (e.g., only diag(ps$bips) is needed)
    ps <- partialSums(cluster,nclust,K)
    csize <- table(factor(cluster,levels=1:nclust))
    Kp <- sweep(ps$ps,1,csize,"/")
    pnorms <- diag(ps$bips)/(csize^2)
    predistances <- sweep(-2*Kp,1,pnorms,"+")
    bmu <- apply(predistances,2,which.min)
    error <- mean(predistances[cbind(bmu,1:length(bmu))]+diag(K))
    list(clusters=bmu,error=error,Kp=t(Kp),pnorms=pnorms)
}

weightedFastKernelKmeansbmu <- function(cluster,nclust,K,weights) {
###FIXME: could be faster (e.g., only diag(ps$bips) is needed)
    ps <- weightedPartialSums(cluster,nclust,K,weights)
    csize <- tapply(weights,factor(cluster,levels=1:nclust),sum)
    Kp <- sweep(ps$ps,1,csize,"/")
    pnorms <- diag(ps$bips)/(csize^2)
    predistances <- sweep(-2*Kp,1,pnorms,"+")
    bmu <- apply(predistances,2,which.min)
    error <- sum(weights*(predistances[cbind(bmu,1:length(bmu))]+diag(K)))/sum(weights)
    list(clusters=bmu,error=error,Kp=t(Kp),pnorms=pnorms)
}

batchkmeans.kernelmatrix <- function(data,ncenters,
                                     init=c("prototypes","random","cluster"),
                                     prototypes,weights,max.iter=75,
                                     verbose=FALSE,keepdata=TRUE,...) {
    if(missing(prototypes)) {
        prototypes <- convex.prototypes.random(data,ncenters,init,weights)
    } else {
        if(ncol(prototypes)!=ncol(data)) {
            stop("'prototypes' and 'data' have different dimensions")
        }
        if(nrow(prototypes)!=ncenters) {
            stop("'prototypes' and 'ncenters' are not compatible")
        }
    }
    if(!missing(weights)) {
        if(length(weights)!=nrow(data)) {
            stop("'weights' and 'data' have different dimensions")
        }
    } else {
        weights <- NULL
    }
    classif <- rep(NA,nrow(data))
    errors <- rep(NA,max.iter)
    ## initialisation round
    bmus <- kernelsombmu(prototypes,data,weights)
    classif <- bmus$clusters
    errors[1] <- bmus$error
    if(verbose) {
        print(paste(1,bmus$error))
    }
    for(i in 2:max.iter) {
        if(is.null(weights)) {
            bmus <- fastKernelKmeansbmu(classif,ncenters,data)
        } else {
            bmus <- weightedFastKernelKmeansbmu(classif,ncenters,data,weights)
        }
        nclassif <- bmus$clusters
        noChange <- identical(classif,nclassif)
        classif <- nclassif
        errors[i] <- bmus$error
        if(verbose) {
            print(paste(i,errors[i]))
        }
        if(noChange) {
            break;
        }
    }
    ## for latter use
    if(is.null(weights)) {
        csize <- table(factor(classif,levels=1:ncenters))
        for(i in 1:ncenters) {
            prototypes[i,] <- (classif==i)/csize[i]
        }
    } else {
        csize <- tapply(weights,factor(classif,levels=1:ncenters),sum)
        for(i in 1:ncenters) {
            mask <- classif==i
            prototypes[i,mask] <- weights[mask]/csize[i]
            prototypes[i,!mask] <- 0
        }
    }
###FIXME: this is slow
    Kp <- tcrossprod(data,prototypes)
    pnorms <- double(nrow(prototypes))
    for(i in 1:length(pnorms)) {
        pnorms[i] <- c(prototypes[i,]%*%Kp[,i])
    }
    if(!noChange) {
        if(is.null(weights)) {
            bmus <- fastKernelKmeansbmu(classif,ncenters,data)
        } else {
            bmus <- weightedFastKernelKmeansbmu(classif,ncenters,data,weights)
        }
        classif <- bmus$clusters
        errors <- c(errors,bmus$errors)
        print(paste("warning: can't reach a stable configuration after",max.iter,"iterations"))
    } else {
        errors <- errors[!is.na(errors)]
    }
    res <- list(prototypes=prototypes,classif=classif,errors=errors,Kp=Kp,pnorms=pnorms)
    if(keepdata) {
        res$data  <- data
        res$weights <- weights
    }
    class(res) <- c("kernelbatchkmeans","batchkmeans")
    res
}

