as.kernelmatrix.matrix <- function(data,...) {
### FIXME: check for positivity
    if(!is.matrix(data)) {
        stop("'data' is not a matrix")
    }
    if(!isSymmetric(data)) {
        stop("'data' is not symmetric")
    }
    class(data) <- c("kernelmatrix",class(data))
    data
}

as.kernelmatrix.dist <- function(data,...) {
### FIXME: check for positivity
    result <- -0.5*double.centering(as.matrix(data^2,diag=0))
    class(result) <- c("kernelmatrix",class(result))
    result
}

as.dist.kernelmatrix <- function(x,FUN=NULL) {
    predist <- sweep(-2*x,1,diag(x),"+")
    predist <- sweep(predist,2,diag(x),"+")
    class(predist) <- c("matrix")
    as.dist(sqrt(predist))
}

predict.kernelsom <- function(object,newdata,newdatanorms,with.secondwinner=FALSE,...) {
    if(nrow(newdata)!=ncol(object$prototypes)) {
        stop("'newdata' and 'object$prototypes' have different dimensions")
    }
    if(length(newdatanorms)!=ncol(newdata)) {
        stop("'newdata' and 'newdatanorms' have incompatible sizes")
    }
    nKp <- object$prototypes%*%newdata
    distances <- t(sweep(sweep(-2*nKp,1,object$pnorms,"+"),2,newdatanorms,"+"))
    if(with.secondwinner) {
### FIXME: very suboptimal
        ordered <- apply(distances,1,order)
        winners <- t(ordered[1:2,])
        bmu <- winners[,1]
    } else {
        bmu <- apply(distances,1,which.min)
    }
    error <- mean(relational.keeppositive(distances[cbind(1:length(bmu),bmu)]))
    if(with.secondwinner) {
        list(classif=bmu,error=error,distances=distances,winners=winners)
    } else {
        list(classif=bmu,error=error,distances=distances)
    }
}

kernelsomsecondbmu <- function(prototypes,K) {
    Kp <- tcrossprod(K,prototypes)
    pnorms <- double(nrow(prototypes))
    for(i in 1:length(pnorms)) {
        pnorms[i] <- c(prototypes[i,]%*%Kp[,i])
    }
    predistances <- sweep(-2*Kp,2,pnorms,"+")
### FIXME: very suboptimal
    ordered <- apply(predistances,1,order)
    winners <- t(ordered[1:2,])
    bmu <- winners[,1]
    error <- predistances[cbind(1:length(bmu),bmu)]+diag(K)
    list(clusters=bmu,error=error,pnorms=pnorms,winners=winners)
}

fastKernelsombmu <- function(cluster,nclust,K,nv) {
    ps <- partialSums(cluster,nclust,K)
    csize <- table(factor(cluster,levels=1:nclust))
    normed <- nv%*%csize
    nvnormed <- sweep(nv,1,normed,"/")
    innerProducts <- nvnormed%*%ps$ps
    interm <- tcrossprod(ps$bips,nvnormed)
    pnorms <- double(nclust)
    for(i in 1:nclust) {
        pnorms[i] <- c(nvnormed[i,]%*%interm[,i])
    }
    predistances <- sweep(-2*innerProducts,1,pnorms,"+")
    bmu <- apply(predistances,2,which.min)
    error <- mean(predistances[cbind(bmu,1:length(bmu))]+diag(K))
    list(clusters=bmu,error=error,Kp=t(innerProducts),pnorms=pnorms)
}

weightedFastKernelsombmu <- function(cluster,nclust,K,nv,weights) {
    ps <- weightedPartialSums(cluster,nclust,K,weights)
    csize <- tapply(weights,factor(cluster,levels=1:nclust),sum)
    csize[is.na(csize)] <- 0
    normed <- nv%*%csize
    nvnormed <- sweep(nv,1,normed,"/")
    innerProducts <- nvnormed%*%ps$ps
    interm <- tcrossprod(ps$bips,nvnormed)
    pnorms <- double(nclust)
    for(i in 1:nclust) {
        pnorms[i] <- c(nvnormed[i,]%*%interm[,i])
    }
    predistances <- sweep(-2*innerProducts,1,pnorms,"+")
    bmu <- apply(predistances,2,which.min)
    error <- sum(weights*(predistances[cbind(bmu,1:length(bmu))]+diag(K)))/sum(weights)
    list(clusters=bmu,error=error,Kp=t(innerProducts),pnorms=pnorms)
}


fastKernelsom.lowlevel.R <- function(somgrid,K,prototypes,
                                     assignment,radii,weights,maxiter,
                                     kernel,normalised,cut,verbose) {
    oldClassif <- rep(NA,nrow(K))
    errors <- vector("list",length(radii))
    ## a round of initialisation is needed
    nv <- neighborhood(somgrid,radii[1],kernel,normalised=normalised)
    bmus <- kernelsombmu(prototypes,K,weights)
    classif <- bmus$clusters
    errors[[1]] <- bmus$error
    if(verbose) {
        print(paste(1,1,bmus$error))
    }
    for(i in 1:length(radii)) {
        if(i==1) {
            iterations <- 2:maxiter
        } else {
            iterations <- 1:maxiter
        }
        for(j in iterations) {
            ## assignment
            if(assignment == "single") {
                if(is.null(weights)) {
                    bmus <-  fastKernelsombmu(classif,somgrid$size,K,nv)
                } else {
                    bmus <-  weightedFastKernelsombmu(classif,somgrid$size,K,nv,weights)
                }
            } else {
                stop(paste(assignment,"is not implemented for kernel SOM"))
            }
            nclassif <- bmus$clusters
            noChange <- identical(classif,nclassif)
            hasLoop <- identical(oldClassif,nclassif)
            oldClassif <- classif
            classif <- nclassif
            error <- bmus$error
            if(verbose) {
                print(paste(i,j,error))
            }
            errors[[i]] <- c(errors[[i]],error)
            ## there is no representation phase!
            if(noChange | hasLoop | j==maxiter) {
                if(verbose) {
                    if(noChange) {
                        print(paste("radius:",radii[i],"iteration",j,
                                    "is stable, decreasing radius"))
                    } else {
                        print(paste("radius:",radii[i],"iteration",j,
                                    "oscillation detected, decreasing radius"))
                    }
                }
                if(i==length(radii)) {
                    ## the fitting is done
                    break;
                }
                ## preparing the loop with the next radius
                nv <- neighborhood(somgrid,radii[i+1],kernel,
                                   normalised=normalised)
            }
            ## break the loop if the partition is stable or when we have
            ## an oscillating behaviour
            if(noChange || hasLoop) {
                break;
            }
        }
        if(!noChange && verbose) {
            print(paste("warning: can't reach a stable configuration with radius",i))
        }
    }
    if(!noChange) {
        ## final assignment
        if(is.null(weights)) {
            bmus <-  fastKernelsombmu(classif,somgrid$size,K,nv)
        } else {
            bmus <-  weightedFastKernelsombmu(classif,somgrid$size,K,nv,weights)
        }
        classif <- bmus$clusters
        errors[[length(radii)]] <- c(errors[[length(radii)]],bmus$error)
    }
    ## for latter use
### FIXME: shouldn't that be moved before the last assignment?
    if(is.null(weights)) {
        nvcl <- nv[,classif]
    } else {
        nvcl <- sweep(nv[,classif],2,weights,"*")
    }
    normed <- rowSums(nvcl)
    prototypes <- sweep(nvcl,1,normed,"/")
    Kp <- tcrossprod(K,prototypes)
    pnorms <- double(nrow(prototypes))
    for(i in 1:length(pnorms)) {
        pnorms[i] <- c(prototypes[i,]%*%Kp[,i])
    }
    res <- list(somgrid=somgrid,prototypes=prototypes,classif=classif,
                errors=unlist(errors),Kp=Kp,pnorms=pnorms)
    class(res) <- c("kernelsom","som")
    res
}

kernelsombmu <- function(prototypes,K,weights) {
    Kp <- tcrossprod(K,prototypes)
    pnorms <- double(nrow(prototypes))
    for(i in 1:length(pnorms)) {
        pnorms[i] <- c(prototypes[i,]%*%Kp[,i])
    }
    predistances <- sweep(-2*Kp,2,pnorms,"+")
    bmu <- apply(predistances,1,which.min)
    if(missing(weights) || is.null(weights)) {
        error <- mean(predistances[cbind(1:length(bmu),bmu)]+diag(K))
    } else {
        error <- sum(weights*(predistances[cbind(1:length(bmu),bmu)]+diag(K)))/sum(weights)
    }
    list(clusters=bmu,error=error,Kp=Kp,pnorms=pnorms)
}

kernelsom.lowlevel.R <- function(somgrid,K,prototypes,
                                 assignment,radii,weights,maxiter,
                                 kernel,normalised,cut,verbose) {
    oldClassif <- rep(NA,nrow(K))
    classif <- rep(NA,nrow(K))
    errors <- vector("list",length(radii))
    if(missing(weights) || is.null(weights)) {
        weights <- rep(1,nrow(K))
    }
    nv <- neighborhood(somgrid,radii[1],kernel,normalised=normalised)
    for(i in 1:length(radii)) {
        for(j in 1:maxiter) {
            ## assignment
            if(assignment == "single") {
                bmus <- kernelsombmu(prototypes,K,weights)
            } else {
                stop(paste(assignment,"is not implemented for kernel SOM"))
            }
            nclassif <- bmus$clusters
            noChange <- identical(classif,nclassif)
            hasLoop <- identical(oldClassif,nclassif)
            oldClassif <- classif
            classif <- nclassif
            error <- bmus$error
            if(verbose) {
                print(paste(i,j,error))
            }
            errors[[i]] <- c(errors[[i]],error)
            ## there is no representation phase!
            if(noChange | hasLoop | j==maxiter) {
                if(verbose) {
                    if(noChange) {
                        print(paste("radius:",radii[i],"iteration",j,
                                    "is stable, decreasing radius"))
                    } else {
                        print(paste("radius:",radii[i],"iteration",j,
                                    "oscillation detected, decreasing radius"))
                    }
                }
                if(i==length(radii)) {
                    ## the fitting is done
                    break;
                }
                ## preparing the loop with the next radius
                nv <- neighborhood(somgrid,radii[i+1],kernel,
                                   normalised=normalised)
            }
            ## always update the weights
            nvcl <-  sweep(nv[,classif],2,weights,"*")
            normed <- rowSums(nvcl)
            prototypes <- sweep(nvcl,1,normed,"/")
            ## break the loop if the partition is stable or when we have
            ## an oscillating behaviour
            if(noChange || hasLoop) {
                break;
            }
        }
        if(!noChange && verbose) {
            print(paste("warning: can't reach a stable configuration with radius",i))
        }
    }
    if(!noChange) {
        ## final assignment
        bmus <- kernelsombmu(prototypes,K,weights)
        classif <- bmus$clusters
        errors[[length(radii)]] <- c(errors[[length(radii)]],bmus$error)
    }
    ## for latter use
    if(is.null(weights)) {
        nvcl <- nv[,classif]
    } else {
        nvcl <- sweep(nv[,classif],2,weights,"*")
    }
    res <- list(somgrid=somgrid,prototypes=prototypes,classif=classif,
                errors=unlist(errors),Kp=bmus$Kp,pnorms=bmus$pnorms)
    class(res) <- c("kernelsom","som")
    res
}

batchsom.kernelmatrix <- function(data,somgrid,init=c("pca","random"),
                                  prototypes,
                                  assignment=c("single","heskes"),
                                  radii=somradii(somgrid),
                                  weights,maxiter=75,
                                  kernel=c("gaussian","linear"),normalised,
                                  cut=1e-7,verbose=FALSE,keepdata=TRUE,...) {
    ## process parameters and perform a few sanity checks
    if(verbose) {
        print(match.call())
    }
    assignment <- match.arg(assignment)
    if(missing(normalised)) {
        normalised <- assignment=="heskes"
    }
    kernel <- match.arg(kernel)
    theKernel <- switch(kernel,"gaussian"=kernel.gaussian,"linear"=kernel.linear)
    if(class(somgrid)!="somgrid") {
        stop("'somgrid' is not of somgrid class")
    }
    if(!missing(weights)) {
        if(length(weights)!=nrow(data)) {
            stop("'weights' and 'data' have different dimensions")
        }
    } else {
        weights <- NULL
    }
    if(missing(prototypes)) {
        ## initialisation based on the value of init
        init <- match.arg(init)
        args <- list(...)
        params <- c(list(data=data,somgrid=somgrid,weights=weights),args)
        if(init=="random") {
            prototypes <- do.call("sominit.random",params)
        } else {
            initresults <- do.call("sominit.pca",params)
            prototypes <- initresults$prototypes
        }
    } else {
        if(ncol(prototypes)!=ncol(data)) {
            stop("'prototypes' and 'data' have different dimensions")
        }
        if(nrow(prototypes)!=somgrid$size) {
            stop("'prototypes' and 'somgrid' are not compatible")
        }
    }
    ## distances?
    if(is.null(somgrid$dist)) {
        somgrid$dist <- as.matrix(dist(somgrid$pts,method="Euclidean"),diag=0)
    }
    pre <- fastKernelsom.lowlevel.R(somgrid,data,prototypes,assignment,
                                    radii,weights,maxiter,theKernel,
                                    normalised,cut,verbose)
    pre$assignment <- assignment
    pre$kernel <- kernel
    pre$normalised <- normalised
    pre$radii <- radii
    if(keepdata) {
        pre$data  <- data
        pre$weights <- weights
    }
    pre
}

sominit.random.kernelmatrix <- function(data,somgrid,
                                        method=c("prototypes","random","cluster"),weights,...) {
### FIXME: share this with relationalsom
    method <- match.arg(method)
    nb.data <- nrow(data)
    nb <- somgrid$size
    if(method=="prototypes" || (method=="cluster" && nb>=nb.data)) {
        protos <- matrix(0,ncol=nb.data,nrow=nb)
        if(missing(weights) || is.null(weights)) {
            protos[cbind(1:nb,sample(nb.data,size=nb,replace=nb>nb.data))] <- 1
        } else {
            protos[cbind(1:nb,sample(nb.data,size=nb,replace=nb>nb.data,prob=weights))] <- 1
        }
    } else if(method=="random") {
        protos <- matrix(runif(nb.data*nb),ncol=nb.data,nrow=nb)
        protos <- sweep(protos,1,rowSums(protos),"/")
    } else {
        ## nb <nb.data
        clusters <- cut(sample(nb.data),nb,labels=FALSE,include.lowest=TRUE)
        protos <- matrix(0,ncol=nb.data,nrow=nb)
        protos[cbind(clusters,1:nb.data)] <- 1
        protos <- sweep(protos,1,rowSums(protos),"/")
    }
    protos
}

sominit.pca.kernelmatrix <- function(data, somgrid, ...) {
### FIXME: data weights support
    ## we do something very close to MDS and kernel PCA
    ## first double centering
    D.c <- double.centering(data)
    ## then eigenanalysis
    D.eigen <- eigen(D.c,symmetric=T)
    ## normalize (positive eigen values only)
    positive <- D.eigen$values>0
### FIXME: check that we indeed have positive eigne values!
    D.eigen$vectors[,positive] <- sweep(D.eigen$vectors[,positive],2,sqrt(D.eigen$values[positive]),"/")
    ## and center
    D.eigen$vectors <- scale(D.eigen$vectors,scale=FALSE)
    ## compute standard deviations
    sdev <- sqrt(D.eigen$values[positive]/nrow(data))

    ## the more detailled axis is assigned to the axis with the largest
    ## standard deviation
    if (somgrid$xdim>=somgrid$ydim) {
        x.ev <- 1
        y.ev <- 2
    } else {
        x.ev <- 2
        y.ev <- 1
    }
    if(somgrid$topo=="hexagonal") {
        xspan <- somgrid$xdim - 1
        if(somgrid$ydim>1) {
            xspan <- xspan+0.5
        }
        x <- seq(from=-2*sdev[x.ev],by=4*sdev[x.ev]/xspan,length.out=somgrid$xdim)
    } else {
        x <- seq(from=-2*sdev[x.ev],to=2*sdev[x.ev],length.out=somgrid$xdim)
    }
    y <- seq(from=2*sdev[y.ev],to=-2*sdev[y.ev],length.out=somgrid$ydim)
    base <- as.matrix(expand.grid(x = x, y = y))
    ## correction for hexagonal grids
    if(somgrid$topo=="hexagonal") {
        base[,1] <- base[,1]+rep(c(0,2*sdev[x.ev]/xspan),each=somgrid$xdim,length.out=nrow(base))
    }
    ## map back the grid to the dissimilarity space
    ## and centers at the same time
    prototypes <- tcrossprod(base,D.eigen$vectors[,c(x.ev,y.ev)])+1/nrow(data)
    list(prototypes=prototypes,D.c=D.c,D.eigen=D.eigen,sdev=sdev)
}

summary.kernelsom <- function(object,...)
{
    cat("\nKernel Self-Organising Map\n\n")
    NextMethod()
}

print.kernelsom <- function(x,...)
{
    cat("\nKernel Self-Organising Map\n\n")
    NextMethod()
}


as.matrix.kernelsom <- function(x,...) {
### FIXME: check for negative values?
    predist <- -2*x$prototypes%*%x$Kp
    thedist <- sweep(predist,1,x$pnorms,"+")
    thedist <- sweep(thedist,2,x$pnorms,"+")
    sqrt(thedist)
}

as.dist.kernelsom <- function(x,FUN=NULL) {
    as.dist(as.matrix(x))
}
