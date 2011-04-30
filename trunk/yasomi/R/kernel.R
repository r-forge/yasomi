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


fastKernelsom.lowlevel.R <- function(somgrid,K,prototypes,weights,control,
                                     verbose) {
    oldClassif <- rep(NA,nrow(K))
    errors <- vector("list",length(control$radii))
    ## a round of initialisation is needed
    nv <- neighborhood(somgrid,control$radii[1],control$kernel.fun,normalised=control$normalised)
    bmus <- kernelsombmu(prototypes,K,weights)
    classif <- bmus$clusters
    errors[[1]] <- bmus$error
    if(verbose) {
        print(paste(1,1,bmus$error))
    }
    for(i in 1:length(control$radii)) {
        if(i==1) {
            iterations <- 2:control$max.iter
        } else {
            iterations <- 1:control$max.iter
        }
        for(j in iterations) {
            ## assignment
            if(control$assignment == "single") {
                if(is.null(weights)) {
                    bmus <-  fastKernelsombmu(classif,somgrid$size,K,nv)
                } else {
                    bmus <-  weightedFastKernelsombmu(classif,somgrid$size,K,nv,weights)
                }
            } else {
                stop(paste(control$assignment,"is not implemented for kernel SOM"))
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
            if(noChange | hasLoop | j==control$max.iter) {
                if(verbose) {
                    if(noChange) {
                        print(paste("radius:",control$radii[i],"iteration",j,
                                    "is stable, decreasing radius"))
                    } else {
                        print(paste("radius:",control$radii[i],"iteration",j,
                                    "oscillation detected, decreasing radius"))
                    }
                }
                if(i==length(control$radii)) {
                    ## the fitting is done
                    break;
                }
                ## preparing the loop with the next radius
                nv <- neighborhood(somgrid,control$radii[i+1],
                                   control$kernel.fun,
                                   normalised=control$normalised)
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
        errors[[length(control$radii)]] <- c(errors[[length(control$radii)]],bmus$error)
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

fastKernelsom.lowlevelcontinuous.R <-
    function(somgrid,K,prototypes,weights,control,verbose) {
    oldClassif <- rep(NA,nrow(K))
    errors <- rep(NA,length(control$radii))
    ## a round of initialisation is needed
    bmus <- kernelsombmu(prototypes,K,weights)
    classif <- bmus$clusters
    errors[1] <- bmus$error
    if(verbose) {
        print(paste(1,bmus$error))
    }
    ## no representation phase is needed
    nv <- neighborhood(somgrid,control$radii[1],control$kernel.fun,normalised=control$normalised)
    for(i in 2:length(control$radii)) {
        ## assignment
        if(control$assignment == "single") {
            if(is.null(weights)) {
                bmus <-  fastKernelsombmu(classif,somgrid$size,K,nv)
            } else {
                bmus <-  weightedFastKernelsombmu(classif,somgrid$size,K,nv,weights)
            }
        } else {
            stop(paste(control$assignment,"is not implemented for kernel SOM"))
        }
        nclassif <- bmus$clusters
        noChange <- identical(classif,nclassif)
        classif <- nclassif
        errors[i] <- bmus$error
        if(verbose) {
            print(paste(i,errors[i]))
        }
        nv <- neighborhood(somgrid,control$radii[i],
                           control$kernel.fun,
                           normalised=control$normalised)
    }
    ## for latter use
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
    if(!noChange) {
        ## final assignment
        if(is.null(weights)) {
            bmus <-  fastKernelsombmu(classif,somgrid$size,K,nv)
        } else {
            bmus <-  weightedFastKernelsombmu(classif,somgrid$size,K,nv,weights)
        }
        classif <- bmus$clusters
        errors <- c(errors,bmus$error)
    }
    res <- list(somgrid=somgrid,prototypes=prototypes,classif=classif,
                errors=errors,Kp=Kp,pnorms=pnorms)
    class(res) <- c("kernelsom","som")
    res
}

kernelsombmu <- function(prototypes,K,weights) {
###FIXME: rename and share with kmeans
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
                                  weights,
                                  mode = c("continuous","stepwise"),
                                  min.radius, max.radius, steps,
                                  decrease = c("power", "linear"), max.iter,
                                  kernel = c("gaussian", "linear", "zeroone"),
                                  normalised,
                                  assignment = c("single", "heskes"),
                                  cut = 1e-07,
                                  verbose=FALSE,keepdata=TRUE,...) {
    ## process parameters and perform a few sanity checks
    if(class(somgrid)!="somgrid") {
        stop("'somgrid' is not of somgrid class")
    }
    the.call <- match.call()
    if(verbose) {
        print(the.call)
    }
    the.call[[1]] <- batchsom.control
    control <- eval(the.call,envir = parent.frame())
    control$kernel.fun <- switch(control$kernel,"gaussian"=kernel.gaussian,"linear"=kernel.linear,"zeroone"=kernel.zeroone)
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
    pre <- switch(control$mode,
                  "stepwise"=fastKernelsom.lowlevel.R(somgrid,data,
                  prototypes,weights,control,verbose),
                  "continuous"=fastKernelsom.lowlevelcontinuous.R(somgrid,data,
                  prototypes,weights,control,verbose))
    pre$control <- control
    if(keepdata) {
        pre$data  <- data
        pre$weights <- weights
    }
    pre
}

sominit.random.kernelmatrix <- function(data,somgrid,
                                        method=c("prototypes","random","cluster"),weights,...) {
    method <- match.arg(method)
    if(missing(weights)) {
        weights <- NULL
    }
    convex.prototypes.random(data,somgrid$size,method,weights)
}

sominit.pca.kernelmatrix <- function(data, somgrid, ...) {
### FIXME: data weights support
    ## we do something very close to MDS and kernel PCA
    ## first double centering
    D.c <- double.centering(data)
    ## then eigenanalysis
    D.eigen <- eigen(D.c,symmetric=T)
    min.eigenvalue <- min(D.eigen$values)
    if(min.eigenvalue < -sqrt(.Machine$double.eps)) {
        warning(paste("Non negligible negative eigenvalue ",min.eigenvalue))
    }
    ## normalize (positive eigen values only)
    positive <- D.eigen$values>0
    if(sum(positive)<2) {
        stop("The centered kernel matrix must have at least two positive eigenvalues")
    }
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
