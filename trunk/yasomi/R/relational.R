# relational SOM

relational.keeppositive <- function(val,with.warning=TRUE) {
    valmin <- min(val)
    if(with.warning && valmin < -sqrt(.Machine$double.eps)) {
        warning(paste("Non negligible negative minimal values (",valmin,") discared",sep=""))
    }
    val[val<0] <- 0
    val
}

relational.sqrt <- function(val,with.warning=TRUE) {
    sqrt(relational.keeppositive(val,with.warning))
}

sominit.random.dist <- function(data,somgrid,
                                method=c("prototypes","random","cluster"),weights,...) {
    method <- match.arg(method)
    if(missing(weights)) {
        weights <- NULL
    }
    convex.prototypes.random(data,somgrid$size,method,weights)
}

double.centering <- function(D) {
    D.rowmean <- rowMeans(D)
    D.mean <- mean(D.rowmean)
    sweep(sweep(D,1,D.rowmean,"-"),2,D.rowmean,"-")+D.mean
}

normsFromDist <- function(D) {
    D.rowmean <- rowMeans(D)
    D.mean <- mean(D.rowmean)
    -0.5*(diag(D)+D.mean-2*D.rowmean)
}

sominit.pca.dist <- function(data, somgrid, ...) {
### FIXME: data weights support
    D <- as.matrix(data^2,diag=0)
    ## we do something very close to MDS and kernel PCA
    ## first double centering
    D.c <- -0.5*double.centering(D)
    ## then eigenanalysis
    D.eigen <- eigen(D.c,symmetric=T)
    ## normalize (positive eigen values only)
    positive <- D.eigen$values>0
    D.eigen$vectors[,positive] <- sweep(D.eigen$vectors[,positive],2,sqrt(D.eigen$values[positive]),"/")
    ## compute standard deviations
    sdev <- sqrt(D.eigen$values[positive]/nrow(D))

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
    prototypes <- tcrossprod(base,D.eigen$vectors[,c(x.ev,y.ev)])
    list(prototypes=prototypes,D=D,D.c=D.c,D.eigen=D.eigen,sdev=sdev)
}

fastRelationalBMU.R <- function(cluster,nclust,diss,nv) {
    ps <- partialSums(cluster,nclust,diss)
    csize <- table(factor(cluster,levels=1:nclust))
    normed <- nv%*%csize
    nvnormed <- sweep(nv,1,normed,"/")
    Dalpha <- nvnormed%*%ps$ps
    interm <- tcrossprod(ps$bips,nvnormed)
    nf <- double(nclust)
    for(i in 1:nclust) {
        nf[i] <- 0.5*c(nvnormed[i,]%*%interm[,i])
    }
    ## the alternate C version is barely faster
    ## nf <- 0.5*nfPS(ps$bips,nvnormed,nclust)
    distances <- sweep(Dalpha,1,nf,"-")
    bmu <- apply(distances,2,which.min)
    error <- mean(relational.keeppositive(distances[cbind(bmu,1:length(bmu))]))
    list(clusters=bmu,error=error,Dalpha=Dalpha,nf=nf)
}

weightedFastRelationalBMU.R <- function(cluster,nclust,diss,nv,weights) {
    ps <- weightedPartialSums(cluster,nclust,diss,weights)
    csize <- tapply(weights,factor(cluster,levels=1:nclust),sum)
    csize[is.na(csize)] <- 0
    normed <- nv%*%csize
    nvnormed <- sweep(nv,1,normed,"/")
    Dalpha <- nvnormed%*%ps$ps
    interm <- tcrossprod(ps$bips,nvnormed)
    nf <- double(nclust)
    for(i in 1:nclust) {
        nf[i] <- 0.5*c(nvnormed[i,]%*%interm[,i])
    }
    ## the alternate C version is barely faster
    ## nf <- 0.5*nfPS(ps$bips,nvnormed,nclust)
    distances <- sweep(Dalpha,1,nf,"-")
    bmu <- apply(distances,2,which.min)
    error <- sum(weights*relational.keeppositive(distances[cbind(bmu,1:length(bmu))]))/sum(weights)
    list(clusters=bmu,error=error,Dalpha=Dalpha,nf=nf)
}

partialSums <- function(cluster,nclust,diss) {
    datasize <- as.integer(length(cluster))
    ## partial_sums modifies only the bisums parameter: DUP=FALSE is safe here
    result <- .C("partial_sums",
                 as.integer(cluster-1),
                 datasize,
                 as.integer(nclust),
                 as.double(diss),
                 sums=double(nclust*datasize),
                 bisums=double(nclust^2),DUP=FALSE)
    list(ps=matrix(result$sums,nrow=nclust,ncol=datasize),
         bips=matrix(result$bisums,nrow=nclust,ncol=nclust))
}

weightedPartialSums <- function(cluster,nclust,diss,weights) {
    datasize <- as.integer(length(cluster))
    ## partial_sums modifies only the bisums parameter: DUP=FALSE is safe here
    result <- .C("weighted_partial_sums",
                 as.integer(cluster-1),
                 datasize,
                 as.integer(nclust),
                 as.double(diss),
                 as.double(weights),
                 sums=double(nclust*datasize),
                 bisums=double(nclust^2),DUP=FALSE)
    list(ps=matrix(result$sums,nrow=nclust,ncol=datasize),
         bips=matrix(result$bisums,nrow=nclust,ncol=nclust))
}

nfPS <- function(bips,nvnormed,nclust) {
    ## th_bips_h modifies only the nf parameter: DUP=FALSE is safe here
    .C("th_bips_h",as.double(bips),as.double(nvnormed),as.integer(nclust),
       nf=double(nclust),DUP=FALSE)$nf
}

relationalbmu.R <- function(prototypes,diss,weights) {
    ## first compute the base distances
    Dalpha <- tcrossprod(diss,prototypes)
    ## then the normalisation factor
    ## can we do this faster?
    nf <- double(nrow(prototypes))
    for(i in 1:length(nf)) {
        nf[i] <- 0.5*c(prototypes[i,]%*%Dalpha[,i])
    }
    distances <- sweep(Dalpha,2,nf,"-")
    clusters <- apply(distances,1,which.min)
    if(missing(weights) || is.null(weights)) {
        error <- mean(relational.keeppositive(distances[cbind(1:length(clusters),clusters)]))
    } else {
        error <- sum(weights*relational.keeppositive(distances[cbind(1:length(clusters),clusters)]))/sum(weights)
    }
    list(clusters=clusters,error=error,Dalpha=Dalpha,nf=nf)
}

predict.relationalsom <- function(object,newdata,with.secondwinner=FALSE,...) {
    if(!inherits(newdata,"crossdist")) {
        stop("'newdata' is not a crossdist object")
    }
    if(nrow(newdata)!=ncol(object$prototypes)) {
        stop("'newdata' and 'object$prototypes' have different dimensions")
    }
    rdist <- sweep(object$prototypes%*%(newdata^2),1,object$nf,"-")
    if(with.secondwinner) {
### FIXME: very suboptimal
        ordered <- apply(rdist,2,order)
        winners <- t(ordered[1:2,])
        clusters <- winners[,1]
    } else {
        clusters <- apply(rdist,2,which.min)
    }
    error <- mean(relational.keeppositive(rdist[cbind(clusters,1:length(clusters))]))
    if(with.secondwinner) {
        list(classif=clusters,error=error,distances=rdist,winners=winners)
    } else {
        list(classif=clusters,error=error,distances=rdist)
    }
}


extended.relationalbmu.R <- function(prototypes,diss,norms,weights) {
    ## we use here the extended formula when rowSums(prototypes)!=1
    Dalpha <- tcrossprod(diss,prototypes)
    ## then the normalisation factor
    ## can we do this faster?
    nf <- double(nrow(prototypes))
    for(i in 1:length(nf)) {
        nf[i] <- 0.5*c(prototypes[i,]%*%Dalpha[,i])
    }
    sums <- 1-rowSums(prototypes)
    protonorms <- prototypes%*%norms
    distances <- sweep(Dalpha,2,nf+sums*protonorms,"-")+norms%o%sums
    clusters <- apply(distances,1,which.min)
    if(missing(weights) || is.null(weights)) {
        error <- mean(relational.keeppositive(distances[cbind(1:length(clusters),clusters)]))
    } else {
        error <- sum(weights*relational.keeppositive(distances[cbind(1:length(clusters),clusters)]))/sum(weights)
    }
    list(clusters=clusters,error=error,Dalpha=Dalpha,nf=nf)
}

relationalsecondbmu.R <- function(Dalpha,nf) {
    distances <- sweep(Dalpha,2,nf,"-")
    ## very suboptimal
    ordered <- apply(distances,1,order)
    winners <- t(ordered[1:2,])
    error <- relational.keeppositive(distances[cbind(1:nrow(winners),winners[,1])])
    list(winners=winners,error=error)
}

fastRelationalsom.lowlevelcontinuous.R <- function(somgrid,diss,prototypes,
                                                   weights,
                                                   control,verbose,extended,
                                                   data.norms) {
    oldClassif <- rep(NA,nrow(diss))
    errors <- rep(NA,length(control$radii))
    ## a round of initialisation is needed
    if(extended) {
        bmus <- extended.relationalbmu.R(prototypes,diss,data.norms,weights)
    } else {
        bmus <- relationalbmu.R(prototypes,diss,weights)
    }
    classif <- bmus$clusters
    errors[1] <- bmus$error
    if(verbose) {
        print(paste(1,bmus$error))
    }
    nv <- neighborhood(somgrid,control$radii[1],control$kernel.fun,normalised=control$normalised)
    for(i in 2:length(control$radii)) {
        ## assignment
        if(control$assignment == "single") {
            if(is.null(weights)) {
                bmus <- fastRelationalBMU.R(classif,somgrid$size,diss,nv)
            } else {
                bmus <- weightedFastRelationalBMU.R(classif,somgrid$size,diss,nv,weights)
            }
        } else {
            stop(paste(control$assignment,"is not implemented for relational SOM"))
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
    if(verbose) {
        print("computing final prototypes, Dalpha and nf")
    }
    prototypes <- sweep(nvcl,1,normed,"/")
    Dalpha <- tcrossprod(diss,prototypes)
    ## this is slow
    nf <- 0.5*diag(prototypes%*%Dalpha)
    if(!noChange) {
        ## final assignment
        if(is.null(weights)) {
            bmus <- fastRelationalBMU.R(classif,somgrid$size,diss,nv)
        } else {
            bmus <- weightedFastRelationalBMU.R(classif,somgrid$size,diss,nv,weights)
        }
        classif <- bmus$clusters
        errors <- c(errors,bmus$error)
    }
    res <- list(somgrid=somgrid,prototypes=prototypes,classif=classif,
                errors=errors,Dalpha=Dalpha,nf=nf)
    class(res) <- c("relationalsom","som")
    res
}

fastRelationalsom.lowlevel.R <- function(somgrid,diss,prototypes,weights,
                                         control,verbose,extended,
                                         data.norms) {
    oldClassif <- rep(NA,nrow(diss))
    errors <- vector("list",length(control$radii))
    nv <- neighborhood(somgrid,control$radii[1],control$kernel.fun,normalised=control$normalised)
    ## a round of initialisation is needed
    if(extended) {
        bmus <- extended.relationalbmu.R(prototypes,diss,data.norms,weights)
    } else {
        bmus <- relationalbmu.R(prototypes,diss,weights)
    }
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
                    bmus <- fastRelationalBMU.R(classif,somgrid$size,diss,nv)
                } else {
                    bmus <- weightedFastRelationalBMU.R(classif,somgrid$size,diss,nv,weights)
                }
            } else {
                stop(paste(control$assignment,"is not implemented for relational SOM"))
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
            if(noChange | hasLoop | j == control$max.iter) {
                if(verbose) {
                    if(noChange) {
                        print(paste("radius:",control$radii[i],"iteration",j,
                                    "is stable, decreasing radius"))
                    } else if(hasLoop) {
                        print(paste("radius:",control$radii[i],"iteration",j,
                                    "oscillation detected, decreasing radius"))
                    }
                }
                if(i<length(control$radii)) {
                    ## preparing the loop with the next radius
                    nv <- neighborhood(somgrid,control$radii[i+1],
                                       control$kernel.fun,
                                       normalised=control$normalised)
                }
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
            bmus <- fastRelationalBMU.R(classif,somgrid$size,diss,nv)
        } else {
            bmus <- weightedFastRelationalBMU.R(classif,somgrid$size,diss,nv,weights)
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
    if(verbose) {
        print("computing final prototypes, Dalpha and nf")
    }
    prototypes <- sweep(nvcl,1,normed,"/")
    Dalpha <- tcrossprod(diss,prototypes)
    ## this is slow
    nf <- 0.5*diag(prototypes%*%Dalpha)
    res <- list(somgrid=somgrid,prototypes=prototypes,classif=classif,
                errors=unlist(errors),Dalpha=Dalpha,nf=nf)
    class(res) <- c("relationalsom","som")
    res
}

batchsom.dist <- function(data,somgrid,init=c("pca","random"),prototypes,
                          weights,
                          mode = c("continuous","stepwise"),
                          min.radius, max.radius, steps,
                          decrease = c("power", "linear"), max.iter,
                          kernel = c("gaussian", "linear"), normalised,
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
#    the.call[[1]] <- as.name("batchsom.control")
    the.call[[1]] <- batchsom.control
    control <- eval(the.call,envir = parent.frame())
    control$kernel.fun <- switch(control$kernel,"gaussian"=kernel.gaussian,"linear"=kernel.linear)
    if(!missing(weights)) {
        if(length(weights)!=nrow(data)) {
            stop("'weights' and 'data' have different dimensions")
        }
    } else {
        weights <- NULL
    }
    ## diss is initialized in this code
    data.norms <- NULL
    if(missing(prototypes)) {
        if(verbose) {
            print("Initializing prototypes")
        }
        ## initialisation based on the value of init
        init <- match.arg(init)
        args <- list(...)
        params <- c(list(data=data,somgrid=somgrid,weights=weights),args)
        if(init=="random") {
            prototypes <- do.call("sominit.random",params)
            extended <- FALSE
            diss <- as.matrix(data^2,diag=0)
        } else {
            initresults <- do.call("sominit.pca",params)
            prototypes <- initresults$prototypes
            extended <- TRUE
            data.norms <- diag(initresults$D.c)
            diss <- initresults$D
        }
    } else {
        diss <- as.matrix(data^2,diag=0)
        if(ncol(prototypes)!=ncol(diss)) {
            stop("'prototypes' and 'diss' have different dimensions")
        }
        if(nrow(prototypes)!=somgrid$size) {
            stop("'prototypes' and 'somgrid' are not compatible")
        }
        if(!isTRUE(all.equal(rowSums(prototypes),rep(1,nrow(prototypes))))) {
            if(verbose) {
                print("'prototypes' rows do not sum to one, using extended relational formula")
            }
            extended <- TRUE
            data.norms <- normsFromDist(diss)
        } else {
            extended <- FALSE
        }
    }
    ## distances?
    if(is.null(somgrid$dist)) {
        somgrid$dist <- as.matrix(dist(somgrid$pts,method="Euclidean"),diag=0)
    }
    pre <- switch(control$mode,
                  "stepwise"=fastRelationalsom.lowlevel.R(somgrid,diss,
                  prototypes,weights,control,verbose,extended,data.norms),
                  "continuous"=fastRelationalsom.lowlevelcontinuous.R(somgrid,
                  diss,prototypes,weights,control,verbose,extended,data.norms))
    pre$control <- control
    if(keepdata) {
        pre$data  <- data
        pre$weights <- weights
    }
    pre
}

summary.relationalsom <- function(object,...)
{
    cat("\nRelational Self-Organising Map\n\n")
    NextMethod()
}

print.relationalsom <- function(x,...)
{
    cat("\nRelational Self-Organising Map\n\n")
    NextMethod()
}


as.matrix.relationalsom <- function(x,...) {
    predist <- x$prototypes%*%x$Dalpha
    thedist <- sweep(predist,1,x$nf,"-")
    thedist <- sweep(thedist,2,x$nf,"-")
    relational.sqrt(thedist)
}

as.dist.relationalsom <- function(x,FUN=NULL) {
    as.dist(as.matrix(x))
}
