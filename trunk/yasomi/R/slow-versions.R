neighborhood <- function(somgrid,T,kernel,normalised=TRUE) {
    raw <- kernel(somgrid$dist/T)
    if(normalised) {
        sweep(raw,1,rowSums(raw),"/")
    } else {
        raw
    }
}

kernel.gaussian <- function(x) {
    exp(-(3*x)^2/2)
}

kernel.linear <- function(x) {
    pre <- sapply(x,function(x) {
        if(x<1) {1-x} else {0}
    })
    if(is.matrix(x)) {
        matrix(pre,ncol=ncol(x),nrow=nrow(x))
    } else {
        pre
    }
}

bmu.R <- function(prototypes,data,weights=NULL) {
    distances <- dist(prototypes,data)
    clusters <- apply(distances,2,which.min)
    if(is.null(weights)) {
        error <- mean(distances[cbind(clusters,1:length(clusters))])
    } else {
        error <- sum(weights*distances[cbind(clusters,1:length(clusters))])/sum(weights)
    }
    list(clusters=clusters,error=error)
}

bmu.heskes.R <- function(prototypes,data,nv,weights=NULL) {
    distances <- dist(prototypes,data)
    dnv <- nv%*%(distances^2)
    clusters <- apply(dnv,2,which.min)
    if(is.null(weights)) {
        error <- mean(distances[cbind(clusters,1:length(clusters))])
    } else {
        error <- sum(weights*distances[cbind(clusters,1:length(clusters))])/sum(weights)
    }
    list(clusters=clusters,error=error)
}

batchsom.R <- function(data,somgrid,init=c("pca","random"),prototypes,
                       assignment=c("single","heskes"),radii=somradii(somgrid),
                       weights,maxiter=75,
                       kernel=c("gaussian","linear"),normalised,
                       cut=1e-7,verbose=FALSE,keepdata=TRUE,...) {
    ## process parameters and perform a few sanity checks
    if(class(somgrid)!="somgrid") {
        stop("'somgrid' is not of somgrid class")
    }
    assignment <- match.arg(assignment)
    if(missing(normalised)) {
        normalised <- assignment=="heskes"
    }
    kernel <- match.arg(kernel)
    theKernel <- switch(kernel,"gaussian"=kernel.gaussian,"linear"=kernel.linear)
    if(!missing(weights)) {
        if(length(weights)!=nrow(data)) {
            stop("'weights' and 'data' have different dimensions")
        }
    } else {
        weights=NULL
    }
    if(missing(prototypes)) {
        ## initialisation based on the value of init
        init <- match.arg(init)
        args <- list(...)
        params <- c(list(data=data,somgrid=somgrid),list(...))
        prototypes <- switch(init,
                             "pca"=do.call("sominit.pca",params)$prototypes,
                             "random"=do.call("sominit.random",params))
    } else {
        if(ncol(prototypes)!=ncol(data)) {
            stop("'prototypes' and 'data' have different dimensions")
        }
    }
    ## distances?
    if(is.null(somgrid$dist)) {
        somgrid$dist <- as.matrix(dist(somgrid$pts,method="Euclidean"),diag=0)
    }
    pre <- batchsom.lowlevel.R(somgrid,data,weights,prototypes,assignment,
                               radii,maxiter,theKernel,normalised,cut,verbose)
    pre$assignment <- assignment
    pre$kernel <- kernel
    pre$normalised <- normalised
    pre$radii <- radii
    if(keepdata) {
        pre$data  <- data
    }
    pre
}

batchsom.lowlevel.R <- function(somgrid,data,dataweights,prototypes,
                                assignment,radii,maxiter,kernel,
                                normalised,cut,verbose) {
    data <- as.matrix(data)
    classif <- rep(NA,nrow(data))
    errors <- vector("list",length(radii))
    for(i in 1:length(radii)) {
        nv <- neighborhood(somgrid,radii[i],kernel,normalised=normalised)
        for(j in 1:maxiter) {
            if(assignment == "single") {
                bmus <- bmu.R(prototypes,data,dataweights)
            } else {
                bmus <- bmu.heskes.R(prototypes,data,nv,dataweights)
            }
            nclassif <- bmus$clusters
            noChange = identical(classif,nclassif)
            classif <- nclassif
            error <- bmus$error
            if(verbose) {
                print(paste(i,j,error))
            }
            errors[[i]] <- c(errors[[i]],error)
            if(is.null(dataweights)) {
                weights <- nv[,classif]
            } else {
                weights <- sweep(nv[,classif],2,dataweights,"*")
            }
            normed <- rowSums(weights)
            mask <- (1:length(normed))[(normed/length(normed))>cut]
            prototypes[mask,] <- sweep(weights%*%data,1,normed,"/")[mask,]
            if(noChange && j >=2) {
                if(verbose) {
                    print(paste("radius:",radii[i],"iteration",j,"is stable, decreasing radius"))
                }
                break;
            }
        }
        if(!noChange && verbose) {
            print(paste("warning: can't reach a stable configuration with radius",i))
        }
    }
    res <- list(somgrid=somgrid,prototypes=prototypes,classif=classif,
                errors=unlist(errors))
    class(res) <- c("somnum","som")
    res
}

