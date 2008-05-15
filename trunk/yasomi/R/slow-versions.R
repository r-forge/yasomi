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

bmu.R <- function(prototypes,data) {
    result=rep(0,nrow(data))
    for(i in 1:length(result)) {
        result[i] <- which.min(rowSums(sweep(prototypes,2,data[i,],"-")^2))
    }
    result
}

bmu.heskes.R <- function(prototypes,data,nv) {
    result=rep(0,nrow(data))
    for(i in 1:length(result)) {
        result[i] <- which.min(nv%*%rowSums(sweep(prototypes,2,data[i,],"-")^2))
    }
    result
}

batchsom.R <- function(data,somgrid,init=c("pca","random"),prototypes,
                       assignment=c("single","heskes"),radii=somradii(somgrid),
                       maxiter=75,
                       kernel=c("gaussian","linear"),normalised,
                       cut=1e-7,verbose=FALSE,keepdata=TRUE,...) {
    ## process parameters and perform a few sanity checks
    assignment <- match.arg(assignment)
    if(missing(normalised)) {
        normalised <- assignment=="heskes"
    }
    kernel <- match.arg(kernel)
    theKernel <- switch(kernel,"gaussian"=kernel.gaussian,"linear"=kernel.linear)
    if(class(somgrid)!="somgrid") {
        stop("'somgrid' is not of somgrid class")
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
    pre <- batchsom.lowlevel.R(somgrid,data,prototypes,assignment,radii,
                               maxiter,theKernel,normalised,cut,verbose)
    pre$assignment <- assignment
    pre$kernel <- kernel
    pre$normalised <- normalised
    pre$radii <- radii
    if(keepdata) {
        pre$data  <- data
    }
    pre
}

batchsom.lowlevel.R <- function(somgrid,data,prototypes,
                                assignment,radii,maxiter,kernel,
                                normalised,cut,verbose) {
    data <- as.matrix(data)
    classif <- rep(NA,nrow(data))
    errors <- vector("list",length(radii))
    for(i in 1:length(radii)) {
        nv <- neighborhood(somgrid,radii[i],kernel,normalised=normalised)
        for(j in 1:maxiter) {
            if(assignment == "single") {
                bmus <- bmu(prototypes,data)
            } else {
                bmus <- bmu.heskes(prototypes,data,nv)
            }
            nclassif <- bmus$clusters
            noChange = identical(classif,nclassif)
            classif <- nclassif
            error <- bmus$error
            if(verbose) {
                print(paste(i,j,error))
            }
            errors[[i]] <- c(errors[[i]],error)
            weights <- nv[,classif]
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

