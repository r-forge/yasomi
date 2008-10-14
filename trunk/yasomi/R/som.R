sominit.pca.default <- function(data,somgrid,...) {
### FIXME: data weights support
    ## we don't scale the data
    data.pca <- prcomp(data)
    ## the more detailled axis is assigned to the first eigenvector
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
        x <- seq(from=-2*data.pca$sdev[x.ev],by=4*data.pca$sdev[x.ev]/xspan,length.out=somgrid$xdim)
    } else {
        x <- seq(from=-2*data.pca$sdev[x.ev],to=2*data.pca$sdev[x.ev],length.out=somgrid$xdim)
    }
    y <- seq(from=2*data.pca$sdev[y.ev],to=-2*data.pca$sdev[y.ev],length.out=somgrid$ydim)
    base <- as.matrix(expand.grid(x = x, y = y))
    ## correction for hexagonal grids
    if(somgrid$topo=="hexagonal") {
        base[,1] <- base[,1]+rep(c(0,2*data.pca$sdev[x.ev]/xspan),each=somgrid$xdim,length.out=nrow(base))
    }
    ## map back the grid to the original space
    mapped <- tcrossprod(base,data.pca$rotation[,c(x.ev,y.ev)])
    ## decentering
    prototypes <- sweep(mapped,2,data.pca$center,"+")
    list(prototypes=prototypes,data.pca=data.pca)
}

sominit.random.default <- function(data,somgrid,
                                   method=c("prototypes","random","cluster"),
                                   weights,...) {
    method <- match.arg(method)
    data <- as.matrix(data)
    nb.data <- nrow(data)
    nb <- somgrid$size
    if(method=="prototypes" || (method=="cluster" && nb>=nb.data)) {
        if(missing(weights) || is.null(weights)) {
            data[sample(nb.data,size=somgrid$size,replace=nb>nb.data),]
        } else {
            data[sample(nb.data,size=somgrid$size,replace=nb>nb.data,prob=weights),]
        }
    } else if(method=="random") {
        result <- matrix(0,ncol=ncol(data),nrow=nb)
        for(i in 1:ncol(data)) {
            therange <- range(data[,i])
            result[,i] <- runif(nb,min=therange[1],max=therange[2])
        }
        result
    } else {
        ## nb <nb.data
        clusters <- cut(sample(nb.data),nb,labels=FALSE,include.lowest=TRUE)
        result <- matrix(0,ncol=ncol(data),nrow=nb)
        for(i in 1:nb) {
            result[i,] <- colMeans(as.matrix(data[clusters==i,],ncol=ncol(data)))
        }
        result
    }
}

bmu <- function(prototypes,data,weights) {
    if(ncol(prototypes)!=ncol(data)) {
        stop("'prototypes' and 'data' have different dimensions")
    }
    if(missing(weights)) {
        weights <- rep(1,nrow(data))
    } else if(length(weights)!=nrow(data)) {
        stop("'weights' and 'data' have different dimensions")
    }
    result <- .C("bmu",
                 as.double(prototypes),
                 as.integer(nrow(prototypes)),
                 as.double(data),
                 as.integer(nrow(data)),
                 as.integer(ncol(prototypes)),
                 as.double(weights),
                 clusters=integer(nrow(data)),
                 error=as.double(0),
                 PACKAGE="yasomi")
    list(clusters=result$clusters+1,error=result$error)
}

secondBmu <- function(prototypes,data) {
    if(ncol(prototypes)!=ncol(data)) {
        stop("'prototypes' and 'data' have different dimensions")
    }
    matrix(.C("second_bmu",
              as.double(prototypes),
              as.integer(nrow(prototypes)),
              as.double(data),
              as.integer(nrow(data)),
              as.integer(ncol(prototypes)),
              clusters=integer(2*nrow(data)),
              PACKAGE="yasomi")$clusters,ncol=2)
}

bmu.heskes <- function(prototypes,data,nv,weights) {
    if(ncol(prototypes)!=ncol(data)) {
        stop("'prototypes' and 'data' have different dimensions")
    }
    if(ncol(nv)!=nrow(nv)) {
        stop("'nv' is not a square matrix")
    }
    if(ncol(nv)!=nrow(prototypes)) {
        stop("'nv' and 'prototypes' have different dimensions")
    }
    if(missing(weights)) {
        weights <- rep(1,nrow(data))
    } else if(length(weights)!=nrow(data)) {
        stop("'weights' and 'data' have different dimensions")
    }
    ## nv must be in row major mode if normalised
    result <- .C("bmu_heskes",
                 as.double(prototypes),
                 as.double(t(nv)),
                 as.integer(nrow(prototypes)),
                 as.double(data),
                 as.integer(nrow(data)),
                 as.integer(ncol(prototypes)),
                 as.double(weights),
                 clusters=integer(nrow(data)),
                 error=as.double(0),
                 PACKAGE="yasomi")
    list(clusters=result$clusters+1,error=result$error)
}

batchsom.default <- function(data,somgrid,init=c("pca","random"),prototypes,
                     assignment=c("single","heskes"),radii=somradii(somgrid),
                     weights,maxiter=75,
                     kernel=c("gaussian","linear"),
                     normalised,
                     cut=1e-7,verbose=FALSE,keepdata=TRUE,...) {
    ## process parameters and perform a few sanity checks
    assignment <- match.arg(assignment)
    if(missing(normalised)) {
        normalised <- assignment=="heskes"
    }
    theRule <- switch(assignment,"single"=0,"heskes"=1)
    kernel <- match.arg(kernel)
    kernelType <- switch(kernel,"gaussian"=0,"linear"=1)
    if(class(somgrid)!="somgrid") {
        stop("'somgrid' is not of somgrid class")
    }
    if(!missing(weights)) {
        if(length(weights)!=nrow(data)) {
            stop("'weights' and 'data' have different dimensions")
        }
    } else {
        weights <- rep(1,nrow(data))
    }
    if(missing(prototypes)) {
        ## initialisation based on the value of init
        init <- match.arg(init)
        args <- list(...)
        params <- c(list(data=data,somgrid=somgrid,weights=weights),args)
        prototypes <- switch(init,
                             "pca"=do.call("sominit.pca",params)$prototypes,
                             "random"=do.call("sominit.random",params))
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
    pre <- batchsom.lowlevel(somgrid,data,prototypes,theRule,radii,weights,
                             maxiter,kernelType,normalised,cut,verbose)
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

newbatchsom.default <- function(data,somgrid,init=c("pca","random"),prototypes,
                                weights,
                                mode = c("stepwise","continuous"),
                                min.radius, max.radius, steps,
                                decrease = c("power", "linear"), max.iter,
                                kernel = c("gaussian", "linear"), normalised,
                                assignment = c("single", "heskes"),
                                cut = 1e-07,
                                verbose=FALSE,keepdata=TRUE,...) {
    if(class(somgrid)!="somgrid") {
        stop("'somgrid' is not of somgrid class")
    }
    the.call <- match.call()
    the.call[[1]] <- as.name("batchsom.control")
    control <- eval(the.call)
    theRule <- switch(control$assignment,"single"=0,"heskes"=1)
    kernelType <- switch(control$kernel,"gaussian"=0,"linear"=1)
    if(!missing(weights)) {
        if(length(weights)!=nrow(data)) {
            stop("'weights' and 'data' have different dimensions")
        }
    } else {
        weights <- rep(1,nrow(data))
    }
    if(missing(prototypes)) {
        ## initialisation based on the value of init
        init <- match.arg(init)
        args <- list(...)
        params <- c(list(data=data,somgrid=somgrid,weights=weights),args)
        prototypes <- switch(init,
                             "pca"=do.call("sominit.pca",params)$prototypes,
                             "random"=do.call("sominit.random",params))
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
                  "stepwise"=batchsom.lowlevel(somgrid,data,prototypes,theRule,
                  control$radii,weights,control$max.iter,kernelType,
                  control$normalised,control$cut,verbose),
                  "continuous"=batchsom.lowlevelcontinuous(somgrid,data,
                  prototypes,theRule,control$radii,weights,kernelType,
                  control$normalised,control$cut,verbose))
    pre$control <- control
    if(keepdata) {
        pre$data  <- data
        pre$weights <- weights
    }
    pre
}

batchsom.lowlevel <- function(somgrid,data,prototypes,
                              assignment,radii,weights,maxiter,
                              kernelType,normalised,cut,verbose) {
    result <- .C("batch_som_optim",
                 proto=as.double(prototypes),
                 as.integer(somgrid$size),
                 as.double(data),
                 as.integer(nrow(data)),
                 as.integer(ncol(data)),
                 as.double(weights),
                 as.integer(assignment),
                 as.double(somgrid$dist),
                 as.double(radii),
                 as.integer(length(radii)),
                 as.integer(maxiter),
                 as.integer(kernelType)[1],
                 as.integer(normalised)[1],
                 as.double(cut)[1],
                 as.integer(verbose),
                 clusters=integer(nrow(data)),
                 errors=as.double(rep(-1,1+length(radii)*maxiter)),
                 PACKAGE="yasomi")
    prototypes <- matrix(result$proto,ncol=ncol(prototypes),
                         dimnames=list(NULL,dimnames(data)[[2]]))
    res <- list(somgrid=somgrid,
                prototypes=prototypes,
                classif=result$cluster+1,
                errors=result$errors[result$errors>=0])
    class(res) <- c("somnum","som")
    res
}

batchsom.lowlevelcontinuous <- function(somgrid,data,prototypes,
                                        assignment,radii,weights,
                                        kernelType,normalised,cut,verbose) {
    result <- .C("batch_som_optim_continuous",
                 proto=as.double(prototypes),
                 as.integer(somgrid$size),
                 as.double(data),
                 as.integer(nrow(data)),
                 as.integer(ncol(data)),
                 as.double(weights),
                 as.integer(assignment),
                 as.double(somgrid$dist),
                 as.double(radii),
                 as.integer(length(radii)),
                 as.integer(kernelType)[1],
                 as.integer(normalised)[1],
                 as.double(cut)[1],
                 as.integer(verbose),
                 clusters=integer(nrow(data)),
                 errors=as.double(rep(-1,1+length(radii))),
                 PACKAGE="yasomi")
    prototypes <- matrix(result$proto,ncol=ncol(prototypes),
                         dimnames=list(NULL,dimnames(data)[[2]]))
    res <- list(somgrid=somgrid,
                prototypes=prototypes,
                classif=result$cluster+1,
                errors=result$errors[result$errors>=0])
    class(res) <- c("somnum","som")
    res
}


colorCode <- function(data,nbcolor) {
    onedimgrid <- somgrid(xdim=nbcolor,ydim=1,topo="rectangular")
    colorsom <- batchsom(data,onedimgrid,radii=somradii(onedimgrid,nb=20))
    colorsom$classif
}

mapToUnit <- function(som,values) {
    if(is.null(nrow(values))) {
        if(length(values)!=length(som$classif)) {
            stop("'values' is not of the same size as the data used to fit the 'som'")
        }
    } else {
        if(nrow(values)!=length(som$classif)) {
            stop("'values' is not of the same size as the data used to fit the 'som'")
        }
    }
    result <- vector("list",nrow(som$prototypes))
    if(is.null(dim(values))) {
        for(i in 1:length(result)) {
            result[[i]] <- values[som$classif==i]
        }
    } else {
        for(i in 1:length(result)) {
            result[[i]] <- values[som$classif==i,]
        }
    }
    result
}

mapFactorToUnit <- function(som,values) {
    if(!is.factor(values)) {
        stop("'values' is not a factor")
    }
    if(length(values)!=length(som$classif)) {
        stop("'values' is not of the same size as data use to fit the 'som'")
    }
    lv <- levels(values)
    result <- matrix(0,nrow=nrow(som$prototypes),ncol=length(lv),
                     dimnames=list(c(),lv))
    coded <- unclass(values)
    for(i in 1:length(som$classif)) {
        result[som$classif[i],coded[i]] <- result[som$classif[i],coded[i]] + 1
    }
    result
}


predict.somnum <- function(object,newdata,...) {
    som <- object
    newdata <- as.matrix(newdata)
    if(ncol(newdata)!=ncol(som$prototypes)) {
        stop("'newdata' and 'object$prototypes' have different dimensions")
    }
    pre <- bmu(som$prototypes,newdata)
    list(classif=pre$clusters,error=pre$error)
}

print.som <- function(x,...)
{
    cat("Parameters:\n")
    cat("       grid: ",x$somgrid$topo," grid of size ",x$somgrid$xdim,"x",
        x$somgrid$ydim," with diameter ",x$somgrid$diam,"\n",sep="")
    cat("     kernel: ",x$kernel,
        if(!x$normalised) {" not normalised"} else {" normalised"},"\n",sep="")
    cat(" assignment:",x$assignment,"\n")
    if(length(x$radii)<6) {
        cat("      radii:",x$radii,"\n")
    } else {
        cat("      radii:",length(x$radii),"values from",max(x$radii),
            "down to",min(x$radii),"\n")
    }
    cat("\n")
    invisible(x)
}

print.somnum <- function(x,...)
{
    cat("\nSelf-Organising Map\n\n")
    NextMethod()
}


summary.som <- function(object,...)
{
    cat("Quantisation error: ",object$errors[length(object$errors)],"\n\n",
        sep="")
    cat("Parameters:\n")
    cat("       grid: ",object$somgrid$topo," grid of size ",
        object$somgrid$xdim,"x",object$somgrid$ydim," with diameter ",
        object$somgrid$diam,"\n",sep="")
    cat("     kernel: ",object$kernel,
        if(!object$normalised) {" not normalised"} else {" normalised"},
        "\n",sep="")
    cat(" assignment:",object$assignment,"\n")
    cat("      radii:",object$radii,"\n")
    invisible()
}

summary.somnum <- function(object,...)
{
    cat("\nSelf-Organising Map\n\n")
    NextMethod()
}

as.dist.somnum <- function(x,FUN=NULL) {
    ## the default Euclidean distance is what we want
    dist(x$prototypes)
}

as.matrix.somnum <- function(x, ...) {
    ## the default Euclidean distance is what we want
    as.matrix(dist(x$prototypes),diag=0)
}
