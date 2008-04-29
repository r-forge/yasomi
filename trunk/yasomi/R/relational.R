# relational SOM

sominit.random.dist <- function(data,somgrid,
                                method=c("prototypes","random","cluster"),...) {
    method <- match.arg(method)
    dim <- nrow(data)
    nb <- somgrid$size
    if(method=="prototypes" || (method=="cluster" && nb>=dim)) {
        protos <- matrix(0,ncol=dim,nrow=nb)
        protos[cbind(1:nb,sample(1:dim,size=nb,replace=nb>dim))] <- 1
    } else if(method=="random") {
        protos <- matrix(runif(dim*nb),ncol=dim,nrow=nb)
        protos <- sweep(protos,1,rowSums(protos),"/")
    } else {
        ## nb <dim
        clusters <- cut(sample(1:dim),nb,labels=FALSE,include.lowest=TRUE)
        protos <- matrix(0,ncol=dim,nrow=nb)
        protos[cbind(clusters,1:dim)] <- 1
        protos <- sweep(protos,1,rowSums(protos),"/")
    }
    protos
}

sominit.pca.dist <- function(data, somgrid, nbsupport=3,
                             type=c("closest","random"),...) {
    type <- match.arg(type)
    ## the distance matrix PCA is implemented in cmdscale
    data.cmd <- cmdscale(data)
    sdev <- sd(data.cmd)

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
    ## this done by finding a barycentric representation for each point of the
    ## grid (this is slow but still very quick compared to cmdscale)
    if(type=="closest") {
        dbd <- dist(base,data.cmd)
        winners <- t(apply(dbd,1,order))[,1:nbsupport]
    }
    prototypes <- matrix(0,ncol=nrow(data.cmd),nrow=nrow(base))
    for(i in 1:nrow(base)) {
        target <- matrix(c(base[i,],1),ncol=1)
        if(type=="closest") {
            winner <- winners[i,]
        } else {
            winner <- sample(1:(nrow(data.cmd)),size=nbsupport)
        }
        H <- rbind(t(data.cmd[winner,]),rep(1,nbsupport))
        if(nbsupport==3) {
            weights <- solve(H,target)
        } else {
            H.svd <- svd(H)
            weights <- H.svd$v%*%diag(1/H.svd$d)%*%crossprod(H.svd$u,target)
        }
        prototypes[i,winner] <- weights
    }
    prototypes
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
    error <- sum(distances[cbind(bmu,1:length(bmu))])
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

nfPS <- function(bips,nvnormed,nclust) {
    ## th_bips_h modifies only the nf parameter: DUP=FALSE is safe here
    .C("th_bips_h",as.double(bips),as.double(nvnormed),as.integer(nclust),
       nf=double(nclust),DUP=FALSE)$nf
}

relationalsecondbmu.R <- function(prototypes,diss) {
    ## first compute the base distances
    Dalpha <- tcrossprod(diss,prototypes)
    ## then the normalisation factor
    ## suboptimal
    nf <- 0.5*diag(prototypes%*%Dalpha)
    distances <- sweep(Dalpha,2,nf,"-")
    ## very suboptimal
    ordered <- apply(distances,1,order)
    winners <- t(ordered[1:2,])
    error <- distances[cbind(1:nrow(winners),winners[,1])]
    list(winners=winners,error=error)
}

fastRelationalsom.lowlevel.R <- function(somgrid,diss,prototypes,
                                     assignment,radii,maxiter,kernel,
                                     normalised,cut,verbose) {
    oldClassif <- rep(NA,nrow(diss))
    errors <- vector("list",length(radii))
    nv <- neighborhood(somgrid,radii[1],kernel,normalised=normalised)
    ## a round of initialisation is needed
    bmus <- relationalbmu.R(prototypes,diss)
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
                bmus <- fastRelationalBMU.R(classif,somgrid$size,diss,nv)
            } else {
                stop(paste(assignment,"is not implemented for relational SOM"))
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
            if(noChange || hasLoop) {
                if(verbose) {
                    if(noChange) {
                        print(paste("radius:",radii[i],"iteration",j,
                                    "is stable, decreasing radius"))
                    } else {
                        print(paste("radius:",radii[i],"iteration",j,
                                    "oscillation detected, decreasing radius"))
                    }
                }
                if(i<length(radii)) {
                    ## preparing the loop with the next radius
                    nv <- neighborhood(somgrid,radii[i+1],kernel,
                                       normalised=normalised)
                }
                break;
            }
        }
        if(!noChange && verbose) {
            print(paste("warning: can't reach a stable configuration with radius",i))
        }
    }
    ## for latter use
    weights <- nv[,classif]
    normed <- rowSums(weights)
    prototypes <- sweep(weights,1,normed,"/")
    Dalpha <- bmus$Dalpha
    nf <- bmus$nf
    res <- list(somgrid=somgrid,prototypes=prototypes,classif=classif,
                errors=unlist(errors),Dalpha=t(Dalpha),nf=nf)
    class(res) <- c("relationalsom","som")
    res
}

batchsom.dist <- function(data,somgrid,init=c("pca","random"),prototypes,
                          assignment=c("single","heskes"),radii,nbRadii=30,
                          maxiter=75,
                          kernel=c("gaussian","linear"),normalised,
                          cut=1e-7,verbose=FALSE,keepdata=TRUE,
                          lowlevel=fastRelationalsom.lowlevel.R,...) {
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
    diss <- as.matrix(data^2,diag=0)
    if(missing(prototypes)) {
        ## initialisation based on the value of init
        init <- match.arg(init)
        args <- list(...)
        params <- c(list(data=data,somgrid=somgrid),list(...))
        prototypes <- switch(init,
                             "pca"=do.call("sominit.pca",params),
                             "random"=do.call("sominit.random",params))
    } else {
        if(ncol(prototypes)!=ncol(diss)) {
            stop("'prototypes' and 'diss' have different dimensions")
        }
    }
    ## distances?
    if(is.null(somgrid$dist)) {
        somgrid$dist <- as.matrix(dist(somgrid$pts,method="Euclidean"),diag=0)
    }
    ## compute radii
    if(missing(radii)) {
        if(kernel=="gaussian") {
            minRadius <- 0.5
        } else {
            minRadius <- 1
        }
        radii <- radius.exp(minRadius,max(minRadius,somgrid$diam/3*2),nbRadii)
    }
    pre <- lowlevel(somgrid,diss,prototypes,assignment,radii,maxiter,theKernel,
                    normalised,cut,verbose)
    pre$assignment <- assignment
    pre$kernel <- kernel
    pre$normalised <- normalised
    pre$radii <- radii
    if(keepdata) {
        pre$data  <- data
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
    themin <- min(thedist)
    if(themin < -sqrt(.Machine$double.eps)) {
        warning(paste("Non negligible negative minimal values (",themin,") removed",sep=""))
    }
    thedist[thedist<0] <- 0
    sqrt(thedist)
}

as.dist.relationalsom <- function(x,FUN=NULL) {
    as.dist(as.matrix(x))
}
