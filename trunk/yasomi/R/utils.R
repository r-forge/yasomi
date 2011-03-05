split.on.prefix <- function(prefix,args) {
    locate.pattern <- paste(prefix,"\\.",sep="")
    with.prefix <- grep(locate.pattern,names(args))
    if(length(with.prefix)>0) {
        new.names <-
            sub(paste(locate.pattern,"(.*)",sep=""),"\\1",names(args)[with.prefix])
        wp <- list()
        wp[new.names] <- args[with.prefix]
        list(with.prefix=wp,others=args[-with.prefix])
    } else {
        list(with.prefix=NULL,others=args)
    }
}

flex.scale <- function(X,targets,global=FALSE) {
    if(global) {
        X.range <- range(X,na.rm=TRUE)
        (X-X.range[1])/diff(X.range)*diff(targets)+targets[1]
    } else {
        apply(X,2,function(x) {x.range <- range(x,na.rm=TRUE); (x-x.range[1])/diff(x.range)*diff(targets)+targets[1]})
    }
}

flex.scale.list <- function(l,targets,global=FALSE) {
    if(global) {
        l.range <- c(min(unlist(sapply(l,function(x) if(nrow(x)>0) min(x,na.rm=TRUE)))),
                     max(unlist(sapply(l,function(x) if(nrow(x)>0) max(x,na.rm=TRUE)))))
        l.scale <- diff(targets)/diff(l.range)
        lapply(l,function(x) (x-l.range[1])*l.scale + targets[1])
    } else {
        l.ranges <- lapply(l,function(x) if (nrow(x)>0) apply(x,2,range,na.rm=TRUE))
        l.mins <- apply(sapply(l.ranges,function(x) if(is.null(x)) rep(Inf,ncol(l[[1]])) else x[1,]),1,min)
        l.maxs <- apply(sapply(l.ranges,function(x) if(is.null(x)) rep(-Inf,ncol(l[[1]])) else x[2,]),1,max)
        l.scales <- diff(targets)/(l.maxs - l.mins)
        lapply(l,function(x) if(nrow(x)>0) sweep(sweep(x,2,l.mins,"-"),2,l.scales,"*")+targets[1] )
    }
}

to.matrix <- function(l) {
    if(is.matrix(l[[1]]) || is.data.frame(l[[1]])) {
        nbrows <- sum(sapply(l,function(x) nrow(x)))
        nbcols <- ncol(l[[1]])
        res <- matrix(NA,ncol=nbcols,nrow=nbrows,dimnames=list(as.character(1:nbrows),colnames(l[[1]])))
        classif <- rep(NA,nbrows)
        pos <- 1
        for(i in seq_along(l)) {
            if(!is.null(l[[i]]) && nrow(l[[i]])>0) {
                res[pos:(pos+nrow(l[[i]])-1),] <- as.matrix(l[[i]])
                classif[pos:(pos+nrow(l[[i]])-1)] <- i
                pos <- pos+nrow(l[[i]])
            }
        }
        list(m=res,classif=classif)
    } else if(is.factor(l[[1]])) {
        res <- t(sapply(l,function(x) tabulate(x,nbins=length(levels(l[[1]])))))
        colnames(res) <- levels(l[[1]])
        classif <- seq_along(l)
        list(m=res,classif=classif,levels=levels(l[[1]]))
    } else {
        stop(paste("cannot handle data type:",class(l[[1]])))
    }
}

default.args <- function(args,replace=FALSE,...) {
    call <- match.call(expand.dots=FALSE)
    for(arg in names(call$...)) {
        if(is.null(args[[arg]]) | replace) {
            args[[arg]] <- call$...[[arg]]
        } 
    }
    args
}

