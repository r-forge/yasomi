som.tunecontrol <- function(somgrid,init="pca",ninit=1,
                            assignment="single",
                            radii=c(2,2/3*somgrid$diam),nradii=10,
                            innernradii=30,maxiter=75,
                            annealing="power",
                            kernel="gaussian",criterion=error.quantisation) {
    if(init=="pca" && ninit>1) {
        warning("PCA based initialization is deterministic")
    }
    list(init=match.arg(init,c("pca","random"),several.ok = TRUE),
         ninit=ninit,
         assignment=match.arg(assignment,c("single","heskes"),several.ok = TRUE),
         radii=seq(from=min(radii),to=max(radii),length.out=nradii),
         innernradii=innernradii,
         maxiter=maxiter,
         annealing=match.arg(annealing,c("power","linear"),several.ok = TRUE),
         kernel=match.arg(kernel,c("gaussian","linear"),several.ok = TRUE),
         criterion=criterion)
}


som.tune <- function(somgrid,data,control=som.tunecontrol(somgrid),verbose=FALSE) {
    nbconf <- length(control$init)*control$ninit*length(control$assignment)*length(control$radii)*length(control$annealing)*length(control$kernel)
    bestSOM <- NULL
    performances <- rep(NA,nbconf)
    bestPerfSoFar <- Inf
    confIndex <- 1
    ## compute distances if they are missing
    if(is.null(somgrid$dist)) {
        somgrid$dist <- as.matrix(dist(somgrid$pts,method="Euclidean"),diag=0)
    }
    ## initialization method
    for(init in control$init) {
        ## assignment
        for(assignment in control$assignment) {
            ## kernel
            for(kernel in control$kernel) {
                if(kernel=="gaussian") {
                    minRadius <- 0.5
                } else {
                    minRadius <- 1
                }
                ## annealing
                for(annealing in control$annealing) {
                    ## initializations
                    for(i in 1:control$ninit) {
                        if(init=="pca") {
                            prototypes=somPCAInit(somgrid,data)
                        } else {
                            prototypes=somRandomInit(somgrid,data)
                        }
                        ## radii
                        for(radius in control$radii) {
                            if(verbose) {
                                print(paste("Configuration ",confIndex,"/",nbconf,sep=""))
                            }
                            radii <- switch(annealing,
                                            "linear"=radius.lin(minRadius,radius,control$innernradii),
                                            "power"=radius.exp(minRadius,radius,control$innernradii))
                            som <- batchsom(somgrid,data,prototypes,assignment,radii,
                                            maxiter=control$maxiter,kernel=kernel)
                            performances[confIndex] <- control$criterion(som,data)
                            if(performances[confIndex]<bestPerfSoFar) {
                                bestSOM <- som
                                bestPerfSoFar <- performances[confIndex]
                            }
                            confIndex <- confIndex + 1
                        }
                    }
                }
            }
        }
    }
    list(best.som=bestSOM,errors=performances)
}
