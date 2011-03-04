rbf.hardy <- function(r,alpha) {
    sqrt(r^2+alpha^2)
}

rbf.gaussian <- function(r,sigma) {
    exp(-((r/sigma)^2)/2)
}

V.zero <- function(d) {
    1 - d
}

V.one <- function(d) {
    1 + d^2 * ( 2 * d - 3)
}

V.two <- function(d) {
    1 + d^3 * ( -10 + d * ( 15 - 6 * d))
}

inbox <- function(box,pts) {
    (pts[,1]>=box[1])&(pts[,2]>=box[2])&(pts[,1]<=box[3])&(pts[,2]<=box[4])
}

dbox <- function(box,pts) {
    L.x <- pts[,1]-box[1]
    L.y <- pts[,2]-box[2]
    U.x <- box[3]-pts[,1]
    U.y <- box[4]-pts[,2]
    1-16*(L.x*U.x)/((box[3]-box[1])^2)*(L.y*U.y)/((box[4]-box[2])^2)
}

boxes <- function(pts,depth,delta) {
    lims <- apply(pts,2,range)
    basesize <- 2^depth
    bounds.x <- seq(lims[1,1],lims[2,1],length=basesize+1)
    bounds.y <- seq(lims[1,2],lims[2,2],length=basesize+1)
    as.matrix(cbind(expand.grid(bounds.x[1:basesize],
                                bounds.y[1:basesize])-delta,
                    expand.grid(bounds.x[2:(basesize+1)],
                                bounds.y[2:(basesize+1)])+delta))
}

gaussianRbfInterp <- function(pts,z,sigma=0.5) {
    ptdist <- as.matrix(dist(pts,method="Euclidean"),diag=0)
    A <- rbf.gaussian(ptdist,sigma)
    coeffs <- try(solve(A,as.vector(z)))
    if(inherits(coeffs,"try-error")) {
        print(A)
        stop(coeffs)
    }
    structure(list(coeffs=coeffs,points=pts,sigma=sigma),
              class="gaussianRbfInterp")
}


predict.gaussianRbfInterp <- function(object,newdata,...) {
    distances <- as.matrix(dist(newdata,object$points,method="Euclidean"))
    class(distances) <- "matrix"
    A <- rbf.gaussian(distances,object$sigma)
    A %*% object$coeffs
}

terrainInterp <- function(pts,z,depth,delta,sigma=0.5) {
    z <- as.vector(z)
    nbBoxes <- (2^depth)^2
    ## build boxes
    myboxes <- boxes(pts,depth,delta)
    ## assign points to boxes
    membership <- matrix(NA,ncol=nbBoxes,nrow=nrow(pts))
    for(i in 1:nbBoxes) {
        membership[,i] <- inbox(myboxes[i,],pts)
    }
    ## compute a rbf model for each box
    models <- vector("list",nbBoxes)
    for(i in 1:nbBoxes) {
        models[[i]] <- gaussianRbfInterp(pts[membership[,i],],
                                         z[membership[,i]],sigma=sigma)
    }
    structure(list(models=models,boxes=myboxes),class="terrainInterp")
}

predict.terrainInterp <- function(object,newdata,...) {
    result <- rep(0,nrow(newdata))
    weights <- rep(0,nrow(newdata))
    nbBoxes <- nrow(object$boxes)
    for(i in 1:nbBoxes) {
        ## get points falling in this box
        members <- inbox(object$boxes[i,],newdata)
        ## predictions for those points
        predicted <- predict(object$models[[i]],newdata[members,])
        ## weights
        delta <- V.two(dbox(object$boxes[i,],newdata[members,]))
        result[members] <- result[members] + delta*predicted
        weights[members] <- weights[members] + delta
    }
    result/weights
}
