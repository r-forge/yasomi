## in quality.R
error.quantisation <- function(som,newdata,...) {
    UseMethod("error.quantisation")
}

error.kaskilagus <- function(som,newdata,...) {
    UseMethod("error.kaskilagus")
}

## in umatrix.R
protoDist <- function(som,i,j,k,l,check) {
    UseMethod("protoDist")
}

## in som.R (slow-versions.R), relational.R (slow-relational.R) and kernel.R

sominit.pca <- function(data,somgrid,...) {
    UseMethod("sominit.pca")
}

sominit.random <- function(data,somgrid,method=c("prototypes","random","cluster"),...) {
    UseMethod("sominit.random")
}

batchsom <- function(data,somgrid,init=c("pca","random"),prototypes,weights,
                     mode = c("continuous","stepwise"),
                     min.radius, max.radius, steps,
                     decrease = c("power", "linear"), max.iter,
                     kernel = c("gaussian", "linear"), normalised,
                     assignment = c("single", "heskes"),
                     cut = 1e-07,
                     verbose=FALSE,keepdata=TRUE,...) {
    UseMethod("batchsom")
}

## in kernel.R

as.kernelmatrix <- function(data,...) {
    UseMethod("as.kernelmatrix")
}


## in annealing.R

batchsom.control <- function(data,somgrid,
                             mode = c("continuous","stepwise"),
                             min.radius, max.radius, steps,
                             decrease = c("power", "linear"), max.iter,
                             kernel = c("gaussian", "linear"),
                             normalised,
                             assignment = c("single", "heskes"),
                             cut = 1e-07,...) {
    UseMethod("batchsom.control")
}

## in kmeans.R

batchkmeans <- function(data,ncenters,
                        init=c("prototypes","random","cluster"),
                        prototypes,weights,max.iter,
                        verbose=FALSE,keepdata=TRUE,...) {
    UseMethod("batchkmeans")
}
