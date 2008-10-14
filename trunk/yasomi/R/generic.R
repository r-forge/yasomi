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

## in som.R (slow-versions.R) and relational.R (slow-relational.R)

sominit.pca <- function(data,somgrid,...) {
    UseMethod("sominit.pca")
}

sominit.random <- function(data,somgrid,method=c("prototypes","random","cluster"),...) {
    UseMethod("sominit.random")
}

batchsom <- function(data,somgrid,init,prototypes,...) {
    UseMethod("batchsom")
}

## in kernel.R

as.kernelmatrix <- function(data,...) {
    UseMethod("as.kernelmatrix")
}


## in annealing.R

batchsom.control <- function(data,somgrid,
                             mode = c("stepwise","continuous"),
                             min.radius, max.radius, steps,
                             decrease = c("power", "linear"), max.iter,
                             kernel = c("gaussian", "linear"),
                             normalised,
                             assignment = c("single", "heskes"),
                             cut = 1e-07,...) {
    UseMethod("batchsom.control")
}
