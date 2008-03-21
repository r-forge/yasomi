## in quality.R
KaskiLagus <- function(som,data) {
    UseMethod("KaskiLagus")
}

## in umatrix.R
protoDist <- function(som,i,j,k,l) {
    UseMethod("protoDist")
}

## in som.R (slow-versions.R) and relational.R (slow-relational.R)

sominit.pca <- function(data,somgrid,...) {
    UseMethod("sominit.pca")
}

sominit.random <- function(data,somgrid,method=c("prototypes","random","cluster"),...) {
    UseMethod("sominit.random")
}

batchsom <- function(data,somgrid,prototypes,...) {
    UseMethod("batchsom")
}
