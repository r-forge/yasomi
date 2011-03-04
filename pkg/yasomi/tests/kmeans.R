library(yasomi)
## first test that the kernel kmeans gives identical results as the standard kmeans
datadim <- ceiling(4+runif(1,min=1,max=16))
datasize <- datadim*100

X <- matrix(nrow=datasize,ncol=datadim,runif(datasize*datadim))
KX <- as.kernelmatrix(tcrossprod(X))

ncenters <- ceiling(sqrt(datasize)+runif(1,min=-2,max=4))

## random initialization by random prototypes
init.points <- sample(datasize,ncenters)

kprototypes.init <- matrix(0,nrow=ncenters,ncol=datasize)
kprototypes.init[cbind(1:ncenters,init.points)] <- 1
prototypes.init <- X[init.points,]

## let's fit the kmeans
km <- batchkmeans(X,ncenters,prototypes=prototypes.init,verbose=TRUE)

## and the kernel version
kkm <- batchkmeans(KX,ncenters,prototypes=kprototypes.init,verbose=TRUE)
kprototypes.final <- kkm$prototypes%*%X

## now compare the results
stopifnot(all.equal(error.quantisation(km),error.quantisation(kkm)))
stopifnot(all.equal(km$errors,kkm$errors))
stopifnot(all.equal(km$classif,kkm$classif))
stopifnot(all.equal(km$prototypes,kprototypes.final,
                    check.attributes = FALSE, check.names = FALSE))


## then test weights support

## new starting point
init.points <- sample(datasize,ncenters)
prototypes.init <- X[init.points,]

## random integer weights
weights <- sample(5,size=nrow(X),replace=T,prob=exp(-seq(0,2,length=5)))

## with weight support
km <- batchkmeans(X,ncenters,prototypes=prototypes.init,weights=weights,verbose=TRUE)

## weighting via replication
X.rep <- X[rep(1:nrow(X),times=weights),]
km.rep <- batchkmeans(X.rep,ncenters,prototypes=prototypes.init,verbose=TRUE)

stopifnot(all.equal(km.rep$prototypes,km$prototypes))
stopifnot(all.equal(km.rep$errors,km$errors))
stopifnot(all.equal(error.quantisation(km),error.quantisation(km.rep)))
stopifnot(all.equal(km$classif[rep(1:nrow(X),times=weights)],km.rep$classif))

## same test for kernel k-means
kprototypes.init <- matrix(0,nrow=ncenters,ncol=datasize)
kprototypes.init[cbind(1:ncenters,init.points)] <- 1

kkm <- batchkmeans(KX,ncenters,prototypes=kprototypes.init,weights=weights,verbose=TRUE)
kprototypes.final <- kkm$prototypes%*%X

stopifnot(all.equal(kprototypes.final,km$prototypes))
stopifnot(all.equal(kkm$errors,km$errors))
stopifnot(all.equal(error.quantisation(km),error.quantisation(kkm)))
stopifnot(all.equal(km$classif,kkm$classif))
