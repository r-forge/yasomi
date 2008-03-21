opar <- par(ask = dev.interactive(orNone = TRUE))

# preparing the data
set.seed(42)
X <- cbind(rnorm(200,mean=2,sd=0.35),rnorm(200,mean=-1,sd=0.35))
Y <- cbind(runif(200,min=-1.5,max=-0.75),runif(200,min=0,max=0.5))
Z <- cbind(rnorm(200,sd=0.15),rnorm(200,sd=0.5))
M <- matrix(c(sin(pi/4),cos(pi/4),-cos(pi/4),sin(pi/4)),ncol=2)
data <- rbind(X,Y,Z%*%M+c(rep(0.25,200),rep(-0.5,200)))
data <- scale(data)

# building the grid
sg <- somgrid(xdim=15,ydim=15,topo="rect")

# training the SOM

somT <- som.tune(data,sg,som.tunecontrol(sg,radii=c(2,sg$diam),nradii=20,criterion=error.kaskilagus),verbose=T)

# and displaying the results

plot(somT)

som <- somT$best.som

plot(som$error,type="l")

plot(data,pch="+",col="red")
lines(grid2lines(som),type="b",pch=20)

code <- colorCode(data,50)
plot(data,col=rainbow(50)[code],pch=20)

componentPlane(som,1)
componentPlane(som,2)

hitMap(som,col="blue")

umatrix(som)

par(opar)
