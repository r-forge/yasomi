opar <- par(ask = dev.interactive(orNone = TRUE))

data <- scale(iris[1:4])

# building the grid

sg <- somgrid(xdim=15,ydim=15,topo="hex")

# training the SOM

somT <- som.tune(data,sg,som.tunecontrol(sg,radii=c(2,sg$diam),nradii=20,criterion=error.kaskilagus),verbose=T)

# and displaying the results

plot(somT$errors,type="h")

som <- somT$best.som

plot(som$error,type="l")

code <- colorCode(data,50)
pairs(data,bg=rainbow(50)[code],pch=21)

plot(som)

plot(som,iris[[5]])

spar <- par(mfrow=c(2,2))
for(i in 1:ncol(data)) {
  componentPlane(som,i)
}
par(spar)

hitMap(som,col="blue")

umatrix(som)

pdist <- prototype.distances(som)

pgrid <- distance.grid(pdist)

filled.contour(pgrid)

par(opar)
