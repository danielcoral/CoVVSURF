ncores <- 3
# ncores <- 20

load("~/Proteo/proteo.Rdata")
# load("/Users/jeromesaracco/Desktop/Proteo/proteo.Rdata")
X <- proteo$x
y <- proteo$y
table(y)
dim(X)

proteo_vsurf <- VSURF.parallel(X, y, clusterType = "FORK", ncores = ncores)

save(proteo_vsurf, file="proteo_vsurf.Rdata")
