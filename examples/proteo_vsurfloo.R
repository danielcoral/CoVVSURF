ncores <- 20

load("~/Proteo/proteo.Rdata")
# load("/Users/jeromesaracco/Desktop/Proteo/proteo.Rdata")
X <- proteo$x
y <- proteo$y
table(y)
dim(X)

sta <- Sys.time()
vsurfloo <- vsurfloo(X, y, nfor = 25, vsurfs.object = "proteo_vsurf_loo", clusterType = "FORK", ncores=ncores)
vsurfloo.time <- Sys.time() - sta

proteo_vsurf_loo <- list(vsurfloo = vsurfloo, time = vsurfloo.time)

save(proteo_vsurf_loo, file="../../results/proteo_vsurf_loo.Rdata")
