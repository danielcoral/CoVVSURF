ncores <- 25

load("~/Proteo/proteo.Rdata")
# load("/Users/jeromesaracco/Desktop/Proteo/proteo.Rdata")
X <- proteo$x
y <- proteo$y
table(y)
dim(X)

load("../../Proteo/proteo_tree.Rdata")
tree <- proteo_tree
plot(tree)
plot(tree,type="index")

sta <- Sys.time()
kval <- 2:2000
cs <- covsurf(X, y, kval=kval, tree=tree, clusterType = "FORK", ncores=ncores)
cs.time <- Sys.time() - sta

plot(kval, cs$oob[, 1], type="l", xlab="Partition cardinal", ylab="OOB error rate")
cs$vsel
cs$csel

# VSURF comparison
vsurf <- VSURF.parallel(X, y, clusterType = "FORK", ncores=ncores)

proteo_covsurf_and_vsurf <- list(covsurf = cs, cs.time = cs.time, kval = kval, vsurf = vsurf)

save(proteo_covsurf_and_vsurf , file="proteo_covsurf_and_vsurf.Rdata")
