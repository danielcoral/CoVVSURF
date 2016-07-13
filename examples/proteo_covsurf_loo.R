ncores <- 25

load("~/Proteo/proteo.Rdata")
# load("/Users/jeromesaracco/Desktop/Proteo/proteo.Rdata")
X <- proteo$x
y <- proteo$y
table(y)
dim(X)

kval <- 50:200

sta <- Sys.time()
csloo <- covsurfloo(X, y, kval = kval, avail.tree = TRUE,
                    prefix = "~/Proteo/proteo_loo_trees/proteo",
#                     prefix = "/Users/jeromesaracco/Desktop/Proteo/proteo_loo_trees/proteo",
                    clusterType = "FORK", ncores=ncores)
csloo.time <- Sys.time() - sta

proteo_covsurf_loo <- list(csloo = csloo, csloo.time = csloo.time)

save(proteo_covsurf_loo, file="proteo_covsurf_loo.Rdata")
