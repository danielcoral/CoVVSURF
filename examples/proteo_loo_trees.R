# load("proteo.Rdata")
load("/Users/jeromesaracco/Desktop/Proteo/proteo.Rdata")

X <- proteo$x
CoVloo(X, prefix = "../../Proteo/proteo_loo_trees/proteo", mc.cores = 24)
