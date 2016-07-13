load("~/Proteo/proteo.Rdata")
# load("/Users/jeromesaracco/Desktop/Proteo/proteo.Rdata")
X <- proteo$x
y <- proteo$y
table(y)
dim(X)

proteo_covsurf_loo <- covsurfloo(X, y, covsurfs.object = "proteo_covsurf_loo",
                    path = "~/Proteo/", mcores = 3)

save(proteo_covsurf_loo[1:2], file="../../results/proteo_covsurf_loo.Rdata")
