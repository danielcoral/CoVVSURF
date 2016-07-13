oob <- matrix(NA, nrow = length(kval), ncol = 2)

sta <- Sys.time()

for (k in kval) {
  oob[which(kval == k), ] <- rfptree(tree, y, k, mtry=k/3, ncores = ncores)
  print(paste("k =", k, "done", sep=" "))
}

comput.time <- Sys.time() - sta

vslasso <- c("p_10405", "p_11709", "p_11778", "p_12040", "p_16331", "p_18172",
                 "p_20402", "p_32826", "p_4722", "p_9981")

lapply(X = vslasso, FUN = function(i) max(abs(cor(subset(x, select = i), z))))

max(abs(round(cor(subset(x, select = vslasso), subset(x, select = vsvsurf)), 2)))
