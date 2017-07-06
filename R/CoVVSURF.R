#' @title Builds a dendrogram with ClustOfVar
#' @name CoV
#' @param X A dataset
#' @importFrom PCAmixdata splitmix
#' @importFrom ClustOfVar hclustvar
#' @export
CoV <- function(X) {
  obj <- PCAmixdata::splitmix(X)
  tree <- ClustOfVar::hclustvar(obj$X.quanti, obj$X.quali)
}

#' @title Performs splitmix and PCAmix
#' @name pcamix
#' @param X A dataset
#' @param ndim Number of dimensions in PCAmix
#'
#' @importFrom PCAmixdata PCAmix
#' @export
pcamix <- function(X, ndim = 5) {
  obj <- PCAmixdata::splitmix(X)
  res <- PCAmixdata::PCAmix(obj$X.quanti, obj$X.quali, ndim = ndim,
                            rename.level = TRUE) 
}

#' @title Combining ClustOfVar and VSURF
#' @name covsurf
#' @description This function selects groups of informative input variables 
#' to predict a response variable.
#' @param X dataframe of input variables
#' @param y vector of responses
#' @param kval vector of number of classes to try
#' @param tree optional tree given by hclustvar
#' @param nse number of standard-deviation to add to select minimum of OOB rate
#' @param ncores number of cores to use for parallel computing
#' @param ... passed to VSURF
#' @return 
#' \item{kopt}{the optimal number of groups of variables}
#' \item{ptree}{the partition in \code{kopt} clusters of the dendrogram of \code{CoV}.}
#' \item{vsurf_ptree}{\code{VSURF} applied to the \code{kopt} synthetic 
#' variables.}
#' \item{vsel}{synthetic variables selected by \code{VSURF}.}
#' \item{csel}{groups of variables selected by \code{VSURF}.}
#' \item{rfsel}{RF applied to  the selected synthetic variables}
#'\item{rfclust}{RF applied to  all the synthetic variables.}
#'\item{oob}{a matrix with mean OOB error (first column) and
#' OOB standard deviation (second column).}
#' @examples 
#' data(don60)
#' kval <- c(2:15, seq(from = 20, to = ncol(X), by = 10))
#' don60covs <- covsurf(X, y, kval)
#' plot(don60covs)
#' @importFrom ClustOfVar cutreevar
#' @importFrom VSURF VSURF
#' @importFrom randomForest randomForest
#' @export
#' @references 
#' Combining clustering of variables and feature selection using random forests: the CoV/VSURF 
#' procedure, Marie Chavent, Robin Genuer, Jerome Saracco, hal-01345840
covsurf <- function(X, y, kval = 2:ncol(X), tree = NULL, nse=1, ncores = 1,...) {

  if (is.null(tree)) {
    tree <- CoV(X)
  }
  
  oob <- matrix(NA, nrow = length(kval), ncol = 2)
    
  for (k in kval) {
    oob[which(kval == k), ] <- rfptree(tree, y, k,ncores=ncores,...)
    print(paste("k =", k, "done", sep=" "))
  }
  
  indkstar <- which.min(oob[, 1])
  indkopt <- min(which(oob[, 1] <= oob[indkstar, 1] + nse*oob[indkstar, 2]))
  kopt <- kval[indkopt]
  
  ptree <- ClustOfVar::cutreevar(tree, kopt)
  xtree <- ptree$scores
  vsurf_ptree <- VSURF(xtree, y, ...)
  vsel <- vsurf_ptree$varselect.interp
  csel <- ptree$var[vsel]
  rfsel <- randomForest::randomForest(xtree[, vsel, drop = FALSE], y)
  rfclust <- randomForest::randomForest(xtree, y)
  
  output <- list(oob=oob,
                 kval=kval,
                 kopt=kopt,
                 ptree=ptree,
                 vsurf_ptree=vsurf_ptree,
                 vsel=vsel,
                 csel=csel,
                 rfsel=rfsel,
                 rfclust=rfclust,
                 y=y)
  class(output) <- "covsurf"
  output
}

#' @title Make predictions for new data with a covsurf object
#' @name predict.covsurf
#' @param object an object of class "covsurf"
#' @param newdata data to predict
#' @param select if TRUE (default): predict only with selected variables,
#' else predict with all synthetic variables
#' @param nfor number of forests to build to obtain a collection of predictors
#' @param mcores number of cores used to parallel computations
#' @importFrom PCAmixdata predict.PCAmix
#' @importFrom parallel mclapply
#' @importFrom stats predict
#' @method predict covsurf
#' @export
predict.covsurf <- function(object, newdata, select = TRUE, nfor = 1,
                            mcores = 1) {
  scores <- stats::predict(object$ptree, PCAmixdata::splitmix(newdata)$X.quanti,
                    PCAmixdata::splitmix(newdata)$X.quali)
  if (nfor == 1) {
    if (select) {
      pred <- stats::predict(object$rfsel, scores)
    }
    
    else {
      pred <- stats::predict(object$rfclust, scores)
    }
  }
  
  else {
    xtree <- object$ptree$scores
    y <- object$y
    
    one_pred_sel <- function(i) {
      stats::predict(randomForest::randomForest(xtree[, vsel, drop = FALSE], y),
              scores)
    }
    
    one_pred <- function(i) {
      stats::predict(randomForest::randomForest(xtree, y), scores)
    }
    
    if (select) {
      vsel <- object$vsurf_ptree$varselect.interp
      predictions <- simplify2array(parallel::mclapply(1:nfor, one_pred_sel,
                                             mc.cores = mcores))
    }
    
    else {
      predictions <- simplify2array(parallel::mclapply(1:nfor, one_pred,
                                             mc.cores = mcores))
    }
  }
}


#' @title VSURF on variable space partition
#' @name vsurfptree
#' @param tree result of hclustvar
#' @param k number of classes
#' @param y output variable
#' @param ... arguments to be passed on VSURF
#' @examples 
#' data(don60)
#' obj <- PCAmixdata::splitmix(X)
#' tree <- ClustOfVar::hclustvar(obj$X.quanti,obj$X.quali)
#' don60ptree2 <- vsurfptree(tree, y, k = 2)
#' @export
vsurfptree <- function(tree, y, k, ...) {
  ptree <- ClustOfVar::cutreevar(tree, k)
  xtree <- ptree$scores
  vsurf_ptree <- VSURF::VSURF(xtree, y, ...)
}

#' @title Random Forests on variable space partition
#' @name rfptree
#' @param tree result of hclustvar
#' @param k number of classes
#' @param y output variable
#' @param nfor number of random forests to build
#' @param ncores number of cores to use for parallel computing
#' @param ... arguments to be passed on VSURF
#' @importFrom stats sd
#' @importFrom utils tail
#' @export
rfptree <- function(tree, y, k, nfor = 25, ncores = 1, ...) {
  ptree <- ClustOfVar::cutreevar(tree, k)
  xtree <- ptree$scores
  
  rf <- function(i, ...) {
    if (is.factor(y))
      out <- utils::tail(randomForest::randomForest(xtree, y,
                             mtry=max(floor(k/3),1), ...)$err.rate[, 1], n=1)
    else
      out <- utils::tail(randomForest::randomForest(xtree, y,
                                                    mtry=max(floor(k/3),1), ...)$mse, n=1)
  }
  
  rfs <- parallel::mclapply(1:nfor, rf, ... , mc.cores = ncores)
  rfs <- unlist(rfs)
  
  output <- c(mean(rfs), stats::sd(rfs))
}

#' @title Construction of trees during leave-one-out procedure
#' @name CoVloo
#' @param X dataframe of input variables
#' @param prefix character string
#' @param ... passed to mclapply
#' @param mc.cores number of cores used for parallel computations
#' @export
CoVloo <- function(X, prefix = NULL, mc.cores = 1, ...) {
  
  CoVsave <- function(i) {
    tree <- CoV(X[-i, ])
    save(tree, file = paste(prefix, "treeloo", i, ".Rdata", sep = ""))
  }
  
  parallel::mclapply(1:nrow(X), CoVsave, ..., mc.cores = mc.cores)
}

#' @title ClustOfVar combine with VSURF during leave-one-out procedure
#' @name covsurfloo
#' @param X dataframe of input variables
#' @param y vector of responses
#' @param kval vector of number of classes to try
#' @param prefix character string
#' @param avail.tree logical value indicating if trees are already available
#' @param covsurfs.object character string
#' @param path path to covsurfs ojects
#' @param nfor number of forests to compute for averaged prediction
#' @param ... passed to covsurf
#' @export
covsurfloo <- function(X, y, kval = 2:ncol(X), prefix = NULL,
                       avail.tree = FALSE, covsurfs.object = NULL,
                       path = NULL, nfor = 100, ...) {
  
  if (!is.null(covsurfs.object)) {
    load(paste(path, covsurfs.object, ".Rdata", sep = ""))
    eval(parse(text = paste("covsurfs <- ", covsurfs.object,
                            "$csloo$covsurfs", sep = "")))
  }
  
  else {
    covsurfs <- vector("list", nrow(X))
  }
  
  pred_covsurf <- data.frame(matrix(NA, nrow = nrow(X), ncol = nfor))
  pred_covrf <- data.frame(matrix(NA, nrow = nrow(X), ncol = nfor))
  
  for (i in 1:nrow(X)) {
    Xl <- X[-i, ]
    yl <- y[-i]
    
    if (is.null(covsurfs.object)) {
      tree <- NULL
      if (avail.tree) {
        load( paste(prefix, "treeloo", i, ".Rdata", sep = "") )
      }
      cs <- covsurf(Xl, yl, kval = kval, tree = tree, ...)
      covsurfs[[i]] <- cs
    }
    
    else {
      cs <- covsurfs[[i]]
      cs$y <- yl
    }
    
    pred_covsurf[i, ] <- stats::predict(cs, X, nfor = nfor, ...)[i, ]
    pred_covrf[i, ] <- stats::predict(cs, X, select = FALSE, nfor = nfor, ...)[i, ]
    
    print(paste("step", i, "of", nrow(X), "done", sep=" "))
  }
  
  error_covsurf <- sapply(pred_covsurf, function(x) mean(x != y))
  error_covrf <- sapply(pred_covrf, function(x) mean(x != y))
  
  errors <- data.frame(covsurf = error_covsurf, covrf = error_covrf)
  preds <- data.frame(covsurf = pred_covsurf[,1], covrf = pred_covrf[,1])
  
  output <- list(errors=errors, preds=preds, covsurfs = covsurfs)
}

#' @title VSURF during leave-one-out procedure
#' @name vsurfloo
#' @param X dataframe of input variables
#' @param y vector of responses
#' @param vsurfs.object optional VSURF objects
#' @param path path to VSURF objects
#' @param nfor number of forests for averaged predictions
#' @param ... further parameters to be passed on VSURF function
#' @export
vsurfloo <- function(X, y, vsurfs.object = NULL, path = NULL, nfor = 100, ...) {
  
  pred_rf <- as.data.frame(matrix(NA, nrow = nrow(X), ncol = nfor))
  pred_rfinterp <- as.data.frame(matrix(NA, nrow = nrow(X), ncol = nfor))
  pred_rfpred <- as.data.frame(matrix(NA, nrow = nrow(X), ncol = nfor))
  
  if (!is.null(vsurfs.object)) {
    load(paste(path, vsurfs.object, ".Rdata", sep = ""))
    eval(parse(text = paste("vsurfs <- ", vsurfs.object,
                            "$vsurfs", sep = "")))
  }
  
  else {
    vsurfs <- vector("list", nrow(X))
  }  
  
  for (i in 1:nrow(X)) {
    Xl <- X[-i, ]
    yl <- y[-i]
    
    if (!is.null(vsurfs.object)){
      vsurf <- vsurfs[[i]]
    }
    
    else {
      vsurf <- VSURF::VSURF(Xl, yl, ...)
      vsurfs[[i]] <- vsurf
    }
    
    vsinterp <- vsurf$varselect.interp
    vspred <- vsurf$varselect.pred
       
    pred_rfinterp[i, ] <- replicate(nfor, stats::predict(
      randomForest(Xl[, vsinterp, drop = FALSE], yl),
      X[i, , drop=FALSE]) )
    pred_rfpred[i, ] <- replicate(nfor, stats::predict(
      randomForest(Xl[, vspred, drop = FALSE], yl),
      X[i, , drop=FALSE]) )
    pred_rf[i, ] <- replicate(nfor, stats::predict(randomForest(Xl, yl),
      X[i, , drop=FALSE]))
    
    print(paste("step", i, "of", nrow(X), "done", sep=" "))
  }
  
  error_interp <- colMeans(pred_rfinterp != y)
  error_pred <- colMeans(pred_rfpred != y)
  error_rf <- colMeans(pred_rf != y)
  
  errors <- data.frame(interp = error_interp, pred = error_pred, rf = error_rf)
  preds <- data.frame(interp = pred_rfinterp[, 1], pred = pred_rfpred[, 1],
                      rf = pred_rf[, 1])
  
  output <- list(errors = errors, preds=preds, vsurfs = vsurfs)
}


#' @title Plot the results from a covsurf object
#' @name plot.covsurf
#' @param x an object of class "covsurf"
#' @param choice choice of the plot to perform: defaut is \code{"covsurf"}, other choices
#' are "vsurf", ... 
#' @param ... to be passed on plot function
#' @importFrom graphics plot
#' @method plot covsurf
#' @export
plot.covsurf <- function(x, choice = "covsurf",...) {
  cs <- x
  if (is.factor(cs$y))
    ylab="OOB error rate"
  else
    ylab="OOB mean square error"
  if (choice == "covsurf") {
    kval <- cs$kval
    graphics::plot(kval, cs$oob[, 1], type="l", xlab="Partition cardinal",
         ylab=ylab, ...)
  }
  
  if (choice == "vsurf") {
    graphics::plot(cs$vsurf_ptree)
  }
}

#' @title Perform covsurf on simulated datasets
#' @name simus_covsurf
#' @param data_name name of the simulated dataset (the dataset will be loaded by 
#' \code{data(data_name)})
#' @param data_test_name name of the test dataset (the dataset will be loaded by 
#' \code{data(data_test_name)})
#' @param  nfor number of forests for averaged predictions
#' @param Rdata.object optional Rdata name
#' @param path path to Rdata object
#' @param mcores number of cores used in parallel computations
#' @param ... passed to VSURF
#' @export
simus_covsurf <- function(data_name, data_test_name = "dontest", nfor = 100, 
                          Rdata.object = NULL, path = NULL, mcores = 1, ...) {
  
  eval(parse(text = paste("data(", data_name, ")", sep="")))
  table(y)
  dim(X)
  
  # covsurf
  if (!is.null(Rdata.object)) {
    load(paste(path, Rdata.object, ".Rdata", sep = ""))
    eval(parse(text = paste("cs <- ", Rdata.object,
                            "$covsurf", sep = "")))
    eval(parse(text = paste("vsurf <- ", Rdata.object,
                            "$vsurf", sep = "")))
    eval(parse(text = paste("cs.time <- ", Rdata.object,
                            "$cs.time", sep = "")))
  }
  
  else {
    sta <- Sys.time()
    cs <- covsurf(X, y, ...)
    cs.time <- Sys.time() - sta
    vsurf <- VSURF::VSURF(X, y, ...)
  }
  
  rf <- randomForest::randomForest(X, y, ntree = 2000, mtry = ncol(X)/3)
  
  # Test error rates
  eval(parse(text = paste("data(", data_test_name, ")", sep="")))
  
  pred_covsurf <- stats::predict(cs, Xt, nfor = nfor, ...)
  pred_covrf <- stats::predict(cs, Xt, nfor = nfor, select = FALSE, ...)
  pred_rf <- stats::predict(rf, Xt)
  
  vsinterp <- vsurf$varselect.interp
  vspred <- vsurf$varselect.pred
  vsurfi <- randomForest::randomForest(X[, vsinterp, drop=FALSE], y)
  vsurfp <- randomForest::randomForest(X[, vspred, drop=FALSE], y)
  pred_vsurfi <- stats::predict(vsurfi, Xt)
  pred_vsurfp <- stats::predict(vsurfp, Xt)

  error_covsurf <- colMeans(pred_covsurf != yt)
  error_covrf <- colMeans(pred_covrf != yt)
  
  error_one_rf <- function(i) {
    err <- utils::tail(randomForest::randomForest(X[, vs, drop=FALSE], y,
                Xt[, vs, drop=FALSE], yt)$test$err.rate[, 1], 1)
  }
  
  vs <- 1:ncol(X)
  error_vsurfi <- unlist(parallel::mclapply(1:nfor, error_one_rf,
                                            mc.cores = mcores))
  vs <- vsinterp
  error_vsurfp <- unlist(parallel::mclapply(1:nfor, error_one_rf,
                                            mc.cores = mcores))
  vs <- vspred
  error_rf <- unlist(parallel::mclapply(1:nfor, error_one_rf,
                                        mc.cores = mcores))
    
  predictions <- data.frame(covsurf = pred_covsurf[, 1],
                            covrf = pred_covrf[, 1],
                            vsurfi = pred_vsurfi, vsurfp = pred_vsurfp,
                            rf = pred_rf)
  errors <- data.frame(covsurf = error_covsurf, covrf = error_covrf,
                       vsurfi = error_vsurfi, vsurfp = error_vsurfp,
                       rf = error_rf)
  
  # save
  simus_covsurf <- list(data_name = data_name, errors = errors, covsurf = cs,
                        covsurf.time = cs.time, rf = rf, vsurf = vsurf,
                        vsurfi = vsurfi, vsurfp = vsurfp,
                        predictions = predictions)
  eval(parse(text = paste(data_name, "_covsurf <- simus_covsurf", sep="")))
  eval(parse(text = paste("save(", data_name, "_covsurf, ",
    "file = \"../../results/", data_name, "_covsurf.Rdata\")", sep="")))
  eval(parse(text = paste("out <- ", data_name, "_covsurf", sep="")))
}


#' Helps to repeat experiments for simulated data
#'
#' @param n Number of observations in simulated data
#' @param nrep Number of repetitions of the simulations
#' @param nfor Number of forests prediction averaged
#' @param ... Arguments passed to \code{\link{simus_covsurf}} function
#'
#' @export
repeated_simus <- function(n = 600, nrep = 50, nfor = 50, ...) {
  
  one_simu <- function(i) {
    eval(parse(text = paste("simu(n = ", n, ", name = \"don", n, "_", i, "\")",
                            sep="")))
    eval(parse(text = paste("onesim <- simus_covsurf(data_name = ", "\"don", n,
                            "_", i, "\", nfor = nfor, ...)$errors", sep ="")))
    print(paste("repetition", i, "done", sep=" "))
    out <- onesim
  }
  
  res <- lapply(1:nrep, one_simu)
  
  errors <- NULL
  for (i in 1:nrep) {
    errors <- rbind(errors, res[[i]])
  }
  
  res$errors <- errors
  
  eval(parse(text = paste("don", n, "_repeated_simus <- res", sep="")))
  eval(parse(text = paste("save(don", n, "_repeated_simus, ",
    "file = \"../../results/don", n, "_repeated_simus.Rdata\")", sep="")))
  eval(parse(text = paste("out <- don", n, "_repeated_simus", sep="")))
}
  
  
  
  
  
  
  