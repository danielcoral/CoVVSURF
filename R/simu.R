#' @title Simulations of data
#' @name simu_reg
#' @param n number of observations
#' @param rho correlation between variables within a group
#' @param sigma standard deviation of the variables in the group of noise
#' @param seed optional seed to pass to set.seed
#' @return the matrix X of explonatory variables and y the numeric variable to predict
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats median
#' @export
#' @references 
#' Combining clustering of variables and feature selection using random forests: the CoV/VSURF 
#' procedure, Marie Chavent, Robin Genuer, Jerome Saracco, hal-01345840
simu_reg <- function(n=600,rho=0.9,sigma=1,seed=NULL)
{
  
  if (!is.null(seed)) set.seed(seed)
  #dans chaque tableau, 3 grpes de variables correlees de taille p1,p2,p3
  p1 <- 3
  Sigma1 <- matrix(rho,p1,p1)
  diag(Sigma1) <-rep(1,p1)
  
  p2 <- 15
  Sigma2<- matrix(rho,p2,p2)
  diag(Sigma2) <-rep(1,p2)
  
  p3 <- 12
  Sigma3<- matrix(rho,p3,p3)
  diag(Sigma3) <-rep(1,p3)
  
  Sigma <-diag(p1+p2+p3)
  Sigma[1:p1,1:p1] <- Sigma1
  Sigma[(p1+1):(p1+p2),(p1+1):(p1+p2)] <- Sigma2
  Sigma[(p1+p2+1):(p1+p2+p3),(p1+p2+1):(p1+p2+p3)] <- Sigma3
  
  #creation de X1,X2 et X3 qui formeront le groupe num, qual, mixte
  mu <- rep(0,p1+p2+p3)
  X1 <- mvtnorm::rmvnorm(n,mu,Sigma)
  X2 <- mvtnorm::rmvnorm(n,mu,Sigma)
  X3<- mvtnorm::rmvnorm(n,mu,Sigma)
  #creation de X4, p4 variables de bruits
  p4 <- 30
  X4 <- mvtnorm::rmvnorm(n,mu,sigma^2*diag(p4))
  
  #binarisation des colonnes de X2 et des colonnes c(3,14,15,16,17,18,27,28,29,30) de X3
  X2tilde <- matrix(0,n,p1+p2+p3)
  for (i in 1:ncol(X2)) 
  {
    X2tilde[which(X2[,i] > stats::median(X2[,i])),i] <- 1
  }
  X3tilde <- X3
  X3tilde[,c(3,14,15,16,17,18,27,28,29,30)] <- 0
  for (i in c(3,14,15,16,17,18,27,28,29,30)) 
  {
    X3tilde[which(X3[,i] > stats::median(X3[,i])),i] <- 1
  } 
  
  #for (i in 1:ncol(X2)) X2tilde[,i]<- as.factor(cut(X2[,i],c(-Inf,quantile(X2[,i],1/2),Inf)))
  #for (i in 1:ncol(X2)) levels(X2tilde[,i])=paste("v",i,"_",1:nlevels(X2tilde[,i]),sep="")
  
  #construction de y
  beta <- c(1,2,3,rep(1/5,5),rep(2/5,5),rep(3/5,5),rep(0,12))
  y <- X1%*%beta+X2tilde%*%beta+X3tilde%*%beta - 9
  
  #création des noms de variable
  X1 <- data.frame(X1)
  X2tilde <- data.frame(X2tilde)
  X3tilde <- data.frame(X3tilde)
  X4 <- data.frame(X4)
  
  for (i in 1:p1) 
  {
    colnames(X1)[i] <- paste("numS",i,sep="")
    colnames(X2tilde)[i] <- paste("categS",i,sep="")
    colnames(X3tilde)[i] <- paste("mixedS",i,sep="")
    colnames(X4)[i] <- paste("noise",i,sep="")
  }
  
  j <- 1
  for (i in (p1+1):(p1+p2) )
  {
    colnames(X1)[i] <- paste("numL",j,sep="")
    colnames(X2tilde)[i] <- paste("categL",j,sep="")
    colnames(X3tilde)[i] <- paste("mixedL",j,sep="")
    colnames(X4)[i] <- paste("noise",i,sep="")
    j <- j+1
  }
  
  j <- 1
  for (i in (p1+p2+1):(p1+p2+p3) )
  {
    colnames(X1)[i] <- paste("numM",j,sep="")
    colnames(X2tilde)[i] <- paste("categM",j,sep="")
    colnames(X3tilde)[i] <- paste("mixedM",j,sep="")
    colnames(X4)[i] <- paste("noise",i,sep="")
    j <- j+1
  }
  
  #construction d'un grand dataframe
  
  for (i in 1:ncol(X2tilde)) X2tilde[,i] <- as.factor(X2tilde[,i])
  for (i in c(3,14,15,16,17,18,27,28,29,30)) X3tilde[,i] <- as.factor(X3tilde[,i])
  X <- data.frame(X1,X2tilde,X3tilde,X4)
  return(list(X=X,y=y))
}



#' @title Simulations of data
#' @name simu_classif
#' @param n number of observations
#' @param rho correlation between variables within a group
#' @param sigma standard deviation of the variables in the group of noise
#' @param threshold to binarize y
#' @param seed optional seed to pass to set.seed
#' @return the matrix X of explonatory variables and y the binary variable to predict
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats median
#' @export
#' @references 
#' Combining clustering of variables and feature selection using random forests: the CoV/VSURF 
#' procedure, Marie Chavent, Robin Genuer, Jerome Saracco, hal-01345840
simu_classif <- function(n=600,rho=0.9,sigma=1,threshold=0.5,seed=NULL)
{
  
  if (!is.null(seed)) set.seed(seed)
  #dans chaque tableau, 3 grpes de variables correlees de taille p1,p2,p3
  p1 <- 3
  Sigma1 <- matrix(rho,p1,p1)
  diag(Sigma1) <-rep(1,p1)
  
  p2 <- 15
  Sigma2<- matrix(rho,p2,p2)
  diag(Sigma2) <-rep(1,p2)
  
  p3 <- 12
  Sigma3<- matrix(rho,p3,p3)
  diag(Sigma3) <-rep(1,p3)
  
  Sigma <-diag(p1+p2+p3)
  Sigma[1:p1,1:p1] <- Sigma1
  Sigma[(p1+1):(p1+p2),(p1+1):(p1+p2)] <- Sigma2
  Sigma[(p1+p2+1):(p1+p2+p3),(p1+p2+1):(p1+p2+p3)] <- Sigma3
  
  #creation de X1,X2 et X3 qui formeront le groupe num, qual, mixte
  mu <- rep(0,p1+p2+p3)
  X1 <- mvtnorm::rmvnorm(n,mu,Sigma)
  X2 <- mvtnorm::rmvnorm(n,mu,Sigma)
  X3<- mvtnorm::rmvnorm(n,mu,Sigma)
  #creation de X4, p4 variables de bruits
  p4 <- 30
  X4 <- mvtnorm::rmvnorm(n,mu,sigma^2*diag(p4))
  
  #binarisation des colonnes de X2 et des colonnes c(3,14,15,16,17,18,27,28,29,30) de X3
  X2tilde <- matrix(0,n,p1+p2+p3)
  for (i in 1:ncol(X2)) 
  {
    X2tilde[which(X2[,i] > stats::median(X2[,i])),i] <- 1
  }
  X3tilde <- X3
  X3tilde[,c(3,14,15,16,17,18,27,28,29,30)] <- 0
  for (i in c(3,14,15,16,17,18,27,28,29,30)) 
  {
    X3tilde[which(X3[,i] > stats::median(X3[,i])),i] <- 1
  } 
  
  #for (i in 1:ncol(X2)) X2tilde[,i]<- as.factor(cut(X2[,i],c(-Inf,quantile(X2[,i],1/2),Inf)))
  #for (i in 1:ncol(X2)) levels(X2tilde[,i])=paste("v",i,"_",1:nlevels(X2tilde[,i]),sep="")
  
  #construction de y
  beta <- c(1,2,3,rep(1/5,5),rep(2/5,5),rep(3/5,5),rep(0,12))
  index <- X1%*%beta+X2tilde%*%beta+X3tilde%*%beta - 9
  prob <- exp(index)/(1+exp(index)) 
  y <- rep(0,n)
  #threshold <- stats::median(prob)
  y[which(prob > threshold)] <- 1
  
  #création des noms de variable
  X1 <- data.frame(X1)
  X2tilde <- data.frame(X2tilde)
  X3tilde <- data.frame(X3tilde)
  X4 <- data.frame(X4)
  
  for (i in 1:p1) 
  {
    colnames(X1)[i] <- paste("numS",i,sep="")
    colnames(X2tilde)[i] <- paste("categS",i,sep="")
    colnames(X3tilde)[i] <- paste("mixedS",i,sep="")
    colnames(X4)[i] <- paste("noise",i,sep="")
  }
  
  j <- 1
  for (i in (p1+1):(p1+p2) )
  {
    colnames(X1)[i] <- paste("numL",j,sep="")
    colnames(X2tilde)[i] <- paste("categL",j,sep="")
    colnames(X3tilde)[i] <- paste("mixedL",j,sep="")
    colnames(X4)[i] <- paste("noise",i,sep="")
    j <- j+1
  }
  
  j <- 1
  for (i in (p1+p2+1):(p1+p2+p3) )
  {
    colnames(X1)[i] <- paste("numM",j,sep="")
    colnames(X2tilde)[i] <- paste("categM",j,sep="")
    colnames(X3tilde)[i] <- paste("mixedM",j,sep="")
    colnames(X4)[i] <- paste("noise",i,sep="")
    j <- j+1
  }
  
  #construction d'un grand dataframe
  
  for (i in 1:ncol(X2tilde)) X2tilde[,i] <- as.factor(X2tilde[,i])
  for (i in c(3,14,15,16,17,18,27,28,29,30)) X3tilde[,i] <- as.factor(X3tilde[,i])
  X <- data.frame(X1,X2tilde,X3tilde,X4)
  y <- as.factor(y)
  return(list(X=X,y=y))
}


#' @title Simulations of data
#' @name simu
#' @param n number of observations
#' @param rho correlation between variables within a group
#' @param sigma standard deviation of the variables in the group of noise
#' @param name name of the data
#' @param threshold to binarize y
#' @param seed optional seed to pass to set.seed
#' @return the matrix X of explonatory variables and y the binary variable to predict
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats median
#' @export
simu <- function(n=600,rho=0.9,sigma=1,name="don",threshold=0.5,seed=NULL)
{

  if (!is.null(seed)) set.seed(seed)
  #dans chaque tableau, 3 grpes de variables correlees de taille p1,p2,p3
  p1 <- 3
  Sigma1 <- matrix(rho,p1,p1)
  diag(Sigma1) <-rep(1,p1)
  
  p2 <- 15
  Sigma2<- matrix(rho,p2,p2)
  diag(Sigma2) <-rep(1,p2)
  
  p3 <- 12
  Sigma3<- matrix(rho,p3,p3)
  diag(Sigma3) <-rep(1,p3)
  
  Sigma <-diag(p1+p2+p3)
  Sigma[1:p1,1:p1] <- Sigma1
  Sigma[(p1+1):(p1+p2),(p1+1):(p1+p2)] <- Sigma2
  Sigma[(p1+p2+1):(p1+p2+p3),(p1+p2+1):(p1+p2+p3)] <- Sigma3
  
  #creation de X1,X2 et X3 qui formeront le groupe num, qual, mixte
  mu <- rep(0,p1+p2+p3)
  X1 <- mvtnorm::rmvnorm(n,mu,Sigma)
  X2 <- mvtnorm::rmvnorm(n,mu,Sigma)
  X3<- mvtnorm::rmvnorm(n,mu,Sigma)
  #creation de X4, p4 variables de bruits
  p4 <- 30
  X4 <- mvtnorm::rmvnorm(n,mu,sigma^2*diag(p4))
  
  #binarisation des colonnes de X2 et des colonnes c(3,14,15,16,17,18,27,28,29,30) de X3
  X2tilde <- matrix(0,n,p1+p2+p3)
  for (i in 1:ncol(X2)) 
  {
  X2tilde[which(X2[,i] > stats::median(X2[,i])),i] <- 1
  }
  X3tilde <- X3
  X3tilde[,c(3,14,15,16,17,18,27,28,29,30)] <- 0
  for (i in c(3,14,15,16,17,18,27,28,29,30)) 
  {
    X3tilde[which(X3[,i] > stats::median(X3[,i])),i] <- 1
  } 
    
  #for (i in 1:ncol(X2)) X2tilde[,i]<- as.factor(cut(X2[,i],c(-Inf,quantile(X2[,i],1/2),Inf)))
  #for (i in 1:ncol(X2)) levels(X2tilde[,i])=paste("v",i,"_",1:nlevels(X2tilde[,i]),sep="")
  
  #construction de y
  beta <- c(1,2,3,rep(1/5,5),rep(2/5,5),rep(3/5,5),rep(0,12))
  index <- X1%*%beta+X2tilde%*%beta+X3tilde%*%beta - 9
  prob <- exp(index)/(1+exp(index)) 
  y <- rep(0,n)
  #threshold <- stats::median(prob)
  y[which(prob > threshold)] <- 1
  
  #création des noms de variable
  X1 <- data.frame(X1)
  X2tilde <- data.frame(X2tilde)
  X3tilde <- data.frame(X3tilde)
  X4 <- data.frame(X4)
  
  for (i in 1:p1) 
  {
    colnames(X1)[i] <- paste("numS",i,sep="")
    colnames(X2tilde)[i] <- paste("categS",i,sep="")
    colnames(X3tilde)[i] <- paste("mixedS",i,sep="")
    colnames(X4)[i] <- paste("noise",i,sep="")
  }
  
  j <- 1
  for (i in (p1+1):(p1+p2) )
  {
    colnames(X1)[i] <- paste("numL",j,sep="")
    colnames(X2tilde)[i] <- paste("categL",j,sep="")
    colnames(X3tilde)[i] <- paste("mixedL",j,sep="")
    colnames(X4)[i] <- paste("noise",i,sep="")
    j <- j+1
  }
  
  j <- 1
  for (i in (p1+p2+1):(p1+p2+p3) )
  {
    colnames(X1)[i] <- paste("numM",j,sep="")
    colnames(X2tilde)[i] <- paste("categM",j,sep="")
    colnames(X3tilde)[i] <- paste("mixedM",j,sep="")
    colnames(X4)[i] <- paste("noise",i,sep="")
    j <- j+1
  }
  
  #construction d'un grand dataframe
 
  for (i in 1:ncol(X2tilde)) X2tilde[,i] <- as.factor(X2tilde[,i])
  for (i in c(3,14,15,16,17,18,27,28,29,30)) X3tilde[,i] <- as.factor(X3tilde[,i])

  if (name == "dontest") {
    Xt <- data.frame(X1, X2tilde, X3tilde, X4)
    yt <- as.factor(y)
    save(yt, Xt, file = paste("data/", name, ".rda", sep = ""))
  }
  
  else {
    X <- data.frame(X1,X2tilde,X3tilde,X4)
    y <- as.factor(y)
    save(y, X, index,file = paste("data/", name, ".rda", sep = ""))
  }
}