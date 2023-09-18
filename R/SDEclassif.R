#' Name: SDEclassif
#' Title: A plug-in type classifier for a multiclass classification of SDE paths.
#'
#' @param Xtrain: a training sample sample which is a matrix or a data.frame whose columns represent the diffusion paths
#' @param ytrain: a second training sample which is a vector containing the labels of diffusion paths
#' @param Xtest: a test sample which is matrix or a data.frame containing in columns the trajectories to be classified
#' @param ytest: this test sample contains the respective labels of paths in the test set Xtest
#' @param NbClass: number of classes, an integer greater or equal to 2
#' @param TimeStep: time step for the discretization of the time interval
#' @param SetK: set of possible values of the dimension of the subspace spanned by the spline basis
#' @param M: degree of the spline basis
#' @param L: a constant allowing to control the coordinate vector of projections of the S.D.E. coefficient on the approximation space
#' @param paths= TRUE or FALSE: if paths=TRUE, then the function returns a plot representing the test data and the data in the predicted classes
#'
#' @return predicted labels
#' @importFrom stats D
#' @import MASS
#' @import qpdf
#' @importFrom sde sde.sim
#' @import fda
#' @importFrom rootSolve multiroot
#' @import data.table
#' @examples
#'## Diffusion model:
#'## b(x) = theta*(mu + nu x)
#'## sigma(x) = 1
#'## Class 1: theta = 1; mu = 1; nu = -1
#'## Class 2: theta = 1; mu = 0; nu = -1
#'## Class 3: theta = 2; mu = 1; nu = -1
#'## Simulation of data:
#'library(sde)
#' T = 1; x0 = 0; NbClass = 2;
#' N = 300      # (number of path)
#' n = 50       # (number of discretization steps)
#' TimeStep = 1/n; M = 3; L = log(N)
#' SetK = 2**(0:5)
#' ## Probabilities of labels:
#' p = rep(1/NbClass, NbClass)
#' N1 = round(N/2); N2 = N - N1
#' ytrain = c(rep(1,N1), rep(2,N2)); ytest = c(rep(1,50), rep(2,50));
#' ## Training sample:
#' X1 = sde.sim(t0=0,T=1,X0=x0,N=n+1,M=N1,theta=c(1,1,1),model="OU"); X1 <- as.matrix(X1)
#' X2 = sde.sim(t0=0,T=1,X0=x0,N=n+1,M=N2,theta=c(1,1,1),model="OU"); X2 <- as.matrix(X2)
#' Xtrain = cbind(X1,X2)
#' ## Test sample
#' Xtest = matrix(ncol=100, nrow=n+2)
#' tX1 = sde.sim(t0=0,T=1,X0=x0,N=n+1,M=50,theta=c(1,1,1),model="OU"); tX1 <- as.matrix(tX1)
#' tX2 = sde.sim(t0=0,T=1,X0=x0,N=n+1,M=50,theta=c(1,1,1),model="OU"); tX2 <- as.matrix(tX2)
#' Xtest = cbind(tX1,tX2)
## Prediction
#' SDEclassif(Xtrain, ytrain, Xtest, ytest, NbClass, TimeStep, SetK, M, L, paths=TRUE)
#'
#' @export
SDEclassif <- function(Xtrain, ytrain, Xtest, ytest, NbClass, TimeStep, SetK, M, L, paths=FALSE){
  #### Estimation of drift functions
  DriftList <- list()        ## List of coordinate vectors
  Kdrift <- rep(NA, NbClass)
  for(i in 1:NbClass){
    iXtrain <- Xtrain[, which(ytrain == i)]
    InfInt <- -log(ncol(iXtrain)); SupInt <- log(ncol(iXtrain));
    iZ <- Zfun(iXtrain, TimeStep)
    Kdrift[i] <- SelectDimDrift(iXtrain, 0.1, SetK, M, InfInt, SupInt, TimeStep, L)
    iB <- Bfun(iXtrain, InfInt, SupInt, Kdrift[i], M)
    DriftList[[i]] <- OptimFun(iB, iZ, Kdrift[i], M, L)
  }
  #### Estimation of the diffusion coefficient
  U <- Ufun(Xtrain, TimeStep)
  InfInt. <- -log(ncol(Xtrain)); SupInt. <- log(ncol(Xtrain));
  Kdiff <- SelectDimDiff(Xtrain, 5, SetK, M, InfInt., SupInt., TimeStep, L)
  B <- Bfun(Xtrain, InfInt., SupInt., Kdiff, M)
  SigmaVector <- OptimFun(B, U, Kdiff, M, L)
  #### Computation of probabilities
  n <- nrow(Xtrain) - 2;
  ExpFspline <- matrix(nrow = ncol(Xtest), ncol = NbClass)
  VectorPi <- matrix(nrow = ncol(Xtest), ncol = NbClass)
  VectorProba <- rep(NA, NbClass)
  for(k in 1:NbClass) VectorProba[k] <- length(which(ytrain == k))/length(ytrain);
  for(j in 1:ncol(Xtest)){
    for(k in 1:NbClass){
      kXtrain <- Xtrain[, which(ytrain == k)]
      InfInt <- -log(ncol(kXtrain)); SupInt <- log(ncol(kXtrain));
      result <- (DriftSpline(Xtest[1:n, j], DriftList[[k]], InfInt, SupInt, Kdrift[k], M, L)/DiffSpline(Xtest[1:n,j], SigmaVector, InfInt., SupInt., Kdiff, M, L))*diff(Xtest[,j])[-1]-(TimeStep/2)*(DriftSpline(Xtest[1:n,j], DriftList[[k]], InfInt, SupInt, Kdrift[k], M, L)*DriftSpline(Xtest[1:n,j], DriftList[[k]], InfInt, SupInt, Kdrift[k], M, L)/DiffSpline(Xtest[1:n,j], SigmaVector, InfInt., SupInt., Kdiff, M, L));
      ExpFspline[j,k] <- VectorProba[k]*exp(sum(result))
    }
    VectorPi[j,] <- ExpFspline[j,]/sum(ExpFspline[j,])
  }
  #### Prediction of labels
  PredClass <- apply(VectorPi, 1, function(x) which.max(x))
  #### Plot
  if(paths==TRUE){
    n = nrow(Xtrain)-2
    graphics::par(mfrow=c(1,2))
    #### Test data
    graphics::plot(Xtest[,1],type='l',col=ytest[1],main="True labels",ylim=range(Xtest),xlab="",ylab="");
    for(i in 2:ncol(Xtest)) graphics::lines(Xtest[,i],type='l',col=ytest[i])
    #### Predicted data
    graphics::plot(Xtest[,1],type='l',col=PredClass[1],main="predicted labels",ylim=range(Xtest),xlab="",ylab="",las=1);
    for(i in 2:ncol(Xtest)) graphics::lines(Xtest[,i],type='l',col=PredClass[i])
  }

  return(PredClass)
}







###############################################################################################################################
###########################################################################################################################
############################## Used functions
###########################################################################################################################
###############################################################################################################################

Zfun<-function(Mx,delta){
  if(ncol(Mx)==1) result <- diff(Mx[,1])[-1]/delta
  if(ncol(Mx)>1){
    result <- diff(Mx[,1])[-1]/delta
    for (i in 2:dim(Mx)[2]){
      result <- c(result, (diff(Mx[,i])[-1])/delta)
    }
  }
  return(result)
}


Ufun<-function(Mx,delta){
  if(ncol(Mx)==1) result <- diff(Mx[,1])[-1]/delta
  if(ncol(Mx)>1){
    result <- diff(Mx[,1])[-1]/delta
    for (i in 2:dim(Mx)[2]){
      result <- c(result, (diff(Mx[,i])[-1])/delta)
    }
  }

  return(delta*(result**2))
}


## The function lab.fun() return vector of 0 and 1

labfun <- function(x, InfInt, SupInt) sapply(x, function(a) (a >= InfInt)*(a <= SupInt))

## Function bound.fun() replace each component of x that is not in the interval [InfInt, SupInt] by InfInt

BoundFun <- function(x, InfInt, SupInt) sapply(x, function(a) (a - InfInt)*(a >= InfInt)*(a <= SupInt) + InfInt)

## The function bspline() evaluates the B-spline basis' functions

bspline <- function(x, InfInt, SupInt, K, M){
  y <- BoundFun(x, InfInt, SupInt)
  id <- diag(labfun(x, InfInt, SupInt))
  dm <- K + M
  bs_basis <- create.bspline.basis(rangeval=c(InfInt, SupInt),nbasis = dm)
  v.basis <- eval.basis(y, basisobj = bs_basis)
  if(length(id)!=0) result <- id%*%v.basis
  if(length(id)==0) result <- 0*v.basis

  return(result)
}


Bfun <- function(Mx, InfInt, SupInt, K, M){
  N = ncol(Mx);  n = nrow(Mx)
  df.X <- as.data.frame(Mx[2:(n-1),])
  l.X <- as.list(df.X)
  B. <- lapply(l.X,function(x) bspline(x, InfInt, SupInt, K, M))
  B.. <- lapply(B.,function(x) as.data.frame(x))
  B <- rbindlist(B..)

  return(as.matrix(B))
}


## Computation of the coordinate vectors
ridge <- function(x, matrix.B, vector.Z, K, M, L){
  upper_bd = (K+M)*L;
  u = ginv(t(matrix.B)%*%matrix.B + x*diag(K+M))%*%t(matrix.B)%*%vector.Z
  result = sum(u^2) - upper_bd

  return(result)
}


OptimFun <- function(matrix.B, vector.Z, K, M, L){
  upper_bound <- sqrt(L*(K+M))
  estimator <- ginv(t(matrix.B)%*%matrix.B)%*%t(matrix.B)%*%vector.Z;
  norm_estimator <- sqrt(sum(estimator^2))
  Px = t(matrix.B) %*% matrix.B
  if(det(Px) != 0 & norm_estimator <= upper_bound){
    result <- estimator
  }else{
    ridge.bis <- function(x) ridge(x,matrix.B,vector.Z,K,M,L)
    root = multiroot(ridge.bis, c(0.01));
    lambda <- root$root
    result <- ginv(t(matrix.B)%*%matrix.B+lambda*diag(K+M))%*%t(matrix.B)%*%vector.Z;
  }

  return(result)
}


## Selection of the dimension - drift coefficient
gamma_pen<-function(x,c,K,M,Z,B,n,N) (1/(n*N))*sum((Z-B%*%x)**2)+c*(log(N))*(K+M)/N;

SelectDimDrift <- function(X, c, SetKspline, M, lower, upper, delta, Lconst){
  N <- ncol(X); n <- nrow(X) - 2; N. <- ncol(D); Z <- Zfun(X, delta)
  B.ls<-lapply(SetKspline, function(x) Bfun(X, lower, upper, x, M))
  a.ls<-lapply(1:length(SetKspline), function(x) OptimFun(B.ls[[x]], Z, SetKspline[x], M, Lconst))
  gpen.vec<-sapply(1:length(SetKspline), function(x) gamma_pen(a.ls[[x]], c, SetKspline[x], M, Z, B.ls[[x]], n, N))
  i.min<-which.min(gpen.vec)
  K.ch<-SetKspline[i.min]

  result=K.ch

  return(result)
}


## drift estimator

sign <- function(x) -1*(x<=0)+1*(x>0)
EstimBound <- function(x, Lconst) sapply(x,function(a) a*(abs(a)<=sqrt(Lconst))+sign(a)*sqrt(Lconst)*(abs(a)>sqrt(Lconst)))

DriftSpline <- function(x, ach, lower, upper, K, M, Lconst){
  lab <- labfun(x, lower, upper)
  idmat <- diag(lab)
  bMat = bspline(x, lower, upper, K, M)                # Use of the bspline functions
  if(length(idmat)!=0){
    bsMat <- idmat%*%bMat;
    bvalues <- bsMat%*%ach
  }
  if(length(idmat) == 0) bvalues <- 0*ach;
  bvalues<-EstimBound(bvalues,Lconst)

  return(as.vector(bvalues))
}


#### Selection of the dimension - diffusion coefficient

gamma_pen2 <- function(x, c, K, M, U, B, n, N) (1/(n*N))*sum((U-B%*%x)**2)+c*log(N)*(K+M)/(N*n);

SelectDimDiff <- function(X,c,set.K,M,lower,upper,delta,L.N){
  N<-ncol(X); n <- nrow(X) - 2; N.<-ncol(D); U<-Ufun(X,delta)
  B.ls<-lapply(set.K,function(x) Bfun(X,lower,upper,x,M))
  L.N<-log(N)
  a.ls<-lapply(1:length(set.K),function(x) OptimFun(B.ls[[x]],U,set.K[x],M,L.N))
  gpen.vec<-sapply(1:length(set.K),function(x) gamma_pen2(a.ls[[x]],c,set.K[x],M,U,B.ls[[x]],n,N))
  i.min<-which.min(gpen.vec)
  K.ch<-set.K[i.min]

  result=K.ch

  return(result)
}

## diffusion coefficient estimator

zerofun <- function(x) sapply(x,function(a) (a-0.01)*(a>0)+0.01)

DiffSpline <- function(x, alpha.ch, lower, upper, K, M, L.N){
  lab <- labfun(x, lower, upper)
  id.mat <- diag(lab)
  bsMat = bspline(x, lower, upper, K, M)                # Use of the bspline functions
  if(length(id.mat)!=0) bs.Mat <- id.mat%*%bsMat
  if(length(id.mat)==0) bs.Mat <- 0*bsMat
  s.values = bs.Mat%*%alpha.ch
  result <- zerofun(s.values); result <- EstimBound(result,L.N)

  return(as.vector(result))
}
