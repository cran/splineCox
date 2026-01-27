#' @title Maximum likelihood estimation for spline copula parameter matrix
#' @description
#' Performs maximum likelihood estimation of the 5x5 spline copula
#' parameter matrix \eqn{R} using an EM-type algorithm based on
#' M-spline basis functions.
#'
#' @param U Numeric vector of length \eqn{n} with values in (0, 1).
#' @param V Numeric vector of length \eqn{n} with values in (0, 1).
#' @param E1 Positive scalar specifying the convergence tolerance
#'   for the outer EM iteration updating \eqn{R}.
#' @param E2 Positive scalar specifying the convergence tolerance
#'   for the inner iteration updating auxiliary parameters.
#' @param R0 Optional 5x5 numeric matrix providing initial values
#'   for \eqn{R}. If \code{NULL}, a data-driven initialization is used.
#'
#' @return A list containing:
#' \describe{
#'   \item{R}{Estimated 5x5 spline copula parameter matrix.}
#'   \item{R0}{Initial matrix used to start the algorithm.}
#'   \item{convergence}{Matrix recording the iteration number and
#'     corresponding log-likelihood values.}
#' }
#'
#' @details
#' The algorithm alternates between updating the spline copula
#' parameter matrix and auxiliary parameters until convergence
#' criteria are satisfied.
#'
#' @examples
#' n <- 100
#' R <- matrix(c(1,0,0,0,0,
#'               0,2,0,0,0,
#'               0,0,2,0,0,
#'               0,0,0,2,0,
#'               0,0,0,0,1)/8, 5, 5, byrow = TRUE)
#' out <- spline.copula.simu(n, R = R)
#' fit <- spline.copula.MLE(out$U, out$V)
#'
#' @importFrom stats optimize
#' 
#' @export

spline.copula.MLE=function(U,V,E1=0.001,E2=0.001,R0=NULL){
  
  n=length(U)
  b=c(1/8,1/4,1/4,1/4,1/8)
  MU=M.spline(U,xi1=0,xi3=1)
  MV=M.spline(V,xi1=0,xi3=1)
  if(is.null(R0)){ R0=(t(MU)%*%MV)*(b%*%t(b))/n }
  
  E.func=function(R){
    Z=matrix(0,5,5)
    for(i in 1:n){ 
      cR=as.numeric(t(MU[i,])%*%R%*%(MV[i,]))
      Z=Z+R*MU[i,]%*%t(MV[i,])/cR
    }
    return(Z/n)
  }
  
  l.func=function(R){ 
    cR=diag(MU%*%R%*%t(MV))
    return(sum(log(cR)))
  }
  
  lambda.func=function(mu,R){
    lambda=rep(0,5)
    E=E.func(R)
    for(k in 1:5){
      f.func=function(x){ (sum(E[k,]/(x+mu))-b[k])^2 }
      sol=optimize(f.func,interval=c(-min(mu),50))
      lambda[k]=sol$minimum
    }
    return(lambda)
  }
  
  mu.func=function(lambda,R){
    mu=rep(0,5)
    E=E.func(R)
    for(l in 1:5){
      f.func=function(x){ (sum(E[,l]/(lambda+x))-b[l])^2 }
      sol=optimize(f.func,interval=c(-min(lambda),50))
      mu[l]=sol$minimum
    }
    return(mu)
  }
  
  M.func=function(R){
    num=0
    Err=10
    mu.old=rep(0,5)
    while (Err > E2) {
      lambda=lambda.func(mu.old,R)
      mu.new=mu.func(lambda,R)
      Err=max(abs(mu.new-mu.old))
      mu.old=mu.new
      num=num+1
      if(num>=100){break}
    }
    mu=mu.new
    E=E.func(R)/(matrix(lambda,5,5,byrow=FALSE)+matrix(mu,5,5,byrow=TRUE))
    num=num
    list(E=E,num=num,mu=mu,lambda=lambda)
  }
  
  #### EM ####
  Err=1
  num=0
  
  R.old=R0
  logl=l.func(R.old)
  while (Err > E1) {
    R.new=M.func(R.old)$E
    Err=max(abs(R.new-R.old))
    logl=c(logl,l.func(R.new))
    R.old=R.new
    num=num+1
    if(num>=100){break}
  }
  
  res=list(R=R.new,R0=R0,convergence=cbind(iteration=0:num,logl=logl))
  res
  
}






