# Fits a logistic group lasso and returns list with optimal value and information
logistic_group_lasso <- function(X,Y, D, lambda1, lambda2, 
                                 jobID = "NULL",verbose = F,
                                 beta_start, b_start, 
                                 NODES, MAX_ITER, 
                                 CONV_CRIT = 1e-06, MAX_TIME = Inf) {
  rho = 1 #---------------- Controls the speed of convergence
  n = length(Y); p = dim(X)[2]
  YX = Y*X
  XY = t(YX)
  # Define functions wrt X,Y,D---------------------------------------------------------------
  # B derivative, B hessian------------------------------------------------------------------
  #Derivative
  b_derivative = function(Xbeta,b)
    sum(-Y/(1+exp(Y*(Xbeta+b))))/n
  #Hessian
  b_hessian = function(Xbeta,b)
    sum(1/(exp(-Y*(Xbeta+b))+exp(Y*(Xbeta+b))+2))/n
  
  # Beta derivative -------------------------------------------------------------------------
  grad_f <- function(Xbeta,b) 
    -crossprod(X,Y/(1+exp(Y*(Xbeta+b))))/n
  # F evaluation-----------------------------------------------------------------------------
  f = function(Xbeta,b)
    sum(log(1+exp(-Y*(Xbeta+b))))/n
  penalty = function(beta,b)
    lambda1*sum(abs(beta)) + lambda2*gl_penalty.c(D%*%beta)

  
  #Proximal and b step-----------------------------------------------------------------------
  proximal <- function(u,lambda,beta_startprox = NULL,tol = 1e-07) {
    if(lambda2>0){
      gl = ADMM_grouplasso_weights(y=u,D,omega1 = lambda*lambda1, omega2 = lambda*lambda2,
                                   beta_start=beta_startprox, TOL = tol,rho = rho)
      return(list(x = gl$best_beta, q = gl$q, r = gl$r))
    }else if(lambda1>0){
      return(list(x = sign(u)*pmax(abs(u)-(lambda1*lambda),0)))
    }else{
      return(list(x=u))
    }
  }
  b_step <- function(Xbeta, b_start = 0) {
    TOL_B = 1e-04; MAX_S_B = 100
    b_n = b_start
    i = 0
    b_deriv = Inf    
    while(abs(b_deriv)>TOL_B & i < MAX_S_B) {
      b_deriv = b_derivative(Xbeta,b_n)
      b_n = b_n - b_deriv/(b_hessian(Xbeta,b_n)+1*(abs(b_deriv/b_hessian(Xbeta,b_n))>100))
      if(abs(b_n)>1000){
        warning("The intercept is too big: ",b_n,"\n  Derivative = ",b_deriv,
                "\n  Hessian = ", b_hessian(Xbeta,b_n))
      }
      i = i+1
    }
    return(b_n)
  }
  
  #--------------------------------------------------------------------------------------
  if(is.null(beta_start))  beta_start = rep(0,p)
  if(is.na(b_start))  b_start = 0

  optimal = fista_line_search(proximal,b_step,f,grad_f,penalty,beta_start,b_start,X,Y,D,jobID,verbose,
                              MAX_STEPS = MAX_ITER, TOLERANCE = CONV_CRIT,MAX_TIME = MAX_TIME)
  return(optimal)
}

library(compiler)
norm_s = function(x) sqrt(sum(x^2))
logistic_group_lasso = cmpfun(logistic_group_lasso, options = list(optimize = 3))
norm_s = cmpfun(norm_s, options = list(optimize = 3))

# Things it should return:
# - Beta, b (obviously)
# - Conv crit, iterations at the end
# - Objective function in each iteration
# - Conv_criteria in each iteration
# - Beta at the end of the path (can differ)
# - tk step length

if(is.loaded("pmax2", PACKAGE="Pmax2.so"))   dyn.unload("Pmax2.so")

dyn.load("Pmax2.so")

Pmax2.C = function(x) {
  n = length(x)
  .C("pmax2", x = as.double(x), n = as.integer(n))$x
} 

soft_thresholding.c = function(x, lambda) {
  n = length(x)
  .C("soft_thresholding", x = as.double(x), lambda = as.double(lambda), n = as.integer(n))$x
}

soft_thresholdingl2 = function(x, lambda) {
  if(norm_s(x)<lambda) {
    return(rep(0,length(x)))
  }
  else
    return(x*(1-lambda/norm_s(x)))
}
gl_penalty.c = function(Db) {
  browser()
  gl = 0;
  .C("group_lasso_node_penalty",Db=as.double(Db),N = as.integer(NODES), gl = as.double(gl))$gl
}
l2_soft_thresholding.c = function(x, lambda) {
  n = length(x)
  .C("l2_soft_thresholding", x = as.double(x), lambda = as.double(lambda), n = as.integer(n))$x
}

l2_Db_soft_thresholding.c = function(Db,lambda) {  
  n = NODES
  .C("l2_node_soft_ghtesholding", 
     Db = as.double(Db), N=as.integer(n), lambda = as.double(lambda))$Db
}