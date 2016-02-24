# Fits a logistic group lasso and returns list with optimal value and information
logistic_union_group_lasso <- function(X,Y,D, lambda1, lambda2, 
                                       jobID = "NULL", verbose = F,
                                       beta_start = NULL, b_start = NA, 
                                       NODES = 264, MAX_ITER = NA, 
                                       CONV_CRIT = 1e-06, MAX_TIME = 10800) {
  rho = 1 #---------------- Controls the speed of convergence
  n = length(Y); p = dim(X)[2]
  YX = Y*X
  XY = t(YX)
  
  # Define functions wrt X,Y,D---------------------------------------------------------------
  # B derivative, B hessian------------------------------------------------------------------
  ## Note that 1/n is included to make them more stable
  #Derivative
  b_derivative = function(Xbeta,b)
    sum(-Y/(1+exp(Y*(Xbeta+b))))/n
  #Hessian
  b_hessian = function(Xbeta,b)
    sum(1/(exp(-Y*(Xbeta+b))+exp(Y*(Xbeta+b))+2))/n
  
  # Beta derivative -------------------------------------------------------------------------
  grad_f <- function(Xbeta,b) 
    D%*%(-crossprod(X,Y/(1+exp(Y*(Xbeta+b))))/n)
  # F evaluation-----------------------------------------------------------------------------
  f = function(Xbeta,b)
    sum(log(1+exp(-Y*(Xbeta+b))))/n
  penalty = function(beta,b)
    lambda1*sum(abs(beta)) + lambda2*gl_penalty.c(beta)
  
  #Proximal and b step-----------------------------------------------------------------------
  proximal <- function(u,lambda,beta_startprox = NULL,tol = 1e-07) {
    browser
    if(lambda2>0){
      return(list(x = l1l2_proximal.c(u, lambda*lambda1, lambda*lambda2)))
    }else if(lambda1>0){
      return(list(x = sign(u)*pmax(as.double(abs(u))-(lambda1*lambda),0)))
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
      b_n = b_n - b_deriv/b_hessian(Xbeta,b_n)
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
  
  optimal = fista_line_search_union(proximal,b_step,f,grad_f,penalty,beta_start,b_start,X,Y,D,jobID,verbose,
                              MAX_STEPS = MAX_ITER, TOLERANCE = CONV_CRIT,MAX_TIME = MAX_TIME)
  return(optimal)
}

library(compiler)
norm_s = function(x) sqrt(sum(x^2))
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
gl_penalty.c = function(b) {
  gl = 0;
  .C("group_lasso_node_penalty",Db=as.double(b),N = as.integer(NODES), gl = as.double(gl))$gl
}
l2_soft_thresholding.c = function(x, lambda) {
  n = length(x)
  .C("l2_soft_thresholding", x = as.double(x), lambda = as.double(lambda), n = as.integer(n))$x
}

l1l2_proximal.c = function(beta, lambda1, lambda2) {
  .C("l1l2_node_soft_thresholding", b = as.double(beta), N=as.integer(NODES), 
     lambda1 = as.double(lambda1), lambda2 = as.double(lambda2))$b
}


l2_Db_soft_thresholding.c = function(Db,lambda) {  
  n = NODES
  .C("l2_node_soft_ghtesholding", 
     Db = as.double(Db), N=as.integer(n), lambda = as.double(lambda))$Db
}