# Fits a logistic group lasso and returns list with optimal value and information
logistic_union_group_lasso_ridge <- function(X,Y,D, lambda1, lambda2, gamma,
                                       jobID = "NULL", verbose = F,
                                       beta_start = NULL, b_start = NA, 
                                       NODES = 264, MAX_ITER = NA, 
                                       CONV_CRIT = 1e-06, MAX_TIME = 10800) {
  #browser()
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
  grad_f <- function(Xbeta,b, beta) 
    D%*%(-Matrix::crossprod(X,Y/(1+exp(Y*(Xbeta+b))))/n) + gamma*beta
  # F evaluation-----------------------------------------------------------------------------
  f = function(Xbeta,b, beta)
    sum(log(1+exp(-Y*(Xbeta+b))))/n + gamma/2*Matrix::crossprod(beta)
  penalty = function(beta,b)
    lambda1*sum(abs(beta)) + lambda2*gl_penalty.c(beta, NODES)
  
  #Proximal and b step-----------------------------------------------------------------------
  proximal <- function(u,lambda,beta_startprox = NULL,tol = 1e-07) {
    browser
    if(lambda2>0){
      return(list(x = l1l2_proximal.c(u, lambda*lambda1, lambda*lambda2, NODES)))
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
  if(is.null(beta_start))  beta_start = rep(0,p) ###2*p??
  if(is.na(b_start))  b_start = 0
  
  optimal = fista_line_search_union_ridge(proximal,b_step,f,grad_f,penalty,beta_start,b_start,X,Y,D,jobID,verbose,
                              MAX_STEPS = MAX_ITER, TOLERANCE = CONV_CRIT,MAX_TIME = MAX_TIME, NODES = NODES)
  return(optimal)
}