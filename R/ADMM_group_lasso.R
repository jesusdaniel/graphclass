# Solving G FLSA with penalty D given and weights w for lasso
# min 1/2 ||y-B||_2^2  + (w'*|B|) + lambda2*||DB||_1
# Delta = D'D
ADMM_grouplasso_weights <- function(y, D, omega1, omega2, 
                               MAX_ITER = 200, TOL = 10^(-7),
                               rho = 1, features=NULL,
                               beta_start = NULL, NODES) {
  require(Matrix)
  n = length(y) 
  m = nrow(D)
  if(is.null(beta_start)){
    beta = y
  }else{
    beta = beta_start
  }
  # Check if soft_thresholding(Y,w) == 0, best solution is zero
  soft_beta = soft_thresholding.c(y,omega1)
  q = soft_beta;  r = D%*%q
  if(max(abs(soft_beta))==0)
    return(list(beta = as.vector(soft_beta), q = q, r = r,
                iter = 0, conv_crit = 0,best_beta = soft_beta))
  qk = q;   rk = r
  u = rep(0,n);  v = rep(0,m)
  iter = 1
  phi_beta_k = 0.5*norm_s(y-beta)^2+omega1*sum(abs(beta))+omega2*gl_penalty.c(D%*%beta,NODES)
  conv_crit = Inf
  sk  = Inf;   resk = Inf
  # MONITOR BEST BETA  ----------------------------------------------------------
  best_beta = beta;  best_phi = phi_beta_k;  is_best_end = T
  #while((conv_crit>1e-04) | (resk > TOL | sk > TOL) & (iter <= MAX_ITER)) {
  while((resk > TOL | sk > TOL) & (iter <= MAX_ITER)) {
    aux = y-u + rho*q + Matrix::crossprod(D,rho*r-v)#<-----------
    beta = (aux)/(1+3*rho)
    #update q
    q = soft_thresholding.c(beta+u/rho,(1/rho)*omega1)
    # update r
    Dbeta = D%*%beta#<-----------
    Dbetavrho = as.vector(Dbeta+ v/rho)#<-----------
    r = l2_Db_soft_thresholding.c(Dbetavrho,omega2/rho, NODES)
    u = u + rho*(beta-q)#<-----------
    v = v + rho*(Dbeta-r)#<-----------    
    # Update convergence criterias --------------------------------------
    phi_beta_k1 = as.numeric( 1/2 * norm_s(beta-y)^2 + 
                                omega1*sum(abs(beta)) +omega2*gl_penalty.c(Dbeta, NODES))
    sk = rho * (max(abs(q-qk))+max(abs(Matrix::crossprod(D,r-rk))))#<-----------
    res1k = norm_s(beta-q);   res2k = norm_s(Dbeta-r)
    resk = res1k + res2k
    qk = q;   rk = r
    conv_crit = abs(phi_beta_k1 - phi_beta_k)/phi_beta_k
    phi_beta_k = phi_beta_k1
    # MONITOR BEST BETA  ----------------------------------------------------------
    if(phi_beta_k1<best_phi) {best_beta = beta;  best_phi = phi_beta_k;  is_best_end = T}else{is_best_end = F}
    #cat("  iter", iter,"-- FLSA\tConv.crit.",conv_crit,"||Dbeta-r||=",res2k,
    #        "||beta-q||",res1k, "||sk||",sk, "phi",phi_beta_k1,"\n")
    iter = iter+1
  }
  # ANALYZE FINAL SOLUTION --------------------------------------
  # selects best beta between final beta, best_beta, and q
  beta_q = beta; beta_q[which(q==0)] = 0
  phi_beta_q = as.numeric(1/2*norm_s(beta_q-y)^2 + omega1*sum(abs(beta_q)) + 
                            omega2*gl_penalty.c(D%*%beta_q, NODES))
  whichm = which.min(c(phi_beta_k1,best_phi,phi_beta_q))
  if(whichm==1) {
    best_beta = beta
  }else if(whichm==3){best_beta = beta_q}
  #cat("  iter", iter, "phi",phi_beta_k1,"phibq",phi_beta_q,"spars. ", 
  #    sum(as.vector(best_beta)==0)/(length(beta)),"sparsq.",sum(q==0),"\n")
  return(list(beta = as.vector(beta), q = q, r = r,iter = iter, 
              conv_crit = conv_crit,best_beta = as.vector(best_beta)))
}
