# Solving G FLSA with penalty D given 
# min 1/2 ||y-B||_2^2  + (w*|B|) + lambda2*||DB||_1
# Delta = D'D
ADMM_fused_lasso <- function(y, D, omega1, omega2, G_inv,
                                    MAX_ITER = 200, TOL = 10^(-7),
                                    rho = 0.1, features=NULL,
                                    beta_start = NULL) {
  #browser()
  require(Matrix)
  n = length(y); m = nrow(D)
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
  phi_beta_k = 0.5*norm_s(y-beta)^2+omega1*sum(abs(beta))+omega2*sum(abs(D%*%beta))
  #t3 = proc.time()[3]
  conv_crit = Inf; sk  = Inf;   resk = Inf
  # MONITOR BEST BETA  ----------------------------------------------------------
  best_beta = beta;  best_phi = phi_beta_k;  is_best_end = T
  #while((conv_crit>1e-10) & (resk > TOL | sk > TOL) & (iter <= MAX_ITER)) {
  while((resk > TOL | sk > TOL) & (iter <= MAX_ITER)) {
    aux = y-u + rho*q + crossprod(D,rho*r-v)#<-----------
    
    beta = G_inv$multiply_by_inverse(aux,rho)
    Dbeta = D%*%beta#<-----------
    #update q
    q = soft_thresholding.c(beta+u/rho,(1/rho)*omega1)
    # update r
    Dbetavrho = as.vector(Dbeta+ v/rho)#<-----------
    
    r = soft_thresholding.c(Dbetavrho, omega2/rho)
    
    u = u + rho*(beta-q)#<-----------
    v = v + rho*(Dbeta-r)#<-----------
    # Update convergence criterias --------------------------------------
    phi_beta_k1 = as.numeric(1/2*norm_s(beta-y)^2 + omega1*sum(abs(beta)) + 
                               omega2*sum(abs(Dbeta)))
    sk = rho * (max(abs(q-qk))+max(abs(crossprod(D,r-rk))))#<-----------
    res1k = norm_s(beta-q);   res2k = norm_s(Dbeta-r)
    resk = res1k + res2k
    qk = q;   rk = r
    conv_crit = abs(phi_beta_k1 - phi_beta_k)/phi_beta_k
    phi_beta_k = phi_beta_k1
    #     ### Update rho ----------------------------------------------------------
    #     sk = mu*norm_s(t(D)%*%(bk1-b))
    #     rk = norm_s(Dbeta-bk1) 
    #if(resk > 10*sk) {
    #  rho = 2*rho
    #}else if(sk > 10*resk) {
    #  rho = rho/2
    #}
    # MONITOR BEST BETA  ----------------------------------------------------------
    if(phi_beta_k1<best_phi) {best_beta = beta;  best_phi = phi_beta_k;  is_best_end = T}else{is_best_end = F}
    #cat("  iter", iter,"-- FLSA\tConv.crit.",conv_crit,"||Dbeta-r||=",res2k,
    #        "||beta-q||",res1k, "||sk||",sk, "phi",phi_beta_k1,"\n")
    iter = iter+1
  }
  # ANALYZE FINAL SOLUTION --------------------------------------
  # selects best beta between final beta, best_beta, and q
  phi_q = as.numeric(1/2*norm_s(q-y)^2 + omega1*sum(abs(q)) + omega2*sum(abs(D%*%q)))
  whichm = which.min(c(phi_beta_k1,best_phi,phi_q))
  if(whichm==1) {
    best_beta = beta
  }else if(whichm==3){best_beta = q}
  #t5 = proc.time()[3]
  #cat("  iter", iter, "phi",phi_beta_k1,"spars. ", sum(as.vector(best_beta)==0)/(length(beta)),"\n")
  return(list(beta = as.vector(beta), q = q, r = r,iter = iter, 
              conv_crit = conv_crit,best_beta = as.vector(best_beta)))
}
