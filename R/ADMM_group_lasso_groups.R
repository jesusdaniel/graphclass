# Solving G FLSA with penalty D given and weights w for lasso
# min 1/2 ||y-B||_2^2  + (w'*|B|) + lambda2*(G_1||D_1B||_2 + ... + G_g||D_gB||_2)
# Delta = D'D
ADMM_grouplasso_weights_groups <- function(y, D_list, omega1, omega2, G_penalty_factors,
                                    MAX_ITER = 200, TOL = 10^(-7),
                                    rho = 1, features=NULL,
                                    beta_start = NULL, NODES) {
  require(Matrix)
  n = length(y) 
  G = length(D_list)
  if(is.null(beta_start)){
    beta = y
  }else{
    beta = beta_start
  }
  # Check if soft_thresholding(Y,w) == 0, best solution is zero
  soft_beta = soft_thresholding.c(y,omega1)
  q = soft_beta;  
  R = tcrossprod(q, rep(1, G))
  if(max(abs(soft_beta))==0)
    return(list(beta = as.vector(soft_beta), q = q, r = R,
                iter = 0, conv_crit = 0,best_beta = soft_beta))
  qk = q;   Rk = R
  u = rep(0,n);  V = matrix(0, nrow = n, ncol = G)
  iter = 1
  beta2 = beta^2
  phi_beta_k = 0.5*norm_s(y-beta)^2+omega1*sum(abs(beta))+omega2*crossprod(G_penalty_factors, sapply(D_list, function(D) sum(sqrt(D%*%beta2))))
  conv_crit = Inf
  sk  = Inf;   resk = Inf
  # MONITOR BEST BETA  ----------------------------------------------------------
  best_beta = beta;  best_phi = phi_beta_k;  is_best_end = T
  #while((conv_crit>1e-04) | (resk > TOL | sk > TOL) & (iter <= MAX_ITER)) {
  while((resk > TOL | sk > TOL) & (iter <= MAX_ITER)) {
    aux = y-u + rho*q + rho*apply(R,1,sum) - apply(V,1, sum)
    beta = aux/(1+(G+1)*rho)
    #update q
    q = soft_thresholding.c(beta+u/rho,(1/rho)*omega1)
    # update r
    BetaVrho <- tcrossprod(beta, rep(1, G)) + V/rho
    R = l2_groups_softhtresholding(BetaVrho, D_list, omega2/rho, G_penalty_factors)
    u = u + rho*(beta-q)#<-----------
    V = V + rho*(tcrossprod(beta, rep(1, G))-R)#<-----------    
    # Update convergence criterias --------------------------------------
    beta2 = beta^2
    phi_beta_k1 = as.numeric( 1/2 * norm_s(beta-y)^2 + 
                                omega1*sum(abs(beta)) + omega2*crossprod(G_penalty_factors, sapply(D_list, function(D) sum(sqrt(D%*%beta2)))))
    sk = rho * (max(abs(q-qk))+max(abs(R- Rk)))#<-----------
    res1k = norm_s(beta-q);   res2k = norm(tcrossprod(beta, rep(1, G))-R)
    resk = res1k + res2k
    qk = q;   Rk = R
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
  beta_q_2 = beta_q^2
  phi_beta_q = as.numeric(1/2*norm_s(beta_q-y)^2 + omega1*sum(abs(beta_q)) + 
                            omega2*crossprod(G_penalty_factors, sapply(D_list, function(D) sum(sqrt(D%*%beta_q_2)))))
  whichm = which.min(c(phi_beta_k1,best_phi,phi_beta_q))
  if(whichm==1) {
    best_beta = beta
  }else if(whichm==3){best_beta = beta_q}
  #cat("  iter", iter, "phi",phi_beta_k1,"phibq",phi_beta_q,"spars. ", 
  #    sum(as.vector(best_beta)==0)/(length(beta)),"sparsq.",sum(q==0),"\n")
  return(list(beta = as.vector(beta), q = q, r = R,iter = iter, 
              conv_crit = conv_crit,best_beta = as.vector(best_beta)))
}