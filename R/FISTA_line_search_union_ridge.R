
fista_line_search_union_ridge <- function(proximal_f,b_step,f,grad_f,penalty, x_start, b_start,
                              X,Y,D,
                              ID = "NULL", verbose = F,
                              MAX_STEPS = 300, TOLERANCE = 1e-06,
                              MAX_TIME = 10800, NODES = 263) {
  crits = c();  fs = c();  tks = c();
  # Initialize -------------------------------------------------------------------
  xk1 = x_start;  xk = x_start
  criterion = Inf
  crit_f = Inf
  k = 0;  tk = 0.125;  beta_step = 1/2
  Xbeta = X%*%crossprod(D,x_start)
  b = b_step(Xbeta,b_start)
  best_beta = x_start;  best_b = b;   
  best_f = f(Xbeta,b, x_start) + penalty(x_start,b);  best_prox = NULL
  newf = as.double(best_f)
  if(verbose) cat(sprintf("---%s - Iteration %d.  Crit_f = -.  F = %.4f\n",ID,0,newf))
  time_start = proc.time()[1]
  crit_f_1  = crit_f
  beta_path  = list()
  beta_path[[1]] = xk1
  is_best_end = FALSE
  while(k<=5 | ((crit_f > TOLERANCE  & criterion > TOLERANCE) & (k < MAX_STEPS) &
                  proc.time()[1]-time_start < MAX_TIME)){
    t1 = proc.time()[3]
    crit_f = newf
    is_best_end = FALSE
    k = k+1
    xk_1 = xk;    xk = xk1
    y = xk + (k-2)/(k+1) * (xk-xk_1)
    ######## line search ###################
    repeat {
      z=y-tk*grad_f(X%*%crossprod(D,y),b, y)
      prox = proximal_f(z,tk,y,tol = 1e-02/k)
      z = prox$x
      Xbeta = X%*%crossprod(D,z)
      Xy = X%*%crossprod(D,y)
      if(as.double(f(Xbeta,b, z)) <= as.double(f(Xy,b,y) +t(grad_f(Xy,b,y))%*%(z-y) + 1/(2*tk) * norm_s(z-y)^2))
        break #tk = 1/beta_step*tk
      #cat("  ...trying new tk\tt =",proc.time()[1]-time_start, "\n")
      tk = beta_step*tk
    }
    
    ######## line search ###################
    b = b_step(Xbeta,b)
    #cat(paste("t2: ",t2-t1, "t4:",t4-t3,"t5:",t4-t5,"t7:",t7-t6,"\n"))
    tks = c(tks,tk)
    xk1 = z;    criterion = norm_s(xk1-xk)
    crits = c(crits,criterion)
    newf = as.double(f(Xbeta,b,z) + penalty(xk1,b))
    if(crit_f>0) #crit_f = newf above
      crit_f = as.double(abs(newf - crit_f)/crit_f)
    if(as.double(newf) < as.double(best_f)) {
      is_best_end = TRUE
      best_beta = xk1;     best_b = b;   best_prox = prox
      best_f = newf
    }
    #cat("2:",time2-time1,"\n",          "3:",time3-time2,"\n",         "4:",time4-time3,"\n",        "5:",time5-time4,"\n")    
    fs = c(fs, newf)
    if(verbose)  cat(sprintf("--%s - Iter %d. Sparse = %d.  Crit_f = %.4f. tk = %.2f. F = %.4f. t=%2f. N=%d.\n",
                             ID,k,sum(xk1==0),crit_f,tk,newf,proc.time()[1]-time_start,active_nodes_union(xk1,NODES)))
    if(as.double(crit_f)==0 | as.double(abs(crit_f-crit_f_1)/crit_f)<0.1)
      tk = 2*tk
    crit_f_1 = crit_f
    beta_path[[k+1]] = xk1
  }
  if(crit_f < TOLERANCE | criterion < TOLERANCE) {
    status = 1
  }else if(k>=MAX_STEPS){
    status = 2
  }else{
    status = 3
  }
  return(list(x = xk1, b=b, crits = crits, fs = fs, L = tk, tks = tks, conv = criterion,
              best_f = best_f, best_beta = best_beta, best_b = best_b, 
              is_best_end = is_best_end, best_prox = best_prox,
              #beta_path = beta_path,
              status = status))
}
