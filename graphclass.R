# Graph class

# Adj_list is a list of symmetric matrices with 0 diagonal
# X is a matrix with rows equal to the upper triangular part of each adjacency matrix
# col(X) edges, row(X) samples
# Params: 
# - beta_start, b_start, MAX_ITER, CONV_CRIT, MAX_TIME
graphclass <- function(X = NULL, Adj_list = NULL, Y = Y, Xtest = NULL, Ytest = NULL,
                       type = "intersection",
           lambda1, lambda2, params = NULL, id = "", verbose = T, D = NULL) {
  source("ADMM_group_lasso.R")
  source("ADMM_fused_lasso.R")
  source("FISTA_line_search.R")
  source("FISTA_line_search_union.R")
  source("Logistic_group_lasso.R")
  source("Logistic_union_group_lasso.R")
  source("Logistic_fused_lasso.R")
  #browser()
  # Dependent data
  if(!is.null(X)) {
    # ncols = N*(N-1)/2
    NODES <- (1+sqrt(1+8*ncol(X)))/2
  }else{
    NODES <- ncol(Adj_list[[1]])
    X <- t(sapply(Adj_list, function(A) A[upper.tri(A)]))
  }
  lambda1 = 2*lambda1/(NODES*(NODES-1))
  lambda2 = lambda2/((NODES)*(NODES-1))
  # Independent variable
  Y <- as.numeric(Y)
  Y <- (Y-min(Y))/(max(Y)-min(Y))
  Y <- 2*Y-1
  # Normalize X for numerical stability
  minX <- min(X); maxX <- max(X)
  alpha_normalization <- (maxX-minX)
  X <- X/alpha_normalization
  lambda1 <- lambda1/alpha_normalization
  lambda2 <- lambda2/alpha_normalization
  # Create D
  if(is.null(D)) {
    if(type=="intersection"| type=="union") {
      D <- construct_D(NODES)  
    }else{
      if(type=="fusion")
        D <- construct_D_fusion(NODES)
      else{
        stop("The value of type should be one between \"intersection\", \"union\" or \"fusion\".")
      }
    }
  }
  # Check which method to use
  # Intersection penalty -----------------------------------------------
  if(type=="intersection") {
    # Check params
    if(is.null(params)){
      beta_start <- rep(0,NODES*(NODES-1)/2)
      b_start <- 0
      MAX_ITER <- 300;    CONV_CRIT <- 1e-05;   MAX_TIME = Inf
    }else{if(is.null(params$beta_start)) {
      beta_start <- rep(0,NODES*(NODES-1)/2)
      b_start <- 0
      MAX_ITER <- params$MAX_ITER;    CONV_CRIT <- params$CONV_CRIT;   MAX_TIME = params$MAX_TIME
    }}
    gl = logistic_group_lasso(X, Y, D,
                              lambda1 = lambda1, lambda2 = lambda2,
                              id, verbose, beta_start = beta_start, b_start = b_start,
                              NODES = NODES, MAX_ITER = MAX_ITER, 
                              CONV_CRIT = CONV_CRIT, 
                              MAX_TIME = MAX_TIME)
    gl_results = list()
    gl_results$beta = gl$best_beta/alpha_normalization;    gl_results$b = gl$b

    # Train fitting
    Yfit <- X%*%gl$best_beta + gl$b
    gl_results$Yfit <- exp(Yfit)/(1+exp(Yfit))
    gl_results$train_error <- 1- sum(diag(table(sign(Yfit),sign(Y))))/length(Y)
    # Test fitting
    if(!is.null(Xtest)) {
      levels(Ytest) <- c(-1, 1)
      Yfit_test <- Xtest%*%gl_results$beta + gl_results$b
      gl_results$Yfit_test <- exp(Yfit_test)/(1+exp(Yfit_test))
      gl_results$test_error = 1-sum(diag(table(sign(Yfit_test),sign(Ytest))))/length(Ytest)
    }
    gl_results$type = "intersection"
    return(gl_results)

  }else{if(type == "union"){
    # Union penalty ---------------------------------------------------
    # Check params
    if(is.null(params)){
      beta_start <- rep(0,NODES*(NODES-1))
      b_start <- 0
      MAX_ITER <- 300;    CONV_CRIT <- 1e-05;   MAX_TIME = Inf
    }else{if(is.null(params$beta_start)) {
      beta_start <- rep(0,NODES*(NODES-1))
      b_start <- 0
      MAX_ITER <- params$MAX_ITER;    CONV_CRIT <- params$CONV_CRIT;   MAX_TIME = params$MAX_TIME
    }}
    gl = logistic_union_group_lasso(X, Y, D,
                              lambda1 = lambda1, lambda2 = lambda2,
                              id, verbose, beta_start = beta_start, b_start = b_start,
                              NODES = NODES, MAX_ITER = MAX_ITER, 
                              CONV_CRIT = CONV_CRIT, 
                              MAX_TIME = MAX_TIME)
    gl_results = list()
    gl_results$beta = as.vector(gl$best_beta)/alpha_normalization;    gl_results$b = gl$b
    # Train fitting
    Yfit <- X%*%(crossprod(D,gl$best_beta )) + gl$b
    gl_results$Yfit <- exp(Yfit)/(1+exp(Yfit))
    gl_results$train_error <- 1-sum(diag(table(as.vector(sign(gl_results$Yfit-0.5)),sign(Y))))/length(Y)
    # Test fitting
    if(!is.null(Xtest)) {
      levels(Ytest) <- c(-1, 1)
      Yfit_test <- Xtest%*%crossprod(D,gl_results$beta) + gl_results$b
      gl_results$Yfit_test <- exp(Yfit_test)/(1+exp(Yfit_test))
      gl_results$test_error = 1-sum(diag(table(as.vector(sign(gl_results$Yfit_test-0.5)),sign(Ytest))))/length(Ytest)
    }
    gl_results$type = "union"
    return(gl_results)
  }else{if(type=="fusion"){
    # Fusion penalty -----------------------------------------------------------
    # Check params
    if(is.null(params)){
      beta_start <- rep(0,NODES*(NODES-1)/2)
      b_start <- 0
      MAX_ITER <- 200;    CONV_CRIT <- 1e-04;   MAX_TIME = Inf
    }else{if(is.null(params$beta_start)) {
      beta_start <- rep(0,NODES*(NODES-1)/2)
      b_start <- 0
      MAX_ITER <- params$MAX_ITER;    CONV_CRIT <- params$CONV_CRIT;   MAX_TIME = params$MAX_TIME
    }}
    gl = logistic_fused_lasso(X, Y, D,
                              lambda1 = lambda1, lambda2 = lambda2,
                              id, verbose, beta_start = beta_start, b_start = b_start,
                              NODES = NODES, MAX_ITER = MAX_ITER, 
                              CONV_CRIT = CONV_CRIT, 
                              MAX_TIME = MAX_TIME)
    gl_results = list()
    gl_results$beta = gl$best_beta/alpha_normalization;    gl_results$b = gl$b
    
    # Train fitting
    Yfit <- X%*%gl$best_beta + gl$b
    gl_results$Yfit <- exp(Yfit)/(1+exp(Yfit))
    gl_results$train_error <- 1- sum(diag(table(sign(Yfit),sign(Y))))/length(Y)
    # Test fitting
    if(!is.null(Xtest)) {
      levels(Ytest) <- c(-1, 1)
      Yfit_test <- Xtest%*%gl_results$beta + gl_results$b
      gl_results$Yfit_test <- exp(Yfit_test)/(1+exp(Yfit_test))
      gl_results$test_error = 1-sum(diag(table(sign(Yfit_test),sign(Ytest))))/length(Ytest)
    }
    gl_results$type = "intersection"
    return(gl_results)
  }else{
    stop("The value of type should be one between \"intersection\" and \"union\"")
  }}}
}


construct_D <- function(nodes = 264) {
  require(Matrix)
  B <- array(0,dim = c(nodes,nodes))
  B[upper.tri(B)] <- 1:(nodes*(nodes-1)/2)
  B <- B + t(B)
  D <- Matrix(0, nrow=nodes*(nodes-1), ncol = nodes*(nodes-1)/2)
  for(i in 1:nodes) {
    D[((i-1)*(nodes-1)+1):(i*(nodes-1)),B[i,-i]] <- diag(nodes-1)
  }
  #d2 <- lapply(1:nodes,function(i) ((i-1)*(nodes-1)+1):(i*(nodes-1)))
  return(D)
}


construct_D_fusion <- function(nodes = 264) {
  require(Matrix)
  U  = cBind(rep(1,nodes-2),-Diagonal(nodes-2))
  W = U
  dim(W)
  for(i in 2:(nodes-3))
    W = rBind(W,cBind(array(0,dim = c(nodes-i-1,i-1)),
                      (U[1:(nodes-i-1),1:(nodes-i)])))
  W = rBind(W,array(c(rep(0,nodes-3),1,-1),dim = c(1,nodes-1)))
  
  aux = array(0,dim = c(nodes,nodes))
  edges = sum(upper.tri(aux))
  aux[upper.tri(aux)] = 1:edges
  aux = aux+t(aux)
  D = Matrix(0,nrow = nrow(W)*nodes,ncol = edges)
  D = Matrix(0,nrow = 0,ncol = edges)
  for(i in 1:nodes) { 
    aux2 = Matrix(0,nrow = nrow(W),ncol = edges)
    aux2[,aux[i,which(aux[i,]!=0)]] = W
    D = rBind(D,aux2)
    #print(i)
  }
  return(D)
}
