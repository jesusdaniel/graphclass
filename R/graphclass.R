# https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html
# install.packages(repos=NULL, "graphclass_1.0.tar.gz")
# R CMD Rd2pdf --pdf --title=graphclass -o graphclass.pdf man/*.Rd
#' Train a graph classifier using regularized logistic regression.
#'
#' @param X A matrix with the training sample, in wich each row represents a vectorized (by column order) upper triangular part of a network.
#' @param Adj_list A list of of symmetric matrices with 0 diagonal for training the classifier
#' @param Y A vector containing the class labels of the training sample (for now only 2 classes are supported).
#' @param Xtest A optional test matrix.
#' @param Ytest Labels of test set.
#' @param type should be either "intersection", "union", "fusion" or "groups". only "intersection" is currently supported.
#' @param lambda penalty parameter \eqn{lambda}, by default is set to 0.
#' @param rho penalty parameter \eqn{rho} controlling sparsity, by default is set to 0.
#' @param gamma ridge parameter (for numerical purposes).
#' @param params A list containing threshold parameters for the algorithm (see details)
#' @param verbose whether output is printed
#' @param D matrix \eqn{D} of the penalty; precomputing it can save time.
#' @param Groups list of lists, where each list correspond to a grouping and each sublist to sets of indexes in X. Each sublist should be a non-overlapping group.
#' @param G_penalty_factors For type "groups", each group is penalized by this factor. Should sum to 1.
#' @return An object containing the trained graph classifier.
#' @examples
#' X = matrix(rnorm(100*34453), nrow = 100)
#' Y = 2*(runif(100) > 0.5) - 1
#' gc = graphclass(X, Y = factor(Y))
#' gc$train_error
# Params: 
# - beta_start, b_start, MAX_ITER, CONV_CRIT, MAX_TIME
graphclass <- function(X = NULL, Y = NULL,...) {
  UseMethod("graphclass")
}


#' @rdname graphclass
#' @export
graphclass.default <- function(X = NULL, Y = NULL,
                       Xtest = NULL, Ytest = NULL,
                       Adj_list = NULL, 
                       type = "intersection",
                       lambda = 0, rho = 0,
                       gamma = 1e-5, 
                       params = NULL, id = "", verbose = F, D = NULL,
                       Groups = NULL,
                       G_penalty_factors = NULL,
                       ...) {
  gl_results = list()
  gl_results$lambda <- lambda
  gl_results$rho <- rho
  gl_results$gamma <- gamma
  # Dependent data
  if(!is.null(X)) {
    # ncols = N*(N-1)/2
    NODES <- (1+sqrt(1+8*ncol(X)))/2
  }else{
    NODES <- ncol(Adj_list[[1]])
    X <- t(sapply(Adj_list, function(A) A[upper.tri(A)]))
  }
  Yoriginal <- Y
  Y <- as.numeric(Y)
  Y <- (Y-min(Y))/(max(Y)-min(Y))
  Y <- 2*Y-1
  gl_results$Ypos_label <- unique(Yoriginal[Y==1])[1]
  gl_results$Yneg_label <- unique(Yoriginal[Y==-1])[1]
  table(Y, Yoriginal)
  # Normalize X for numerical stability
  alpha_normalization <- max(apply(X, 2, function(v) sd(v)))
  X <- X/alpha_normalization
  lambda1 <- lambda*rho/alpha_normalization
  lambda2 <- lambda/alpha_normalization
  gamma <- gamma/alpha_normalization
  # Create D
  if(is.null(D)) {
    if(type=="intersection"| type=="union") {
      D <- construct_D(NODES)  
    }else{
      if(type=="fusion")
        D <- construct_D_fusion(NODES)
      else{if(type=="groups") {
        D_list <- lapply(Groups, function(G) {
          D <- Matrix(0, nrow = length(G), ncol = ncol(X))
          for(i in 1:length(G)) {
            D[i, G[[i]]] <- 1
          }
          D
        })
        if(is.null(G_penalty_factors)) {
          G_penalty_factors <- rep(1/length(D_list), length(D_list))
        }
      }else{
        stop("The value of type should be one 
             between \"intersection\", \"union\" or \"fusion\".")
          }
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
    }else{
      beta_start <- params$beta_start*alpha_normalization
      b_start <- params$b_start
      MAX_ITER <- params$MAX_ITER;    CONV_CRIT <- params$CONV_CRIT;   MAX_TIME = params$MAX_TIME
    }}
    gl = logistic_group_lasso_ridge(X, Y, D,
                                    lambda1 = lambda1, 
                                    lambda2 = lambda2,
                                    gamma = gamma,
                                    id, verbose, 
                                    beta_start = beta_start, 
                                    b_start = b_start,
                                    NODES = NODES, 
                                    MAX_ITER = MAX_ITER, 
                                    CONV_CRIT = CONV_CRIT, 
                                    MAX_TIME = MAX_TIME)
    gl_results$beta = gl$best_beta/alpha_normalization
    gl_results$b = gl$best_b
    
    # Train fitting
    Yfit <- alpha_normalization*(X%*%gl_results$beta) + gl_results$b
    gl_results$Yfit <- exp(Yfit)/(1+exp(Yfit))
    gl_results$train_error <- 1- sum(diag(table(sign(Yfit),sign(Y))))/length(Y)
    class(gl_results) <- "graphclass"
    # Test fitting
    if(!is.null(Xtest)) {
      predict.graphclass(gl_results, Xtest)
      Yfit_test <- Xtest%*%gl_results$beta + gl_results$b
      Ypred <- 1*(Yfit_test>0) -1*(Yfit_test<=0)
      gl_results$Ypred <- predict.graphclass(gl_results, Xtest)
      gl_results$Yfit_test <- predict.graphclass(gl_results, Xtest, type = "prob")
      if(!is.null(Ytest)) {
        gl_results$test_error = predict.graphclass(gl_results, Xtest, type = "error", Ytest = Ytest)
      } 
    }
    gl_results$type = "intersection"
    gl_results$active_nodes = active_nodes(gl_results$beta, NODES = NODES)
    gl_results$subgraph_active_edges_rate = active_edges_rate(gl_results$beta, NODES = NODES)
    gl_results$nonzeros_percentage = sum(gl_results$beta!=0)/length(gl_results$beta)
    return(gl_results)
    
  }else{if(type=="groups") {
    # Check params
    if(is.null(params)){
      beta_start <- rep(0,NODES*(NODES-1)/2)
      b_start <- 0
      MAX_ITER <- 300;    CONV_CRIT <- 1e-05;   MAX_TIME = Inf
    }else{if(is.null(params$beta_start)) {
      beta_start <- rep(0,NODES*(NODES-1)/2)
      b_start <- 0
      MAX_ITER <- params$MAX_ITER;    CONV_CRIT <- params$CONV_CRIT;   MAX_TIME = params$MAX_TIME
    }else{
      beta_start <- params$beta_start*alpha_normalization
      b_start <- params$b_start
      MAX_ITER <- params$MAX_ITER;    CONV_CRIT <- params$CONV_CRIT;   MAX_TIME = params$MAX_TIME
    }}
    gl = logistic_group_lasso_ridge_groups(X, Y, D_list,
                                    lambda1 = lambda1, 
                                    lambda2 = lambda2,
                                    G_penalty_factors = G_penalty_factors,
                                    gamma = gamma,
                                    id, verbose, 
                                    beta_start = beta_start, 
                                    b_start = b_start,
                                    NODES = NODES, 
                                    MAX_ITER = MAX_ITER, 
                                    CONV_CRIT = CONV_CRIT, 
                                    MAX_TIME = MAX_TIME)
    gl_results$beta = gl$best_beta/alpha_normalization
    gl_results$b = gl$best_b
    
    # Train fitting
    Yfit <- alpha_normalization*(X%*%gl_results$beta) + gl_results$b
    gl_results$Yfit <- exp(Yfit)/(1+exp(Yfit))
    gl_results$train_error <- 1- sum(diag(table(sign(Yfit),sign(Y))))/length(Y)
    class(gl_results) <- "graphclass"
    # Test fitting
    if(!is.null(Xtest)) {
      predict.graphclass(gl_results, Xtest)
      Yfit_test <- Xtest%*%gl_results$beta + gl_results$b
      Ypred <- 1*(Yfit_test>0) -1*(Yfit_test<=0)
      gl_results$Ypred <- predict.graphclass(gl_results, Xtest)
      gl_results$Yfit_test <- predict.graphclass(gl_results, Xtest, type = "prob")
      if(!is.null(Ytest)) {
        gl_results$test_error = predict.graphclass(gl_results, Xtest, type = "error", Ytest = Ytest)
      } 
    }
    gl_results$type = "intersection"
    gl_results$active_nodes = active_nodes(gl_results$beta, NODES = NODES)
    gl_results$subgraph_active_edges_rate = active_edges_rate(gl_results$beta, NODES = NODES)
    gl_results$nonzeros_percentage = sum(gl_results$beta!=0)/length(gl_results$beta)
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
    }else{
      beta_start <- params$beta_start*alpha_normalization
      b_start <- params$b_start
      MAX_ITER <- params$MAX_ITER;    CONV_CRIT <- params$CONV_CRIT;   MAX_TIME = params$MAX_TIME
    }}
    gl = logistic_union_group_lasso_ridge(X, Y, D,
                                          lambda1 = lambda1, 
                                          lambda2 = lambda2, 
                                          gamma = gamma,
                                          id, verbose, 
                                          beta_start = beta_start, 
                                          b_start = b_start,
                                          NODES = NODES, 
                                          MAX_ITER = MAX_ITER, 
                                          CONV_CRIT = CONV_CRIT, 
                                          MAX_TIME = MAX_TIME)
    gl_results$beta = gl$best_beta/alpha_normalization
    gl_results$b = gl$best_b
    # Train fitting
    Yfit <- alpha_normalization*(X%*%(crossprod(D,gl_results$beta ))) + gl$b
    gl_results$Yfit <- exp(Yfit)/(1+exp(Yfit))
    gl_results$train_error <- 1-sum(diag(table(as.vector(sign(gl_results$Yfit-0.5)),sign(Y))))/length(Y)
    # Test fitting
    if(!is.null(Xtest)) {
      Yfit_test <- Xtest%*%crossprod(D,gl_results$beta) + gl_results$b
      gl_results$Yfit_test <- exp(Yfit_test)/(1+exp(Yfit_test))
      if(is.factor(Yoriginal)) {
        gl_results$Ypred <- factor(Ypred)
        levels(gl_results$Ypred) <- levels(Yoriginal)
      }
      if(!is.null(Ytest)) {
        levels(Ytest) <- c(-1, 1)
        gl_results$test_error = 1-sum(diag(table(sign(Yfit_test),sign(Ytest))))/length(Ytest)
      }
    }
    gl_results$active_nodes = active_nodes_union(gl_results$beta, NODES = NODES)
    gl_results$subgraph_active_edges_rate = NULL
    gl_results$nonzeros_percentage = sum(gl_results$beta!=0)/length(gl_results$beta)
    gl_results$type = "union"
    class(gl_results) <- "graphclass"
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
    }else{
      beta_start <- params$beta_start*alpha_normalization
      b_start <- params$b_start
      MAX_ITER <- params$MAX_ITER;    CONV_CRIT <- params$CONV_CRIT;   MAX_TIME = params$MAX_TIME
    }}
    print(lambda1)
    print(lambda2)
    gl = logistic_fused_lasso(X, Y, D,
                              lambda1 = lambda1, lambda2 = lambda2,
                              id, verbose, beta_start = beta_start, b_start = b_start,
                              NODES = NODES, MAX_ITER = MAX_ITER, 
                              CONV_CRIT = CONV_CRIT, 
                              MAX_TIME = MAX_TIME)
    gl_results$beta = gl$best_beta/alpha_normalization
    gl_results$b = gl$best_b
    
    # Train fitting
    Yfit <- alpha_normalization*(X%*%gl$best_beta) + gl$b
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
    class(gl_results) <- "graphclass"
    return(gl_results)
  }else{
    stop("The value of type should be one between \"intersection\" and \"union\"")
  }
  }
  }
  }
}




#' Predict function for graph classifier.
#'
#' @rdname graphclass
#' @export
#'
#' @param object trained graphclass object
#' @param newdata matrix of observations to predict. Each row corresponds to a new observation.
#' @param type type of response. class: predicted classes. prob: predicted probabilities. error: misclassification error
#' @param Ytest if type = "error", true classes to compare.
#' @return A vector containing the predicted classes.
#' @examples
#' X = matrix(rnorm(100*34453), nrow = 100)
#' Y = 2*(runif(100) > 0.5) - 1
#' gc = graphclass(X, Y = factor(Y))
#' Xtest = matrix(rnorm(100*34453), nrow = 100)
#' predictions = predict(gc, Xtest)
predict.graphclass <- function(object, newdata, type = "class", Ytest, ...) {
  if (!inherits(object, "graphclass"))  {
    stop("Object not of class 'graphclass'")
  }
  pred <- newdata %*% object$beta + object$b
  Ypred <- sapply(c(pred), function(y) if(y>0) { object$Ypos_label}else{ object$Yneg_label})
  if(type=="class")
    return(Ypred)
  if(type=="pred")
    return(exp(pred) / (1 + exp(pred)))
  if(type == "error")
    return(sum(Ytest!=Ypred)/length(Ypred))
}


#' Plot function for graph classifier.
#'
#' Plots the adjacency matrix of the coefficients network
#' 
#' @rdname graphclass
#' @export
#'
#' @param object trained graphclass object
#' @examples
#' X = matrix(rnorm(100*34453), nrow = 100)
#' Y = 2*(runif(100) > 0.5) - 1
#' gc = graphclass(X, Y = factor(Y))
#' plot(gc)
plot.graphclass <- function(object, ...) {
  if (!inherits(object, "graphclass"))  {
    stop("Object not of class 'graphclass'")
  }
  plot_adjmatrix(object$beta)
}

