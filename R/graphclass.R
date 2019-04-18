# https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html
# install.packages(repos=NULL, "graphclass_1.0.tar.gz")
# R CMD Rd2pdf graphclass
#' Regularized logistic regression classifier for networks.
#' 
#' \code{graphclass} fits a regularized logistic regression to a set of network adjacency matrices with responses, and returns an 
#' object with the classifier.
#' 
#' The function  \code{graphclass} fits a regularized logistic regression to classify a set of network adjacency matrices
#' with \eqn{N} labeled nodes and corresponding responses. The classifier fits a matrix of coefficients \eqn{B\in{R}^{N\times N}},
#' in which \eqn{B_{ij}} indicates the coefficient corresponding to the edge \eqn{(i,j)}.
#' 
#'  The argument \code{type} provides options to choose the penalty function.
#' If \code{type = "intersection"} or \code{"union"}, the penalty corresponds to the node selection penalty defined as
#' \deqn{\Omega(B) = \lambda \left(\sum_{i=1}^N\sqrt{\sum_{j=1}^N B_{ij}^2} + \rho \sum_{i=1}^N\sum_{j=1}^N|B_{ij}|\right).}
#' When \code{type = "intersection"}, a symmetric restriction on  \eqn{B} is enforced, and the penalty promotes subgraph selection.
#' If \code{type = "union"}, the penalty promotes individual node selection.
#' See \insertCite{relion2017network;textual}{graphclass} for more details.
#' 
#' The value \code{type = "groups"} corresponds to a generic  group lasso penalty. The groups of edges have to be specified using the argument \code{Groups} with a list of arrays,
#' in which each element of the list corresponds to a group, and the array indicates the indexes of the variables in that group.
#' The optional argument \code{G_penalty_factors} is an array of  size equal to the number of groups, and can be used to 
#' specify different weights for each group on the penalty (for example, when groups have different sizes).
#' 
#' The optional argument \code{params} is a list that allows to control some internal parameters of the optimization algorithm. 
#' The elements \code{beta_start} and \code{b_start} are  initial values for the optimization algorithm. The value
#' of \code{beta_start} is a vector that indicates the weights of the upper triangular part of \eqn{B}, and \code{b_start}
#' is the initial value of the threshold in the logistic regression. By default, these parameters are set to zero. The elements
#' \code{MAX_ITER} and \code{CONV_CRIT} can be used to change the maximum number of iterations and the convergence criterion in the proximal algorithm
#' for fitting the node selection penalty (see \insertCite{relion2017network;textual}{graphclass}). By default, these values are set to
#' \code{MAX_ITER=300} and \code{CONV_CRIT = 1e-5}.
#' 
#' @rdname graphclass
#' @export
#' 
#' @references
#' \insertRef{relion2017network}{graphclass}
#'
#' @seealso \code{\link{plot.graphclass}}, \code{\link{predict.graphclass}}
#' 
#' @param X A matrix with the training samples, in wich each row represents the vectorized (by column order) upper triangular part of a network adjacency matrix.
#' @param Adj_list A training list of of symmetric adjacency matrices with zeros in the diagonal
#' @param Y A vector containing the class labels of the training samples (only 2 classes are supported for now).
#' @param Xtest Optional argument for providing a matrix containing the test samples, with each row representing an upper-triangular vectorized adjacency matrix.
#' @param Ytest Optional argument containing the labels of test samples.
#' @param type  Type of penalty function. Default is \code{"intersection".}
#' See details.
#' @param lambda penalty parameter \eqn{lambda}, by default is set to 0.
#' @param rho penalty parameter \eqn{rho} controlling sparsity, by default is set to 0.
#' @param gamma ridge parameter (for numerical purposes). Default is \code{gamma = 1e-5}.
#' @param params A list containing internal parameters for the optimization algorithm. See details.
#' @param verbose whether output is printed
#' @param D matrix \eqn{D} used by the penalty to define the groups. This optional argument can be used to pass a precomputed matrix \code{D}, which can be time saving if the method is fitted multiple times. See the function \code{construct_D}.
#' @param Groups list of lists, where each list correspond to a grouping and each sublist to sets of indexes in X. Each sublist should be a non-overlapping group.
#' @param G_penalty_factors For type "groups", each group is penalized by this factor. Should sum to 1.
#' 
#' @return An object containing the trained graph classifier.
#' \item{beta}{Edge coefficients vector of the regularized logistic regression solution.}
#' \item{b}{Intercept value.}
#' \item{Yfit}{Fitted logistic regression probabilities in the train data.}
#' \item{Ypred}{Predicted class for the test samples (if available).}
#' \item{train_error}{Percentage of train samples that are misclassified.}
#' \item{test_error}{Percentage of test samples that are misclassified (if available).}
#' 
#' @examples
#' 
#' # Load COBRE data
#' data(COBRE.data)
#' X <- COBRE.data$X.cobre
#' Y <- COBRE.data$Y.cobre
#' 
#' # An example of the subgraph selection penalty
#' gc = graphclass(X = X, Y = factor(Y), type = "intersection",
#'                lambda = 1e-4, rho = 1, gamma = 1e-5)
#' plot(gc)
#' 
#' 
#' # 5-fold cross validation
#' fold_index <- (1:length(Y) %% 5) + 1
#' 
#' gclist <- list()
#' for(fold in 1:5) {
#'     foldout <- which(fold_index == fold) 
#'     gclist[[fold]] <- graphclass(X = X[-foldout,], Y = factor(Y[-foldout]),
#'                      Xtest = X[foldout,], Ytest = factor(Y[foldout]),
#'                      type = "intersection",
#'                      lambda = 1e-4, rho = 1, gamma = 1e-5,
#'                      D = D263)
#' }
#' # test error on each fold
#' lapply(gclist, function(gc) gc$test_error)
#' 
#' @encoding UTF-8
#' @importFrom Rdpack reprompt
#' @useDynLib Pmax2
graphclass <- function(X = NULL, Y = NULL, 
                       type = c("intersection", "union", "groups", "fusion"),...) {
  UseMethod("graphclass")
}


#' @rdname graphclass
#' @export
graphclass.default <- function(X = NULL, Y = NULL,
                       Xtest = NULL, Ytest = NULL,
                       Adj_list = NULL, 
                       type = c("intersection", "union", "groups", "fusion"),
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
  type <- match.arg(type)
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
        require(Matrix)
        D_list <- lapply(Groups, function(G) {
          D <- Matrix(0, nrow = ncol(X), ncol = ncol(X))
          Matrix::diag(D)[G] <- 1
          D
          })
        if(is.null(G_penalty_factors)) {
          G_penalty_factors <- rep(1/length(D_list), length(D_list))
        }
      }else{
        stop("The value of type should be one 
             between \"intersection\", \"union\" or \"groups\".")
          }
        }
      }
    }
  # Check which method to use
  # Intersection penalty -----------------------------------------------
  if(type=="intersection") {
    beta_start <- rep(0,NODES*(NODES-1)/2)
    b_start <- 0
    MAX_ITER <- 300;    CONV_CRIT <- 1e-05;   MAX_TIME = Inf
    # Check params
    if(!is.null(params)){
      if(!is.null(params$beta_start))
        beta_start <- params$beta_start*alpha_normalization
      if(!is.null(params$b_start))
        b_start <- params$b_start
      if(!is.null(params$MAX_ITER))
        MAX_ITER <- params$MAX_ITER
      if(!is.null(params$CONV_CRIT))
        CONV_CRIT <- params$CONV_CRIT
        
    }
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
    beta_start <- rep(0,NODES*(NODES-1)/2)
    b_start <- 0
    MAX_ITER <- 300;    CONV_CRIT <- 1e-05;   MAX_TIME = Inf
    # Check params
    if(!is.null(params)){
      if(!is.null(params$beta_start))
        beta_start <- params$beta_start*alpha_normalization
      if(!is.null(params$b_start))
        b_start <- params$b_start
      if(!is.null(params$MAX_ITER))
        MAX_ITER <- params$MAX_ITER
      if(!is.null(params$CONV_CRIT))
        CONV_CRIT <- params$CONV_CRIT
    }
    print(CONV_CRIT)
    print(MAX_ITER)
    gl = logistic_group_lasso_ridge_groups(X, Y, Groups,
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
    beta_start <- rep(0,NODES*(NODES-1))
    b_start <- 0
    MAX_ITER <- 300;    CONV_CRIT <- 1e-05;   MAX_TIME = Inf
    # Check params
    if(!is.null(params)){
      if(!is.null(params$beta_start))
        beta_start <- params$beta_start*alpha_normalization
      if(!is.null(params$b_start))
        b_start <- params$b_start
      if(!is.null(params$MAX_ITER))
        MAX_ITER <- params$MAX_ITER
      if(!is.null(params$CONV_CRIT))
        CONV_CRIT <- params$CONV_CRIT
      
    }
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




.onLoad <- function(lib, pkg){
  Rdpack::Rdpack_bibstyles(package = pkg, authors = "LongNames")
  invisible(NULL)
}

