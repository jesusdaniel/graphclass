#' Converts a vector into an adjacency matrix
#' 
#' Given a vector that encodes an adjacency matrix, returns the matrix representation.
#' 
#' @export
#' 
#' @param beta Vectorized adjacency matrix. If the network is undirected,
#' the vector is assumed to represent the upper triangular part of the adjacency 
#' matrix in column major order. For directed networks, the vector contains the
#' entries ordered by column, and excluding the diagonal elements if \code{selfloops = FALSE}.
#' @param  type Specifies whether the vector represents an \code{undirected}
#' or \code{directed} network. Default is \code{undirected}.
#' @param selfloops This parameter indicates whether self loops are included in the
#' entries of a vector. Default is \code{FALSE}.
#' 
#' @return Adjacency matrix. 
#' 
#' @examples
#' 
#' # Obtain the adjacency matrix of a COBRE data subject
#' data(COBRE.data)
#' A <- get_matrix(COBRE.data$X.cobre[1,], type = "undirected")
#' plot_adjmatrix(A)
#' 
#' @encoding UTF-8
#' @importFrom Rdpack reprompt
get_matrix <- function(beta, type=c("undirected", "directed"), selfloops = FALSE) {
  type <- match.arg(type)
  if(type=="undirected" & !selfloops) {
    NODES <- (1+sqrt(1+8*length(beta)))/2
    Adj_matr = array(0,dim = c(NODES, NODES))
    Adj_matr[upper.tri(Adj_matr)] <- beta
    Adj_matr <- Adj_matr + t(Adj_matr)  
  }else{if(type=="directed" & !selfloops) {
    NODES <- (1+sqrt(1+8*length(beta)/2))/2
    Adj_matr = array(T,dim = c(NODES, NODES))
    diag(Adj_matr) = F
    Adj_matr[Adj_matr] = beta
  }else{if(type=="undirected" & selfloops){
    NODES <- (-1+sqrt(1+8*length(beta)))/2
    Adj_matr = array(0,dim = c(NODES, NODES))
    Adj_matr[upper.tri(Adj_matr, diag = TRUE)] <- beta
    Adj_matr <- Adj_matr + t(Adj_matr)  
    diag(Adj_matr) = diag(Adj_matr)/2
  }else{if(type=="directed" & selfloops){
    NODES <- sqrt(length(beta))
    Adj_matr = array(beta,dim = c(NODES, NODES))
  }else{
    stop("The value of type should be one between \"directed\" and \"undirected\"")
  }}}}
  return(Adj_matr)
}

#' Convert matrix to vector.
#' 
#' The function encodes an adjacency matrix into a vector.
#'
#' @param A Adjacency matrix of a network.
#' @param  type Parameter to specify whether the adjacency matrix represents an \code{undirected}
#' or \code{directed} network. Default is \code{undirected}.
#' @param selfloops This parameter indicates whether self loops are included in the
#' entries of a vector. Default is \code{FALSE}.
#' 
#' @return A vector containing the upper triangular part of an
#'  adjacency matrix (if the graph is undirected) or adjacency entry
#'  values ordered by column and excluding the diagonal entries (if the graph
#'  is directed).
#' 
#' @export
#' 
#' @encoding UTF-8
#' @importFrom Rdpack reprompt
matrix_to_vec <- function(A, type=c("undirected", "directed"), selfloops = FALSE) {
  type <- match.arg(type)
  if(type=="undirected") {
    beta <- A[upper.tri(A, diag = selfloops)]
  }else{if(type=="directed") {
    mat <- matrix(T, ncol = ncol(A), nrow = nrow(A))
    diag(mat) <- selfloops
    beta <- as.vector(A[mat])
  }else{
    stop("The value of type should be one between \"directed\" and \"undirected\"")
  }}
  return(beta)
}

#' Percentage of inactive nodes in a network
#' 
#' Calculates the node sparsity of a vectorized adjacency matrix, defined
#' as the average number of nodes with no edges.
#' 
#' @param beta Vectorized adjacency matrix. 
#' @param  type Specifies whether the vector represents an \code{undirected}
#' or \code{directed} network. Default is \code{undirected}.
#' 
#' @return Percentage of inactive nodes in the graph
#' 
#' @examples
#' A <- matrix(0, ncol = 4, nrow = 4)
#' A[2, 1] <- 1
#' A[1, 2] <- 1
#' A[2, 3] <- 1
#' A[3, 2] <- 1
#' vec <- matrix_to_vec(A)
#' node_sparsity(vec)
#' 
#' @export
#' 
#' @encoding UTF-8
#' @importFrom Rdpack reprompt
node_sparsity <- function(beta, type=c("undirected", "directed")) {
  type <- match.arg(type)
  A <- get_matrix(beta, type = type)
  return(sum(apply(A,1,function(v) sum(v!=0))==0)/ncol(A))
}


#' Constructor for penalty matrix
#' 
#' This function constructs an auxiliary matrix to compute the node penalty in the regularized
#' classifier, which is internally used by the optimization algorithm. This matrix can be passed
#' to the function \code{\link{graphclass}} in order to avoid creating it every time this
#' function is called.
#'  Use it for efficiency when running the classifier multiple times (for example, 
#'  in a cross-validation routine).
#'
#' @param nodes Number of nodes in the network, by default is 264 (Power parcellation).
#' @return A sparse matrix
#' 
#' @examples
#' D263 = construct_D(263)
#' 
#' # Load COBRE data
#' data(COBRE.data)
#' X <- COBRE.data$X.cobre
#' Y <- COBRE.data$Y.cobre
#' 
#' # 5-fold cross validation of the subgraph selection penalty
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
#' @export
#' 
#' @encoding UTF-8
#' @importFrom Rdpack reprompt
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

library(compiler)
norm_s = function(x) sqrt(sum(x^2))
norm_s = cmpfun(norm_s, options = list(optimize = 3))


# 
# Pmax2.C = function(x) {
#   n = length(x)
#   .C("pmax2", x = as.double(x), n = as.integer(n))$x
# } 

# soft_thresholding.c = function(x, lambda) {
#   n = length(x)
#   #pmax(abs(x)-lambda, 0)*sign(x)
#   .C("soft_thresholding", x = as.double(x), lambda = as.double(lambda), n = as.integer(n))$x
# }

#http://r.789695.n4.nabble.com/pmin-pmax-slower-than-necessary-for-common-cases-td909143.html
pmax.fast <- function(scalar, vector) 
{ 
  ## Purpose: FAST substitute for pmax(s, v) when length(s) == 1 
  ## Author: Martin Maechler, Date: 21 Aug 1992 (for S) 
  vector[ scalar > vector ] <- scalar 
  vector 
} 

soft_thresholding.c = function(x, lambda) {
   n = length(x)
   #pmax.fast(0, abs(x)-lambda) * sign(x)
   .C("soft_thresholding", #PACKAGE = "graphclass",
      x = as.double(x), lambda = as.double(lambda), n = as.integer(n))$x
}

soft_thresholdingl2 = function(x, lambda) {
  if(norm_s(x)<lambda) {
    return(rep(0,length(x)))
  }
  else
    return(x*(1-lambda/norm_s(x)))
}

l2_Db_soft_thresholding.c = function(Db,lambda, NODES) {  
  n = NODES
  .C("l2_node_soft_ghtesholding", #PACKAGE = "graphclass",
     Db = as.double(Db), N=as.integer(n), lambda = as.double(lambda))$Db
}
gl_penalty.c = function(b, NODES) {
  gl = 0;
  .C("group_lasso_node_penalty", #PACKAGE = "graphclass",
     Db=as.double(b),N = as.integer(NODES), gl = as.double(gl))$gl
}
# l2_soft_thresholding.c = function(x, lambda) {
#   n = length(x)
#   .C("l2_soft_thresholding", x = as.double(x), lambda = as.double(lambda), n = as.integer(n))$x
# }

l1l2_proximal.c = function(beta, lambda1, lambda2, NODES) {
  .C("l1l2_node_soft_thresholding", #PACKAGE = "graphclass",
     b = as.double(beta), N=as.integer(NODES), 
     lambda1 = as.double(lambda1), lambda2 = as.double(lambda2))$b
}

## Fused lasso
invert_D_full <- function(D,N = 264) {
  require(Matrix)
  n = ncol(D)
  whichEdges = array(T,dim = c(N,N))
  diag(whichEdges) = F
  whichEdges[which(upper.tri(whichEdges) & whichEdges,arr.ind = T)] = 1:n
  whichEdges[lower.tri(whichEdges)]=0
  whichEdges = whichEdges + t(whichEdges)
  U = Matrix(apply(whichEdges,1,function(x) {u = rep(0,n);u[x] = 1;return(u)}))
  
  multiply_by_inverse <- function(x,rho){
    kappa = 1/(1+rho*(2*N-1))
    eta = 1/(1-rho*kappa*(N-2))
    res = kappa*(x+rho*kappa*(eta*U%*%(crossprod(U,x)) + 4*eta^2*rho*kappa/(1-rho*kappa*N*eta)*sum(x)))    
    return(res)
  }
  return(list(multiply_by_inverse = multiply_by_inverse))
}


## 
l2_groups_softhtresholding <- function(BetaVrho, D_list, lambda, G_penalty_factors) {
  require(Matrix)
  for(g in 1:length(D_list)) {
    BetaVrho[,g] = max(1- lambda *  G_penalty_factors[g] /sqrt(sum(BetaVrho[,g]^2)),0) * BetaVrho[,g]
  }
  BetaVrho
}
