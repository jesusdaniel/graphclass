#' Returns a matrix from a vectorized network
#'
#' @param beta Vectorized adjacency matrix. 
#' @param  type Parameter to specify whether the vector represents an \code{undirected}
#' or \code{directed} network. Default is \code{undirected}.
#' 
#' @return Adjacency matrix for a vectorized network
#' 
#' # Obtain the adjacency matrix of a COBRE data subject
#' data(COBRE.data)
#' A <- get_matrix(COBRE.data$X.cobre[1,])
#' 
#' @encoding UTF-8
#' @importFrom Rdpack reprompt
get_matrix <- function(beta, type=c("undirected", "directed")) {
  type <- match.arg(type)
  if(type=="undirected") {
    NODES <- (1+sqrt(1+8*length(beta)))/2
    Adj_matr = array(0,dim = c(NODES, NODES))
    Adj_matr[upper.tri(Adj_matr)] <- beta
    Adj_matr <- Adj_matr + t(Adj_matr)  
  }else{if(type=="directed") {
    NODES <- (1+sqrt(1+8*length(beta)/2))/2
    Adj_matr = array(T,dim = c(NODES, NODES))
    diag(Adj_matr) = F
    Adj_matr[Adj_matr] = beta
  }else{
    stop("The value of type should be one between \"directed\" and \"undirected\"")
  }}
  return(Adj_matr)
}

#' Returns a vectorized adjacency matrix
#'
#' @param A Adjacency matrix of a network.
#' @param  type Parameter to specify whether the adjacency matrix represents an \code{undirected}
#' or \code{directed} network. Default is \code{undirected}.
#' 
#' @return Vector represeingt an adjacency matrix
#' 
#' @encoding UTF-8
#' @importFrom Rdpack reprompt
matrix_to_vec <- function(A, type=c("undirected", "directed")) {
  type <- match.arg(type)
  if(type=="undirected") {
    beta <- A[upper.tri(A)]
  }else{if(type=="directed") {
    mat <- matrix(T, ncol = ncol(A), nrow = nrow(A))
    diag(mat) <- F
    beta <- as.vector(A[mat])
  }else{
    stop("The value of type should be one between \"directed\" and \"undirected\"")
  }}
  return(beta)
}
#' Calculates the node sparsity of a vectorized adjacency matrix
#'
#' @param beta Vectorized adjacency matrix. 
#' @return Percentage of inactive nodes in the graph solution
node_sparsity <- function(beta) {
  A <- get_matrix(beta)
  return(sum(apply(A,1,function(v) sum(v!=0))==0)/ncol(A))
}


#' Constructor for D penalty matrix, use it for efficiency when running the classifier multiple times.
#'
#' @param nodes Number of nodes in the network, by default is 264 (Power parcellation).
#' @return A sparse D matrix
#' @examples
#' D = construct_D(100)
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
   pmax.fast(0, abs(x)-lambda) * sign(x)
   #.C("soft_thresholding", x = as.double(x), lambda = as.double(lambda), n = as.integer(n))$x
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
  .C("l2_node_soft_ghtesholding", 
     Db = as.double(Db), N=as.integer(n), lambda = as.double(lambda))$Db
}
gl_penalty.c = function(b, NODES) {
  gl = 0;
  .C("group_lasso_node_penalty",Db=as.double(b),N = as.integer(NODES), gl = as.double(gl))$gl
}
# l2_soft_thresholding.c = function(x, lambda) {
#   n = length(x)
#   .C("l2_soft_thresholding", x = as.double(x), lambda = as.double(lambda), n = as.integer(n))$x
# }

l1l2_proximal.c = function(beta, lambda1, lambda2, NODES) {
  .C("l1l2_node_soft_thresholding", b = as.double(beta), N=as.integer(NODES), 
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
