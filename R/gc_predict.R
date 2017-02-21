#' Predict function for graph classifier
#'
#' @param gc A trained graph classifier object
#' @param X A matrix containing rows with vectorized upper triangular adjacency matrices (column-major order)
#' @param type Indicates the type of response. class: predicted classes. prob: predicted probabilities. error: misclassification error
#' @param Ytest If type = "error", true classes to compare.
#' @return A vector containing the predicted classes.
#' @examples
#' X = matrix(rnorm(100*34453), nrow = 100)
#' Y = 2*(runif(100) > 0.5) - 1
#' gc = graphclass(X, Y = Y)
#' Xtest = matrix(rnorm(100*34453), nrow = 100)
# gc_predict(gc, Xtext) 
gc_predict <- function(gc,X, type = "class", Ytest) {
  pred <- X %*% gc$beta + gc$b
  Ypred <- sapply(c(pred), function(y) ifelse(y>0, gc$Ypos_label, gc$Yneg_label))
  if(type=="class")
    return(Ypred)
  if(type=="pred")
    return(exp(pred) / (1 + exp(pred)))
  if(type == "error")
    return(sum(Ytest!=Ypred)/length(Ypred))
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
