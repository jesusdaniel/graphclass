active_nodes <- function(beta, NODES = 263) {
  Adj_matr = array(0,dim = c(NODES, NODES))
  Adj_matr[upper.tri(Adj_matr)] <- beta
  Adj_matr <- Adj_matr + t(Adj_matr)
  return( sum(apply(Adj_matr,1,sum)!=0))
}

active_edges_rate <- function(beta, NODES = 263) {
  act_nodes = active_nodes(beta, NODES)
  max_edgs = act_nodes*(act_nodes-1)/2
  return(sum(beta!=0)/max_edgs)
}

active_nodes_union <- function(beta, NODES = 263) {
  Adj_matr = array(T,dim = c(NODES, NODES))
  diag(Adj_matr)=F
  Adj_matr[Adj_matr] <- as.double(beta)
  return( sum(apply(Adj_matr,2,sum)!=0))
}

classification_error <- function(X,Y,beta, b) {
  if(is.null(X))
    return(Inf)
  Yfit <- X%*%beta + b
  return(1- sum(diag(table(sign(Yfit),sign(Y))))/length(Y))
}
