# Simulate random networks
sbm <- function(p,q, K, nk) {
  Psi = (p-q)*diag(K)  + q
  A = kronecker(Psi,matrix(rep(1,(nk)^2),nrow = nk))
  Aobs = apply(A,MARGIN = c(1,2),function(u) rbinom(1,1,prob = u))
  Aobs[lower.tri(Aobs)] = 0
  Aobs =  Aobs + t(Aobs)
  diag(Aobs) = 0
  Degree = apply(Aobs,1, sum)
  L = diag(Degree) - Aobs
  return(list(Aobs= Aobs, membership = kronecker(1:K,rep(1,nk)),Laplacian = L))
}

sbm_q1 <- function(p,q, K, nk, q1) {
  Psi = (p-q)*diag(K)  + q
  Psi[2,1] = q1
  Psi[1,2] = q1
  A = kronecker(Psi,matrix(rep(1,(nk)^2),nrow = nk))
  Aobs = apply(A,MARGIN = c(1,2),function(u) rbinom(1,1,prob = u))  
  Aobs[lower.tri(Aobs)] = 0
  Aobs =  Aobs + t(Aobs)
  diag(Aobs) = 0
  Degree = apply(Aobs,1, sum)
  L = diag(Degree) - Aobs
  return(list(Aobs= Aobs, membership = kronecker(1:K,rep(1,nk)),Laplacian = L))
}

sbm_q2 <- function(p,q, K, nk, q1) {
  Psi = (p-q)*diag(K)  + q
  Psi[3,4] = q1
  Psi[4,3] = q1
  A = kronecker(Psi,matrix(rep(1,(nk)^2),nrow = nk))
  Aobs = apply(A,MARGIN = c(1,2),function(u) rbinom(1,1,prob = u))  
  Aobs[lower.tri(Aobs)] = 0
  Aobs =  Aobs + t(Aobs)
  diag(Aobs) = 0
  Degree = apply(Aobs,1, sum)
  L = diag(Degree) - Aobs
  return(list(Aobs= Aobs, membership = kronecker(1:K,rep(1,nk)),Laplacian = L))
}


signal_subgraph <- function(m,s, p, q, N, n = 70) {
  signal_vertex = sample(x = 1:n,size = m,replace = F)
  selected_edges = sample(1:(m*(2*n-m-1)/2),s,replace = F)
  signal_edges = lapply(1:m, function(x) c())
  for(i in 1:m) {
    u = selected_edges[which(selected_edges<= i*n-i*(i+1)/2 & selected_edges >(i-1)*n-(i-1)*(i)/2)]
    signal_edges[[i]] = sapply(u, function(x) x+sum(x>signal_vertex[1:i]))
  }
  prob_matrix = array(p,dim = c(n,n))
  for(i in 1:m){
    prob_matrix[signal_vertex[i], signal_edges[[i]]] = q
    prob_matrix[signal_edges[[i]], signal_vertex[i]] = q
  }
  Alist = lapply(1:N,function(x) {
    Aobs = apply(prob_matrix,MARGIN = c(1,2),function(u) rbinom(1,1,prob = u))
    Aobs[lower.tri(Aobs)] = 0
    Aobs =  Aobs + t(Aobs)
    return(Aobs)
    diag(Aobs) = 0
  })
  return(list(A_sample = Alist, signal_vertex = signal_vertex, signal_edges = signal_edges))
}

model_signal_subgraph <- function(m,s,p,q,n=70) {
  signal_vertex = sample(x = 1:n,size = m,replace = F)
  selected_edges = sample(1:(m*(2*n-m-1)/2),s,replace = F)
  signal_edges = lapply(1:m, function(x) c())
  for(i in 1:m) {
    u = selected_edges[which(selected_edges<= i*n-i*(i+1)/2 & selected_edges >(i-1)*n-(i-1)*(i)/2)]
    signal_edges[[i]] = sapply(u, function(x) x+sum(x>signal_vertex[1:i]))
  }
  prob_matrix = array(p,dim = c(n,n))
  for(i in 1:m){
    prob_matrix[signal_vertex[i], signal_edges[[i]]] = q
    prob_matrix[signal_edges[[i]], signal_vertex[i]] = q
  }
  return(list(signal_vertex = signal_vertex, signal_edges = signal_edges, prob_matrix= prob_matrix))
}
sample_probMatrix <- function(prob_matrix, n = 70, N) {
  Alist = lapply(1:N,function(x) {
    Aobs = apply(prob_matrix,MARGIN = c(1,2),function(u) rbinom(1,1,prob = u))
    Aobs[lower.tri(Aobs)] = 0
    Aobs =  Aobs + t(Aobs)
    return(Aobs)
    diag(Aobs) = 0
  })
  return(list(A_sample = Alist))
}