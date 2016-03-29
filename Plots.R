plot_adjmatrix <- function(beta, type="intersection") {
  #browser()
  require(lattice)
  if(type=="intersection") {
    NODES <- (1+sqrt(1+8*length(beta)))/2
    Adj_matr = array(0,dim = c(NODES, NODES))
    Adj_matr[upper.tri(Adj_matr)] <- beta
    Adj_matr <- Adj_matr + t(Adj_matr)  
  }else{if(type=="union") {
    NODES <- (1+sqrt(1+8*length(beta)/2))/2
    #D <- construct_D(NODES)
    #v = crossprod(D,beta)
    Adj_matr = array(0,dim = c(NODES, NODES))
    delta <- row(Adj_matr) - col(Adj_matr)
    #Adj_matr[upper.tri(Adj_matr)] = as.vector(v)
    Adj_matr[delta!=0] = beta
    #Adj_matr = Adj_matr+t(Adj_matr)
  }else{
    stop("The value of type should be one between \"intersection\" and \"union\"")
  }}
  cuts = 100
  df = (Adj_matr)
  df <- df[,seq(from=ncol(df),to=1,by=-1)] #reversing the columns
  levelplot(df, at = c(seq(min(Adj_matr),-1e-10,length.out =cuts),0,seq(1e-10,max(Adj_matr),length.out =cuts)),
            xlab = "Nodes", ylab = "Nodes",
            col.regions = c(rgb((cuts:0)/cuts, green = 0, blue = 1,red = 0),rgb(red = 0,green = 0, blue = 0),
                            rgb((0:cuts)/cuts, green = 0, blue = 0,red = 1)),
            scales = list(x = list(labels=seq(0,ncol(df),10)),
                          y = list(labels = rev(seq(0,ncol(df),10)))))
}

plot_adj_community <- function(beta, communities, type="intersection") {
  require(lattice)
  if(type=="intersection") {
    NODES <- (1+sqrt(1+8*length(beta)))/2
    Adj_matr = array(0,dim = c(NODES, NODES))
    Adj_matr[upper.tri(Adj_matr)] <- beta
    Adj_matr <- Adj_matr + t(Adj_matr)  
  }else{if(type=="union") {
    NODES <- (1+sqrt(1+8*length(beta)/2))/2
    D <- construct_D(NODES)
    v = crossprod(D,beta)
    Adj_matr = array(0,dim = c(NODES, NODES))
    Adj_matr[upper.tri(Adj_matr)] = as.vector(v)
    Adj_matr = Adj_matr+t(Adj_matr)
  }else{
    stop("The value of type should be one between \"intersection\" and \"union\"")
  }}
  
  df = data.frame(Adj_matr)
  df <- t(df)
  df <- df[,seq(from=ncol(df),to=1,by=-1)] #reversing the columns
  
  cuts = 100
  levelplot(Adj_matr[comms_order,comms_order], at = c(seq(min(Adj_matr),-1e-10,length.out =cuts),0,seq(1e-10,max(Adj_matr),length.out =cuts)),
            xlab = "Nodes", ylab = "Nodes",
            col.regions = c(rgb((cuts:0)/cuts, green = 0, blue = 1,red = 0),rgb(red = 0,green = 0, blue = 0),
                            rgb((0:cuts)/cuts, green = 0, blue = 0,red = 1)))
  comms_order = unlist(sapply(unique(communities), function(x) which(communities==x)))
  length_comms = sapply(unique(communities), function(x) sum(communities==x))
  cuts = 100
  levelplot(Adj_matr[comms_order,comms_order], at = c(seq(min(Adj_matr),-1e-10,length.out =cuts),0,seq(1e-10,max(Adj_matr),length.out =cuts)),
            xlab = "Brain regions", ylab = "Brain regions",
            col.regions = c(rgb((cuts:0)/cuts, green = 0, blue = 1,red = 0),rgb(red = 0,green = 0, blue = 0),
                            rgb((0:cuts)/cuts, green = 0, blue = 0,red = 1)),
            panel = function(...){
              panel.levelplot(...)
              for(u in Reduce('+',length_comms,accumulate = T)) {
                panel.abline(h = u+0.5)
                panel.abline(v = u+0.5) 
              }
            },scales=list(x=list(at=Reduce('+',c(0,length_comms[1:13]),accumulate = T ) + length_comms[1:14]/2, 
                                 labels=c(1:13,-1)),
                          y=list(at=Reduce('+',c(0,length_comms[1:13]),accumulate = T ) + length_comms[1:14]/2, 
                                 labels=c(1:13,-1))))
}

plot_node_degree <- function(beta, type = "intersection") {
  if(type=="intersection") {
    NODES <- (1+sqrt(1+8*length(beta)))/2
    Adj_matr = array(0,dim = c(NODES, NODES))
    Adj_matr[upper.tri(Adj_matr)] <- beta
    Adj_matr <- Adj_matr + t(Adj_matr)  
  }else{if(type=="union") {
    NODES <- (1+sqrt(1+8*length(beta)/2))/2
    D <- construct_D(NODES)
    v = crossprod(D,beta)
    Adj_matr = array(0,dim = c(NODES, NODES))
    Adj_matr[upper.tri(Adj_matr)] = as.vector(v)
    Adj_matr = Adj_matr+t(Adj_matr)
  }else{
    stop("The value of type should be one between \"intersection\" and \"union\"")
  }}
  barplot(apply(Adj_matr!=0,1,sum), main = "Node degree")
}
