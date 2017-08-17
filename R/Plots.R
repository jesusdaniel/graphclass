#' Plot a vectorized adjacency matrix.
#'
#' @param beta Vectorized adjacency matrix. For undirected networks use only upper triangle in column-major order, for directed use both
#' @param type Either intersection for undirected networks, union for directed.
#' @examples
#' B = runif(34453)
#' plot_adjmatrix(B)
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
  print(levelplot(df, at = unique(c(seq(min(Adj_matr),max(-1e-10,min(Adj_matr)),length.out =cuts),0,seq(1e-10,max(Adj_matr),length.out =cuts))),
            xlab = "Nodes", ylab = "Nodes",
            col.regions = c(rgb((cuts:0)/cuts, green = 0, blue = 1,red = 0),rgb(red = 0,green = 0, blue = 0),
                            rgb((0:cuts)/cuts, green = 0, blue = 0,red = 1)),
            scales = list(tck = c(1,0), 
                          x = list(at=seq(0,ncol(df),round(NODES/100)*10)),
                          y = list(at = NODES-round(NODES/100)*10- seq(0,ncol(df),round(NODES/100)*10),
                                     labels = (seq(round(NODES/100)*10,ncol(df),round(NODES/100)*10))))))
                            
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
  print(levelplot(Adj_matr[comms_order,comms_order], at = c(seq(min(Adj_matr),-1e-10,length.out =cuts),0,seq(1e-10,max(Adj_matr),length.out =cuts)),
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
                                 labels=c(1:13,-1)))))
}

plot_node_degree <- function(beta, type = "intersection") {
  if(type=="intersection") {
    NODES <- (1+sqrt(1+8*length(beta)))/2
    Adj_matr = array(0,dim = c(NODES, NODES))
    Adj_matr[upper.tri(Adj_matr)] <- beta
    Adj_matr <- Adj_matr + t(Adj_matr)  
  }else{if(type=="union") {
    NODES <- (1+sqrt(1+8*length(beta)/2))/2
    Adj_matr = array(T,dim = c(NODES, NODES))
    diag(Adj_matr) = F
    Adj_matr[Adj_matr] = beta
  }else{
    stop("The value of type should be one between \"intersection\" and \"union\"")
  }}
  print(barplot(apply(Adj_matr!=0,2,sum), main = "Node degree"))
}

plot_node_value <- function(node_vals, 
                            communities = NULL, type = "pval", community_labels = c(1:13,-1),main="") {
  node_vals <- ifelse(node_vals ==0, 0.01/length(node_vals), node_vals)
  top_size <- -log(0.01/length(node_vals))
  barplot(-log(node_vals), ylab="-log(P-values)", xlab="Brain regions", border="white",
          ylim=c(0,top_size), axes=FALSE, space=0, col="grey50")
  lengths_rois = sapply(communities, length)
  pos_rois = Reduce('+',c(0,lengths_rois),accumulate = T )
  at_tick = Reduce('+',c(0,lengths_rois[1:(length(communities)-1)]),accumulate = T ) + lengths_rois/2
  
  sapply(2*(0:((length(communities)-1)/2))+1, 
         function(i) rect(pos_rois[i], 0, pos_rois[i+1], -log(0.01/length(node_vals)), border = "gray80", col = "gray80"))
  barplot(-log(node_vals), ylim = c(0,top_size), main = main,
          border="black",ylab="", xlab="", axes=FALSE, space=0, col="gray50",add=T)
  axis(2, pos=0, at = c(-log(0.05), -log(0.05/length(node_vals))), labels = c("-log(0.05)", "-log(Bonf)"))
  axis(1, at = at_tick, pos=0, labels  =community_labels)
  segments(sum(lengths_rois),top_size, sum(lengths_rois),0)
  segments(0,top_size, sum(lengths_rois),top_size)
  segments(0,0,0,top_size)
  segments(0,0,sum(lengths_rois),0)
  abline(h = -log(0.05), lty = 2)
  abline(h = -log(0.05/length(node_vals)), lty = 2, col = "red")
}

plot_node_value_2 <- function(node_vals, 
                            communities = NULL, type = "pval", community_labels = c(1:13,-1), main="") {
  nodes_order <- unlist(communities)
  top_size <- 1
  lengths_rois <- sapply(communities, length)
  pos_rois <- Reduce('+',c(0,lengths_rois),accumulate = T )
  at_tick2 = Reduce('+',c(0,lengths_rois[1:(length(communities)-1)]),accumulate = T ) + lengths_rois/2
  at_tick = Reduce('+',c(0,lengths_rois[1:(length(communities))]),accumulate = T )
  cols <- c()
  for(j in 1:length(lengths_rois)){
    cols <- c(cols, rep(ifelse(j%%2==1, "grey10", "grey60"), lengths_rois[j]))
  }
  
  barplot(node_vals[nodes_order], ylab="", xlab="Brain systems", border="white",
          ylim=c(0,top_size), axes=FALSE, space=0, col="grey50")
  u <- sapply(2*(0:((length(communities)-1)/2))+1, 
         function(i) rect(pos_rois[i], 0, pos_rois[i+1], 1, border = "gray80", col = "gray80"))
  barplot(node_vals[nodes_order], ylim = c(0,top_size), main = main,
          border="black",ylab="", xlab="", axes=FALSE, space=0, col=cols,add=T)

  axis(1, at = at_tick, pos=0, labels  = FALSE)
  axis(1, at = at_tick2, tick = FALSE, pos=0, labels  =community_labels)
  axis(2, pos=0)
  segments(sum(lengths_rois),top_size, sum(lengths_rois),0)
  segments(0,top_size, sum(lengths_rois),top_size)
  segments(0,0,0,top_size)
  segments(0,0,sum(lengths_rois),0)
}


#' Returns a matrix from a vectorized network
#'
#' @param beta Vectorized adjacency matrix. 
#' @param type Either intersection for undirected networks, union for directed.
#' @return Adjacency matrix for a vectorized network
get_matrix <- function(beta, type="intersection") {
  if(type=="intersection") {
    NODES <- (1+sqrt(1+8*length(beta)))/2
    Adj_matr = array(0,dim = c(NODES, NODES))
    Adj_matr[upper.tri(Adj_matr)] <- beta
    Adj_matr <- Adj_matr + t(Adj_matr)  
  }else{if(type=="union") {
    NODES <- (1+sqrt(1+8*length(beta)/2))/2
    Adj_matr = array(T,dim = c(NODES, NODES))
    diag(Adj_matr) = F
    Adj_matr[Adj_matr] = beta
  }else{
    stop("The value of type should be one between \"intersection\" and \"union\"")
  }}
  return(Adj_matr)
}

#' Returns node sparsity of a vector
#'
#' @param beta Vectorized adjacency matrix. 
#' @return Percentage of inactive nodes in the graph solution
node_sparsity <- function(beta) {
  A <- get_matrix(beta)
  return(sum(apply(A,1,function(v) sum(v!=0))==0)/ncol(A))
}

#' Plot a vectorized adjacency matrix with cells divisions
#'
#' @param edge_values Vectorized adjacency matrix. Only undirected networks are supported for now.
#' @param communities Community of each node
#' @param community_labels Name of each community that will appear on the plot.
#' @param main Title of the plot
#' @param type Either "real" for valued networks, "prob" for [0,1] valued networks or "prob_cells" for equal value on each cell
plot_square_adj_mat <- function(edge_values, 
                                communities = NULL, type = "real", 
                                community_labels = c(1:13,-1), 
                                main= "", cut_at, sel_cells) {
  edge_values <- edge_values[,seq(from=ncol(edge_values),to=1,by=-1)] #reversing the columns
  require(lattice)
  cuts = 1000
  if(type == "real") {
    at_cuts = c(seq(min(edge_values,0)-1e-08,0,length.out =cuts),seq(1e-016,max(edge_values),length.out =cuts))
    col_cuts = c(rgb((cuts:0)/cuts, green = 0, blue = 0,red = 1),rgb(red = 0,green = 0, blue = 0),
             rgb((0:cuts)/cuts, green = 0, blue = 1,red = 0))
  }else{
    if(type == "pval") {
      at_cuts = seq(0,1,length.out = cuts)
      col_cuts = c(rgb((cuts:0)/cuts, green = 0, blue = 0,red = 1))
    }else{
      if(type == "prob"|type == "prob_cells") {
        at_cuts = seq(0,1,length.out = cuts)
        col_cuts = c(rgb(((0:cuts)/cuts)^3, green = 0, blue = 1,red = 0))
      }else{
        if(type == "prob2") {
          at_cuts = seq(0,1,length.out = cuts)
          col_cuts = c(rgb((seq(0,cuts*cut_at,length.out = cuts*cut_at)/cuts)^3, green = 0, blue = 1,red = 0), 
                       rgb((seq(cuts*cut_at+1e-16,cuts,length.out = cuts*(1-cut_at))/cuts)^2, green = 0, blue = 0,red = 1))
        }
      }
    }
  }
  if(!is.null(communities)) {
    lengths_rois = sapply(communities, length)
    scales_list = list(tck = c(0,0),
                       x=list(at=Reduce('+',c(0,lengths_rois[1:(length(lengths_rois)-1)]),accumulate = T ) + lengths_rois/2, 
                              labels=community_labels),
                       y=list(at=Reduce('+',c(0,rev(lengths_rois)[1:(length(lengths_rois)-1)]),accumulate = T ) + rev(lengths_rois)/2, 
                              labels=rev(community_labels)))
    panel_func = function(...){ panel.levelplot(...)
      if(type=="prob_cells") {
        select_list <- which(sel_cells,arr.ind = T)
        for(cell in 1:nrow(select_list)) {
          fill_block(select_list[cell,1],select_list[cell,2], communities)
        }
      }
        
      for(u in Reduce('+',lengths_rois,accumulate = T)) {
        panel.abline(v = u+0.5)
      }
      for(u in Reduce('+',rev(lengths_rois),accumulate = T)) {
        panel.abline(h = u+0.5)
      }
    }
    plot_sm =levelplot(edge_values, 
                       xlab = "Brain systems", ylab = "Brain systems",
                       at = at_cuts,
                       col.regions = col_cuts,
                       panel = panel_func,
                       scales=scales_list, main = main)
  }else{
    plot_sm =levelplot(edge_values, main = main,
                       xlab = "Brain systems", ylab = "Brain systems",
                       at = at_cuts,
                       col.regions = col_cuts)
  }
  print(plot_sm)
}



fill_block <- function(whichX, whichY, communities) {
  xnodes = 1:length(communities[[whichX]]) + (whichX>1)*length(unlist(communities[1:(whichX-1)]))
  ynodes = 1:length(communities[[whichY]]) + (whichY>1)*length(unlist(communities[1:(whichY-1)]))
  
  fx1 <- xnodes[which((xnodes+min(ynodes))%%3 == 0)]
  fy1 <- rep(min(ynodes), length(fx1))
  fy2 <-  sapply(fx1, function(j) min(j-min(xnodes)+min(ynodes), max(ynodes)))
  fx2 <- sapply(1:length(fx1), function(i) fx1[i] +fy1[i] - fy2[i])
  
  gy1 <- ynodes[which((ynodes+max(xnodes))%%3 == 0)]
  gx1 <- rep(max(xnodes), length(gy1))
  gy2 <-  sapply(gy1, function(j) min(j+max(xnodes)-min(xnodes), max(ynodes)))
  gx2 <- sapply(1:length(gx1), function(i) gx1[i] +gy1[i] - gy2[i])
  
  fx1 <- c(fx1, gx1)
  fx2 <- c(fx2, gx2)
  fy1 <- c(fy1, gy1)
  fy2 <- c(fy2, gy2)
  for (i in seq(length(fx1)))
  {
    panel.linejoin(x = c(fx1[i]+0.5,fx2[i]-0.5), y= length(unlist(communities))-c(fy1[i]-0.5,fy2[i]+0.5)+1,
                   col="black")
  }
}