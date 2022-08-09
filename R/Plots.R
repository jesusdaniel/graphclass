#' Plots an adjacency matrix.
#'
#' Draws a plot of the adjacency matrix represented by a vector or a matrix
#' using the function \code{levelplot} of the \code{lattice} package.
#' 
#' @export
#' 
#' @param edgevalues Edge values of the adjacency matrix. The argument can be either a square adjacency matrix, or a 
#' vectorized adjacency matrix. For undirected networks, the vector should
#' contain the upper triangular part in column-major order. For directed networks, the vector should contain all
#' entries of the matrix in column major order.
#' @param type If \code{edgevalues} is a vector, this parameter specifies whether the vector represents 
#' an \code{undirected} or \code{directed} network. Default is \code{undirected}.
#' @param edgetype This parameter specifies the type of edge values for the scale color of the plot. 
#' For real-valued edges, \code{type = "real"}.
#' If the edges are between 0 and 1, \code{type = "prob"}. For binary edges, \code{type = "binary"}.
#' @param communities Optional argument to specify a list in which each element contains an array indicating
#' the indexes of the nodes on each community.
#' @param community_labels Labels for each community. The array should have the same length than \code{communities}.
#' @param main Title of the plot
#' @param axislabel Label for the axes.
#' @param colorlims An array with two elements indicating the minimum and maximum value in the color bar.
#' @param selfloops If the \code{edgevalues} parameter is a vector, indicates whether self loops are included in the
#' entries of a vector. Default is \code{FALSE}.
#' 
#' @return An objected returned by the function \code{levelplot} of the \code{lattice} package.
#' 
#' @examples
#' 
#' # Plot the adjacency matrix of a COBRE data subject
#' data(COBRE.data)
#' X1 <- COBRE.data$X.cobre[1,]
#' 
#' # Plot adjacency matrix with nodes in default order
#' plot_adjmatrix(X1)
#' 
#' # Plot adj. matrix divided by communities
#' data(power.parcellation)
#' # Node assignments (note that node 75 is missing on COBRE)
#' node.assignments <- power.parcellation$Master.Assignment[-75]
#' communities = lapply(c(1:13, -1), function(x) which(node.assignments==x))
#' plot_adjmatrix(X1, communities = communities, 
#'       community_labels = c(1:13, -1), axislabel = "Brain systems")
#' 
#' @encoding UTF-8
#' @importFrom Rdpack reprompt
plot_adjmatrix <- function(edgevalues, type=c("undirected", "directed"),
                           edgetype = c("real", "prob", "binary"),
                           communities = NULL, 
                           community_labels = NULL, 
                           main= "", axislabel = "Nodes",
                           colorlims = NULL,
                           selfloops = FALSE) 
  {
  
  type <- match.arg(type)
  require(lattice)
  require(Matrix)
  edgetype <- match.arg(edgetype)
  
  if(is.null(dim(edgevalues))) {
    Adj_matr <- as.matrix(get_matrix(edgevalues, type, selfloops))
  }else{
    Adj_matr <- as.matrix(edgevalues)
  }
  NODES <- ncol(Adj_matr)
  
  tckspace <- max(1, round(NODES/5, max(0, floor(log10(NODES/5)))))
  
  cuts <- 1000
  colorkey <- TRUE
  #edgetype = real----------------------------------------
  atneg <- 0
  atpos <- 0
  col.regions.neg <- rgb(red = 1,green = 1, blue = 1)
  col.regions.pos <- rgb(red = 1,green = 1, blue = 1)
  if(is.null(colorlims)) {
    min_Adj <- min(Adj_matr[!is.na(Adj_matr)])
    max_Adj <- max(Adj_matr[!is.na(Adj_matr)])
  }else{
    min_Adj <- colorlims[1]
    max_Adj <- colorlims[2]
  }
  if(min_Adj <0 ) {
    atneg <- seq(min_Adj ,max(-1e-16,min(Adj_matr)), length.out =cuts)
    col.regions.neg <- rgb((cuts:0)/cuts, green = 0, blue = 1,red = 0)
  }
  if(max_Adj > 0) {
    atpos <- seq(1e-16,max_Adj,length.out =cuts)
    col.regions.pos <- rgb((0:cuts)/cuts, green = 0, blue = 0,red = 1)
  }
  atval = unique(c(atneg, 0, atpos))
  col.vals = (c(col.regions.neg, rgb(red = 1,green = 1, blue = 1), col.regions.pos))
  #edgetype = prob----------------------------------------
  if(edgetype == "prob") {
    atval = seq(0,1,length.out =cuts)
    col.vals = rgb((0:cuts)/cuts, green = 0, blue = 1,red = 0)
  }
  #edgetype = binary -------------------------------------
  if(edgetype == "binary") {
    colorkey <- FALSE
    Adj_matr <- 1*(Adj_matr != 0)
    atval <- seq(0,1,length.out =cuts)
    col.vals <- rgb((0:cuts)/cuts, green = 0, blue = 1,red = 0)
  }
  
  # Plot community groupings --------------------------------
  if(!is.null(communities)) {
    lengths_coms = sapply(communities, length)
    scales_list = list(tck = c(0,0),
                       x=list(at=Reduce('+',c(0,lengths_coms[1:(length(lengths_coms)-1)]),accumulate = T ) + 
                                lengths_coms/2, 
                              labels=community_labels),
                       y=list(at=Reduce('+',c(0,rev(lengths_coms)[1:(length(lengths_coms)-1)]),accumulate = T ) + 
                                rev(lengths_coms)/2, 
                              labels=rev(community_labels)))
    panel_func = function(...){ panel.levelplot(...)
      if(type=="prob_cells") {
        select_list <- which(sel_cells,arr.ind = T)
        for(cell in 1:nrow(select_list)) {
          fill_block(select_list[cell,1],select_list[cell,2], communities)
        }
      }
      
      for(u in Reduce('+',lengths_coms,accumulate = T)) {
        panel.abline(v = u+0.5)
      }
      for(u in Reduce('+',rev(lengths_coms),accumulate = T)) {
        panel.abline(h = u+0.5)
      }
    }
    nodeorder = unlist(communities)
    Adj_matr <- Adj_matr[nodeorder, nodeorder]
    Adj_matr <- Adj_matr[,seq(from=NODES,to=1,by=-1)] #reversing the columns
    levelplot(Adj_matr, at = atval,
              xlab = axislabel, ylab = axislabel,
              main = main,
              colorkey = colorkey,
              col.regions = col.vals,
                       panel = panel_func,
                       scales=scales_list)
  }else{
    Adj_matr <- Adj_matr[,seq(from=NODES,to=1,by=-1)] #reversing the columns
    levelplot(Adj_matr, at = atval,
              xlab = axislabel, ylab = axislabel,
              col.regions = col.vals,
              main = main,
              colorkey = colorkey,
              scales = list(tck = c(1,0), 
                            x = list(at=seq(0,NODES, by = tckspace)),
                            y = list(at = NODES + 1 - tckspace - seq(0,ncol(Adj_matr), by = tckspace),
                                     labels = (seq(tckspace,ncol(Adj_matr),tckspace)))))
  }
}




plot_node_degree_distribution <- function(beta, type = "intersection") {
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
