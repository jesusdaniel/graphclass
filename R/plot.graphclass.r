#' Plot function for graph classifier.
#'
#' Plots the adjacency matrix of the coefficients network
#' 
# @rdname graphclass
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
