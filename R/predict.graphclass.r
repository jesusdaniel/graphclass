#' Predict function for graph classifier.
#'
# @rdname graphclass
#' @export
#'
#' @param object trained graphclass object
#' @param newdata matrix of observations to predict. Each row corresponds to a new observation.
#' @param type type of response. class: predicted classes. prob: predicted probabilities. error: misclassification error
#' @param Ytest if type = "error", true classes to compare.
#' @return A vector containing the predicted classes.
#' @examples
#' X = matrix(rnorm(100*34453), nrow = 100)
#' Y = 2*(runif(100) > 0.5) - 1
#' gc = graphclass(X, Y = factor(Y))
#' Xtest = matrix(rnorm(100*34453), nrow = 100)
#' predictions = predict(gc, Xtest)
predict.graphclass <- function(object, newdata, type = "class", Ytest, ...) {
  if (!inherits(object, "graphclass"))  {
    stop("Object not of class 'graphclass'")
  }
  pred <- newdata %*% object$beta + object$b
  Ypred <- sapply(c(pred), function(y) if(y>0) { object$Ypos_label}else{ object$Yneg_label})
  if(type=="class")
    return(Ypred)
  if(type=="pred")
    return(exp(pred) / (1 + exp(pred)))
  if(type == "error")
    return(sum(Ytest!=Ypred)/length(Ypred))
}