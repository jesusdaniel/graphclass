#' Prediction function for graph classifier
#' 
#' Given an object of type \code{graphclass}, this function obtains the fitted classes for a new data.
#' @export
#' 
#' @param gc A trained graph classifier object
#' @param X A matrix containing rows with vectorized upper triangular adjacency matrices (column-major order)
#' @param type Indicates the type of response: \code{class} for class predictions, \code{prob} for predicted probabilities, 
#' \code{error} for misclassification error (if \code{Ytest} is provided).
#' @param Ytest If type = "error", true classes to compare.
#' @return A vector containing the predicted classes.
#' @examples
#' data(COBRE.data)
#' X <- COBRE.data$X.cobre
#' Y <- COBRE.data$Y.cobre
#' 
#' #Split into train and test
#' test <- runif(length(Y)) <= 0.1
#' gc <- graphclass(X = X[!test, ], Y = factor(Y[!test]), type = "intersection",
#'                lambda = 1e-4, rho = 1, gamma = 1e-5)
#'                
#' gc.test <- predict(gc, X[test, ]) 
#' predict(gc, X[test, ], "error", Ytest = Y[test])
predict.graphclass <- function(object, newdata, 
                               type = c("class", "prob", "error"), 
                               Ytest, ...) {
  type <- match.arg(type)
  if (!inherits(object, "graphclass"))  {
    stop("Object not of class 'graphclass'")
  }
  pred <- newdata %*% object$beta + object$b
  Ypred <- sapply(c(pred), function(y) if(y>0) { object$Ypos_label}else{ object$Yneg_label})
  if(type=="class")
    return(Ypred)
  if(type=="prob")
    return(exp(pred) / (1 + exp(pred)))
  if(type == "error")
    return(sum(Ytest!=Ypred)/length(Ypred))
}