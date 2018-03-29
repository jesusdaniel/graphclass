#################################### Signal- subgraph

train_signalsubgraph_cv <- function(X, Y, fold_list = NULL, folds = 10, grid_params, algorithm_parameters, verbose = F,
                        truepos = NULL, truenodes = NULL, threshold_bin = NULL) {
  Xbin = 1*(X > threshold_bin)
  if(is.null(fold_list) & folds > 0) {
    require(cvTools)
    fold_partition = cvFolds(n = nrow(X), K = folds)
    fold_list <- lapply(1:fold_partition$K, function(j) 
      fold_partition$subsets[which(fold_partition$which==j)])
  }
  foldmatrix = sapply(fold_list, function(x) {u = rep(0, length(Y)); u[x] = 1; return(u)})
  Xadj = apply(Xbin, 1, get_matrix)
  NODES <- (1+sqrt(1+8*ncol(X)))/2
  Xadj = array(Xadj, dim = c(NODES, NODES, length(Y)))
  require(R.matlab)
  writeMat(X = Xadj, Y = 1*(Y>0), foldmatrix = foldmatrix,  truepos = get_matrix(truepos), truenodes = 1*truenodes,
           con =  paste("SignalSubgraph/SS_train_file_", algorithm_parameters$fileid,".mat", sep = ""))
  require(matlabr)
  params = paste("[",Reduce(function(x,y) paste(x,";",y),sapply(grid_params, function(x) paste(x[1],", ",x[2]))),"]")
  
  code = c("cd 'SignalSubgraph'",
           paste("load('SS_train_file_", algorithm_parameters$fileid, ".mat')", sep=""),
           paste("[errors, sp, nodesp,TPR, FPR, nodeTPR, nodeFPR] = train_signalsubgraph(X, Y, ", params, ", foldmatrix, truepos,truenodes);"),
           paste("save('SS_train_result_",algorithm_parameters$fileid, ".txt', 'errors', 'sp', 'nodesp',",
                 "'TPR', 'FPR','nodeTPR','nodeFPR', '-ascii')", sep=""))
  res = run_matlab_code(code)
  while(!file.exists(paste("SignalSubgraph/SS_train_result_",algorithm_parameters$fileid, ".txt", sep=""))){}
  results <- as.matrix(read.table(paste("SignalSubgraph/SS_train_result_",algorithm_parameters$fileid, ".txt", sep=""),header=F))
  file.remove(paste("SignalSubgraph/SS_train_result_",algorithm_parameters$fileid, ".txt", sep=""))
  file.remove(paste("SignalSubgraph/SS_train_file_",algorithm_parameters$fileid,".mat", sep = ""))
  errors = results[1:length(grid_params),]
  sparsity_table = results[(length(grid_params)+1):(2*length(grid_params)),]
  node_sparsity_table = results[(2*length(grid_params)+1):(3*length(grid_params)),]
  
  list(trained_classifiers = list(X = X, Y = Y, fileid = algorithm_parameters$fileid, grid_params = grid_params, threshold_bin = threshold_bin),
       errors = errors,
       cv_errors = apply(errors, 1, mean),
       sd_errors = apply(errors, 1, sd),
       sparsity_table = sparsity_table,
       node_sparsity_table = node_sparsity_table,
       grid_params = grid_params)
}

train_signalsubgraph_select <- function(X, Y, grid_params, truepos = NULL, truenodes, threshold_bin = NULL, fileid = "") {
  if(!is.null(threshold_bin))
    X = 1*(X>threshold_bin)
  Xadj = apply(X, 1, get_matrix)
  NODES <- (1+sqrt(1+8*ncol(X)))/2
  Xadj = array(Xadj, dim = c(NODES, NODES, length(Y)))
  require(R.matlab)
  writeMat(X = Xadj, Y = 1*(Y>0), truepos = get_matrix(truepos), truenodes = 1*truenodes,
           con =  paste("SignalSubgraph/select_variables_file_", fileid,".mat", sep = ""))
  require(matlabr)
  
  params = paste("[",Reduce(function(x,y) paste(x,";",y),sapply(grid_params, function(x) paste(x[1],", ",x[2]))),"]")
  
  code = c("cd 'SignalSubgraph'",
           paste("load('select_variables_file_",fileid, ".mat')", sep=""),
           paste("[TPR, FPR, nodeTPR, nodeFPR, sparsity, nodesparsity] = select_variables(", params, ", X,Y, truepos,truenodes);"),
           paste("save('select_result_",fileid, ".txt', 'TPR', 'FPR','nodeTPR','nodeFPR', 'sparsity', 'nodesparsity', '-ascii')", sep="")) 
  res = run_matlab_code(code)
  while(!file.exists(paste("SignalSubgraph/select_result_",fileid, ".txt", sep=""))){}
  results <- as.matrix(read.table(paste("SignalSubgraph/select_result_",fileid, ".txt", sep=""),header=F))
  file.remove(paste("SignalSubgraph/select_result_",fileid, ".txt", sep=""))
  file.remove(paste("SignalSubgraph/select_variables_file_", fileid,".mat", sep = ""))
  return(list(TPR = results[1,], FPR = results[2,], nodeTPR = results[3,], nodeFPR = results[4,], sparsity = results[5,],
              node_sparsity = results[6,]))
}

test_signalsubgraph <- function(train_obj, X, Y) {
  Xbin = 1*(train_obj$X>train_obj$threshold_bin)
  Xadj = apply(Xbin, 1, get_matrix)
  NODES <- (1+sqrt(1+8*ncol(train_obj$X)))/2
  Xadj = array(Xadj, dim = c(NODES, NODES, length(train_obj$Y)))
  Xbin_test = 1*(X>train_obj$threshold_bin)
  Xadj_test = apply(Xbin_test, 1, get_matrix)
  Xadj_test = array(Xadj_test, dim = c(NODES, NODES, length(Y)))
  require(R.matlab)
  writeMat(X = Xadj, Y = 1*(train_obj$Y>0), Xtest = Xadj_test, Ytest = 1*(Y>0),
           con =  paste("SignalSubgraph/File_test_signalsubgraph_", train_obj$fileid,".mat", sep = ""))
  require(matlabr)
  
  params = paste("[",Reduce(function(x,y) paste(x,";",y),sapply(train_obj$grid_params, 
                                                                function(x) paste(x[1],", ",x[2]))),"]")
  code = c("cd 'SignalSubgraph'",
           paste("load('File_test_signalsubgraph_", train_obj$fileid,".mat')", sep=""),
           paste("params = mat2cell(", params, ", ones(", length(train_obj$grid_param), ",1),2)",sep="" ),
           paste("errors = test_signal_subgraph(X, Y, params, Xtest,Ytest);"),
           paste("save('test_result_",train_obj$fileid, ".txt', 'errors', '-ascii')", sep="")) 
  res = run_matlab_code(code)
  while(!file.exists(paste("SignalSubgraph/test_result_",train_obj$fileid, ".txt", sep=""))){}
  results <- as.matrix(read.table(paste("SignalSubgraph/test_result_",train_obj$fileid, ".txt", sep=""),header=F))
  file.remove(paste("SignalSubgraph/test_result_",train_obj$fileid, ".txt", sep=""))
  file.remove(paste("SignalSubgraph/File_test_signalsubgraph_", train_obj$fileid,".mat", sep = ""))
  return(results)
}
###################################### SVML1

# Grid_params: 
# 1. l1
train_svml1 <- function(X, Y, grid_params, algorithm_parameters, verbose = F,
                        truepos = NULL, truenodes = NULL) {
  require(penalizedSVM)
  trained_classifiers <- lapply(grid_params, function(l1) 
    svm.fs(x = X, y = factor(Y), cross.outer= 0, grid.search = "discrete",  
           fs.method = "1norm", lambda1.set = l1, lambda2.set = c(0), parms.coding = "none",
           inner.val.method = "cv", cross.inner = 0, verbose = verbose))
  sparse <- sapply(trained_classifiers, function(svml1) length(svml1$model$xind))
  beta = rep(0, ncol(X))
  node_sparse <- sapply(trained_classifiers, function(svml1) {
    beta[svml1$model$xind] = 1
    node_sparsity(beta)
  })
  if(!is.null(truepos)) {
    TPR <- sapply(trained_classifiers, function(svml1) {
      beta[svml1$model$xind] = 1
      sum(beta*truepos) / sum(truepos)
    })
    FPR = sapply(trained_classifiers, function(svml1) {
      beta[svml1$model$xind] = 1
      sum(beta*(1-truepos)) / sum(1-truepos)
    })
    nodeTPR <- sapply(trained_classifiers, function(svml1) {
      beta[svml1$model$xind] = 1
      sum((apply(get_matrix(beta),1,sum)>0)*truenodes) / sum(truenodes)
    })
    nodeFPR <- sapply(trained_classifiers, function(svml1) {
      beta[svml1$model$xind] = 1
      sum((apply(get_matrix(beta),1,sum)>0)*(1-truenodes)) / sum(1-truenodes)
    })
  }else{
    TPR = NULL; nodeTPR = NULL
    FPR = NULL; nodeFPR = NULL
  }
  list(trained_classifiers = trained_classifiers,
       sparsity = sparse,
       node_sparsity = node_sparse,
       grid_params = grid_params,
       TPR = TPR, FPR = FPR, nodeTPR=nodeTPR, nodeFPR=nodeFPR)
}

test_svml1 <- function(train_obj, X, Y) {
  require(penalizedSVM)
  sapply(train_obj$trained_classifiers, function(svml1) {
    predict(svml1, newdata = X, newdata.labels = factor(Y))$error
  })
}
###################################### Naive Bayes

# Grid_params: list of number of variables
train_naivebayes <- function(X, Y, grid_params, algorithm_parameters, verbose = F, truepos = NULL, truenodes = NULL) {
  #source("https://bioconductor.org/biocLite.R")
  #biocLite("genefilter")
  require(genefilter)
  tt <-  colttests(X, factor(Y), tstatOnly = T)[,2]
  ordvars <- order(abs(tt), decreasing = T)
  require(naivebayes)
  trained_classifiers <- lapply(grid_params, function(numvars) {
    list(naivebayesclass = naive_bayes(x = X[,ordvars[1:numvars]], y = factor(Y)),
         list_vars = ordvars[1:numvars])
  })
  sparse <- unlist(grid_params)
  beta = rep(0, ncol(X))
  node_sparse <- sapply(grid_params, function(numvars) {
    beta[ordvars[1:numvars]] = 1
    node_sparsity(beta)
  })
  if(!is.null(truepos)) {
    TPR <- sapply(grid_params, function(numvars) {
      beta[ordvars[1:numvars]] = 1
      sum(beta*truepos) / sum(truepos)
    })
    FPR = sapply(grid_params, function(numvars) {
      beta[ordvars[1:numvars]] = 1
      sum(beta*(1-truepos)) / sum(1-truepos)
    })
    nodeTPR <- sapply(grid_params, function(numvars) {
      beta[ordvars[1:numvars]] = 1
      sum((apply(get_matrix(beta),1,sum)>0)*truenodes) / sum(truenodes)
    })
    nodeFPR <- sapply(grid_params, function(numvars) {
      beta[ordvars[1:numvars]] = 1
      sum((apply(get_matrix(beta),1,sum)>0)*(1-truenodes)) / sum(1-truenodes)
    })
  }else{
    TPR = NULL
    FPR = NULL
    nodeTPR = NULL
    nodeFPR = NULL
  }
  list(trained_classifiers = trained_classifiers, 
       sparsity = sparse,
       node_sparsity = node_sparse,
       grid_params = grid_params,
       TPR = TPR, FPR = FPR,
       nodeTPR = nodeTPR, nodeFPR = nodeFPR)
}

train_naivebayes_select <- function(X, Y, grid_params, truepos = NULL, truenodes) {
  require(genefilter)
  tt <-  colttests(X, factor(Y), tstatOnly = T)[,2]
  ordvars <- order(abs(tt), decreasing = T)
  beta = rep(0, ncol(X))
  TPR <- sapply(grid_params, function(numvars) {
    beta[ordvars[1:numvars]] = 1
    sum(beta*truepos) / sum(truepos)
  })
  FPR = sapply(grid_params, function(numvars) {
    beta[ordvars[1:numvars]] = 1
    sum(beta*(1-truepos)) / sum(1-truepos)
  })
  nodeTPR <- sapply(grid_params, function(numvars) {
    beta[ordvars[1:numvars]] = 1
    sum((apply(get_matrix(beta),1,sum)>0)*truenodes) / sum(truenodes)
  })
  nodeFPR <- sapply(grid_params, function(numvars) {
    beta[ordvars[1:numvars]] = 1
    sum((apply(get_matrix(beta),1,sum)>0)*(1-truenodes)) / sum(1-truenodes)
  })
  return(list(TPR = TPR, FPR = FPR, nodeTPR = nodeTPR, nodeFPR = nodeFPR))
}
test_naivebayes <- function(train_obj, X, Y) {
  require(naivebayes)
  sapply(train_obj$trained_classifiers, function(nb) {
    Yhat <- predict(nb$naivebayesclass, newdata = X[, nb$list_vars], type = "class")
    sum(Yhat != factor(Y)) / length(Y)
  })
}

################################################################################## SVM
# Grid_params: 
# 1. Cost
# alg params : kernel (linear)
train_svm <- function(X, Y, grid_params, algorithm_parameters, verbose = F) {
  require(e1071)
  trained_classifiers <- lapply(grid_params, function(cost) 
    svm(x = X, y = factor(Y), kernel = algorithm_parameters, cost = cost))
  list(trained_classifiers = trained_classifiers, grid_params = grid_params, 
       algorithm_parameters = algorithm_parameters)
}

test_svm <- function(train_obj, X, Y) {
  require(e1071)
  sapply(train_obj$trained_classifiers, function(svm) {
    Yhat <- predict(svm, newdata = X, type = "class")
    sum(Yhat != factor(Y)) / length(Y)
  })
}


############################################################################ Random forest
# Grid_params: 
# 1. number of variables
# 2. number of trees
train_randomforest <- function(X, Y, grid_params, algorithm_parameters, verbose = F) {
  require(randomForest)
  trained_classifiers <- lapply(grid_params, function(pars) 
    randomForest(x = X, y = factor(Y), mtry = pars[1], ntree = pars[2]))
  list(trained_classifiers = trained_classifiers, grid_params = grid_params)
}

test_randomforest <- function(train_obj, X, Y) {
  require(randomForest)
  sapply(train_obj$trained_classifiers, function(rf) {
    Yhat <- predict(rf, newdata = X, type = "class")
    sum(Yhat != factor(Y)) / length(Y)
  })
}


###################################### GC

# Grid_params: 
# 1. lambda 2. rho 3. gamma
# algorithm parameters :  D
train_graphclass <- function(X, Y, grid_params, algorithm_parameters, verbose = F, truepos = NULL, truenodes= NULL) {
  require(graphclass)
  trained_classifiers <- lapply(grid_params, function(params) {
    if(verbose) print("GC parameter done!")
    graphclass(X = X, Y = Y, Xtest = X, Ytest = Y, type = "intersection",
               lambda = params[1], rho = params[2],  gamma = params[3],
               params = list(MAX_ITER =200, CONV_CRIT = 1e-05, MAX_TIME = Inf),
               id = "", verbose = F,       # verbose = T, then prints progress on every iteration, id identifies job when multiple jobs are running
               D = algorithm_parameters$D)})
  sparse <- sapply(trained_classifiers, function(gc) sum(gc$beta!=0))
  node_sparse <- sapply(trained_classifiers, function(gc) node_sparsity(gc$beta))
  if(!is.null(truepos)) {
    TPR <- sapply(trained_classifiers, function(gc) {
      beta = gc$beta!=0
      sum(beta*truepos) / sum(truepos)
    })
    FPR = sapply(trained_classifiers, function(gc) {
      beta = gc$beta!=0
      sum(beta*(1-truepos)) / sum(1-truepos)
    })
    nodeTPR <- sapply(trained_classifiers, function(gc) {
      beta = gc$beta!=0
      sum((apply(get_matrix(beta),1,sum)>0)*truenodes) / sum(truenodes)
    })
    nodeFPR <- sapply(trained_classifiers, function(gc) {
      beta = gc$beta!=0
      sum((apply(get_matrix(beta),1,sum)>0)*(1-truenodes)) / sum(1-truenodes)
    })
  }else{
    TPR = NULL; nodeFPR = NULL
    FPR = NULL; nodeTPR = NULL
  }
  list(trained_classifiers = trained_classifiers,
       sparsity = sparse,
       node_sparsity = node_sparse,
       grid_params = grid_params,
       TPR = TPR, FPR = FPR, nodeTPR = nodeTPR, nodeFPR = nodeFPR)
}

test_graphclass <- function(train_obj, X, Y) {
  require(graphclass)
  sapply(train_obj$trained_classifiers, function(gc) {
    predict(gc, newdata = X, Ytest = Y, type = "error")
  })
}



###################################### Relaxed GC
# Grid_params: 
# 1. lambda 2. rho 3. gamma
# algorithm parameters :  D
train_relaxed_graphclass <- function(X, Y, grid_params, algorithm_parameters, verbose = F) {
  require(graphclass)
  trained_classifiers <- lapply(grid_params, function(params)  {
    gc = graphclass(X = X, Y = Y, Xtest = X, Ytest = Y, type = "intersection",
               lambda = params[1], rho = params[2],  gamma = params[3],
               params = list(MAX_ITER =200, CONV_CRIT = 1e-05, MAX_TIME = Inf),
               id = "", verbose = F,       # verbose = T, then prints progress on every iteration, id identifies job when multiple jobs are running
               D = algorithm_parameters$D)
    coeffs = which(gc$beta!=0)
    require(glmnet)
    relaxed_gc <- glmnet(x = X[,coeffs], y = factor(Y), family = "binomial", alpha = 0, lambda = params[3])
    return(list(gc = gc, relaxed_gc = relaxed_gc))
  })
  sparse <- sapply(trained_classifiers, function(tc) sum(tc$gc$beta!=0))
  node_sparse <- sapply(trained_classifiers, function(tc) node_sparsity(tc$gc$beta))
  list(trained_classifiers = trained_classifiers,
       sparsity = sparse,
       node_sparsity = node_sparse,
       grid_params = grid_params)
}

test_relaxed_graphclass <- function(train_obj, X, Y) {
  require(graphclass)
  sapply(train_obj$trained_classifiers, function(tc) {
    Yhat <- predict(tc$relaxed_gc, newx = X[, which(tc$gc$beta!=0)], type = "class")
    sum(Yhat != factor(Y)) / length(Y)
  })
}



###################################### Glmnet
# Grid_params: 
# 1. lambda 
# algorithm parameters :  alpha
train_glmnet <- function(X, Y, grid_params, algorithm_parameters, verbose = F, truepos = NULL, truenodes = truenodes) {
  require(glmnet)
  trained_classifiers <- lapply(grid_params, function(params)  {
    gc = glmnet(x = X,y = Y, family = "binomial", alpha = algorithm_parameters$alpha, lambda = params)
    coeffs = which(gc$beta!=0)
    require(glmnet)
    return(list(gc = gc))
  })
  sparse <- sapply(trained_classifiers, function(tc) sum(tc$gc$beta!=0))
  node_sparse <- sapply(trained_classifiers, function(tc) node_sparsity(tc$gc$beta))
  if(!is.null(truepos)) {
    TPR <- sapply(trained_classifiers, function(tc) {
      beta = tc$gc$beta!=0
      sum(beta*truepos) / sum(truepos)
    })
    FPR = sapply(trained_classifiers, function(tc) {
      beta = tc$gc$beta!=0
      sum(beta*(1-truepos)) / sum(1-truepos)
    })
    nodeTPR <- sapply(trained_classifiers, function(tc) {
      beta = tc$gc$beta!=0
      sum((apply(get_matrix(beta),1,sum)>0)*truenodes) / sum(truenodes)
    })
    nodeFPR <- sapply(trained_classifiers, function(tc) {
      beta = tc$gc$beta!=0
      sum((apply(get_matrix(beta),1,sum)>0)*(1-truenodes)) / sum(1-truenodes)
    })
  }else{
    TPR = NULL; nodeTPR = NULL
    FPR = NULL; nodeFPR = NULL
  }
  list(trained_classifiers = trained_classifiers,
       sparsity = sparse,
       node_sparsity = node_sparse,
       grid_params = grid_params, TPR = TPR, FPR = FPR,
       nodeTPR = nodeTPR, nodeFPR = nodeFPR)
}

train_glmnet_select <- function(X, Y, grid_params = NULL, truepos = NULL, truenodes = truenodes, alpha = 1) {
  require(glmnet)
  glmns <- glmnet(x = X, y = factor(Y), lambda = grid_params, family = "binomial", alpha = alpha)
  TPR <- apply(glmns$beta, 2, function(beta) {
    sum((beta!=0)*truepos) / sum(truepos)
  })
  FPR = apply(glmns$beta, 2, function(beta) {
    sum((beta!=0)*(1-truepos)) / sum(1-truepos)
  })
  nodeTPR <- apply(glmns$beta, 2, function(beta) {
    sum((apply(get_matrix(abs(beta)),1,sum)>0)*truenodes) / sum(truenodes)
  })
  nodeFPR <- apply(glmns$beta, 2, function(beta) {
    sum((apply(get_matrix(abs(beta)),1,sum)!=0)*(1-truenodes)) / sum(1-truenodes)
  })
  return(list(TPR = TPR, FPR = FPR,nodeTPR = nodeTPR, nodeFPR = nodeFPR))
}

test_glmnet <- function(train_obj, X, Y) {
  require(graphclass)
  sapply(train_obj$trained_classifiers, function(tc) {
    Yhat <- predict(tc$gc, newx = X, type = "class")
    sum(Yhat != factor(Y)) / length(Y)
  })
}



###################################### Global summaries
# Grid_params: 
# 1. lambda 
# algorithm parameters :  alpha
train_glm <- function(X, Y, grid_params, algorithm_parameters, verbose = F) {
  Y = ifelse(Y==max(as.vector(Y)),1,0)
  trained_classifiers <- glm(formula = Y ~ X, family = "binomial")
  list(trained_classifiers = trained_classifiers)
}

test_glm <- function(train_obj, X, Y) {
  #browser()
  Yhat <- 1*((cbind(1, X) %*% train_obj$trained_classifiers$coefficients)>0)
  Y = ifelse(Y==max(as.vector(Y)),1,0)
  sum(Yhat != (Y)) / length(Y)
}
