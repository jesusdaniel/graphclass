### Cross-validation
cross_validation_function <- function(X, Y, parameters_grid, methodname = "",
                                      train_classifier, test_classifier,
                                      algorithm_parameters = "",
                                folds = 5, fold_list = NULL, 
                                parallel = F, num_clusters = 5, windows = F, 
                                nested_cv = F, save_files = F, filename = "", sparsity_results=F) {
  NODES <- (1+sqrt(1+8*ncol(X)))/2

  n <- nrow(X)
  if(is.null(fold_list) & folds > 0) {
    require(cvTools)
    fold_partition = cvFolds(n = nrow(X), K = folds)
    fold_list <- lapply(1:fold_partition$K, function(j) 
      fold_partition$subsets[which(fold_partition$which==j)])
  }else{
    folds <- length(fold_list)
  }
  
  if(nested_cv) {
    combination_folds <- combinat::combn(1:folds,2)
    nestedfolds_subsets <- lapply(1:ncol(combination_folds), function(i) {
      sort((1:n)[c(-fold_list[[combination_folds[1,i]]], -fold_list[[combination_folds[2,i]]])])
    })
  }
  
  ################################################################################
  # Parallel ---------------------------------------------------------------------
  ################################################################################
  if(parallel) {
    cat("Trying parallel with",num_clusters,"cores\n")
    if(!windows) library(Rmpi)
    library(parallel)
    ifelse(windows, cl <- makeCluster(num_clusters),   # use this for computations on Windows
           cl <- makeCluster(num_clusters ,type="MPI"))    # use this for computations on Flux
    
    clusterEvalQ(cl, require(graphclass))
    
    clusterExport(cl, varlist =  list("train_classifier", "test_classifier", "fold_list", "node_sparsity",
                                      "parameters_grid", "algorithm_parameters",
                                      "X", "Y"), envir=environment())
    #clusterExport(cl = cl, varlist = list("X", "Y", algorithm_parameters), envir=environment())
    if(nested_cv) {
      cat("Trying nested...\n")
      
      nested_cross_validation <- parLapplyLB(cl, 1:length(nestedfolds_subsets), function(i)
        train_classifier(X[nestedfolds_subsets[[i]], ], Y[nestedfolds_subsets[[i]]], 
                         parameters_grid, algorithm_parameters, verbose = F))
      nested_errors <- array(NA, dim = c(length(fold_list), length(fold_list), length(parameters_grid)))
      #nested_errors
      for(i in 1:length(fold_list)) {
        #subfold
        for(j in (1:length(fold_list))[-i]) {
          u1 = min(i,j); u2 = max(i,j)
          nested_errors[i,j, ] = test_classifier(nested_cross_validation[[intersect(which(combination_folds[1,]==u1),which(combination_folds[2,]==u2))]], 
                                                 X[fold_list[[j]], ], Y[fold_list[[j]]])
        }
      }
      cat("Nested over!!\n")
    }
    cross_validation <- parLapplyLB(cl, 1:length(fold_list), function(i)
      train_classifier(X[-fold_list[[i]],], Y[-fold_list[[i]]],
                       parameters_grid, algorithm_parameters, verbose = F))
    clusterExport(cl,varlist =  "cross_validation",  envir=environment())
    cross_validation_errors <- parSapplyLB(cl, 1:length(fold_list), function(i) 
      test_classifier(cross_validation[[i]], X[fold_list[[i]], ], Y[fold_list[[i]]]))
    stopCluster(cl)
  }else{
    if(nested_cv) {
      cat("Trying nested...\n")
      nested_cross_validation <- lapply(1:length(nestedfolds_subsets), function(i)
        train_classifier(X[nestedfolds_subsets[[i]], ], Y[nestedfolds_subsets[[i]]], 
                         parameters_grid, algorithm_parameters, verbose = F))
      nested_errors <- array(NA, dim = c(length(fold_list), length(fold_list), length(parameters_grid)))
      #nested_errors
      for(i in 1:length(fold_list)) {
        #subfold
        for(j in (1:length(fold_list))[-i]) {
          u1 = min(i,j); u2 = max(i,j)
          nested_errors[i,j, ] = test_classifier(nested_cross_validation[[intersect(which(combination_folds[1,]==u1),which(combination_folds[2,]==u2))]], 
                                                 X[fold_list[[j]], ], Y[fold_list[[j]]])
        }
      }
      cat("Nested over!!\n")
    }
    cross_validation <- lapply(1:length(fold_list), function(i)
      train_classifier(X[-fold_list[[i]],], Y[-fold_list[[i]]],
                       parameters_grid, algorithm_parameters, verbose = F))
    cross_validation_errors <- sapply(1:length(fold_list), function(i) 
      test_classifier(cross_validation[[i]], X[fold_list[[i]], ], Y[fold_list[[i]]]))
  }
  ################################################################################
  # Construct cv matrix  ---------------------------------------------------------
  ################################################################################
  if(length(parameters_grid)==1) {
    cv_errors <- mean(cross_validation_errors)
    sd_errors <- sd(cross_validation_errors)
    errors_table <- cross_validation_errors
  }else{
    cv_errors <- rowMeans(cross_validation_errors)
    require(genefilter)
    sd_errors <- rowSds(cross_validation_errors)
    errors_table <- cross_validation_errors
    dimnames(errors_table) = list(parameters = sapply(1:length(parameters_grid), function(v) paste('param',v,sep="")),
                                  fold = sapply(1:length(fold_list), function(v) paste('f',v,sep="")))
    
  }
  if(sparsity_results) {
    sparsity_table <- sapply(cross_validation, function(cvs) cvs$sparsity)
    node_sparsity_table <- sapply(cross_validation, function(cvs) cvs$node_sparsity)
    dimnames(sparsity_table) = list(parameters = sapply(1:length(parameters_grid), function(v) paste('param',v,sep="")),
                                    fold = sapply(1:length(fold_list), function(v) paste('f',v,sep="")))
    dimnames(node_sparsity_table) = list(parameters = sapply(1:length(parameters_grid), function(v) paste('param',v,sep="")),
                                         fold = sapply(1:length(fold_list), function(v) paste('f',v,sep="")))
    results <- list(methodname = methodname, errors = errors_table, cv_errors = cv_errors,
                    sd_errors = sd_errors, parameters_grid = parameters_grid,
                    sparsity_table = sparsity_table, node_sparsity_table = node_sparsity_table)
  }else{
    results <- list(methodname = methodname, errors = errors_table, cv_errors = cv_errors,
                    sd_errors = sd_errors, parameters_grid = parameters_grid)  
  }
  
  if(nested_cv) {
    results$nested_errors = nested_errors
    results$nested_cverrors = apply(nested_errors, c(1,3), mean, na.rm = T)
  }
  return(results) 
    
}
