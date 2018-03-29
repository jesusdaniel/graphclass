source("Cross-validation-function.R")
source("train-test-functions-cv.R")

#---- Load network dataset 
#   X matrix with rows for individuals and columns for edges,
#     taking the upper triangular part of the adjacency matrix)
#   Y class labels
#load("/data/jarroyor/data/COBRE/Data/COBRE_Xranks")
#load("../GRAPHCLASS_2/Data/COBRE_Xranks")
#load("../GRAPHCLASS_4/COBRE_Folds_1jun_1_10.RData")
#Standardize edges
Xnorm <- apply(Xranks, 2, function(v) (v-mean(v))/sd(v))

# load folds for cross validation (list with entries containing an array of indexes)
# load("/data/jarroyor/data/COBRE/Data/COBRE_Folds_1jun_1_10.RData")


#install.packages("penalizedSVM", repos='http://cran.us.r-project.org')
#library(penalizedSVM)

### SVM L1
parameters_grid <- lapply(10^seq(-1,4,length.out = 31), function(x) x)
svmL1_cv <- cross_validation_function(X = Xnorm, Y = factor(Y), parameters_grid = parameters_grid, 
                                      methodname = "penalizedSVM-L1",
                                      train_classifier = train_svml1, test_classifier = test_svml1,
                                      folds = 1, fold_list = folds_list, 
                                      parallel = F, num_clusters = 10, windows = T, 
                                      nested_cv = F, save_files = F, filename = "", sparsity_results=T) 

save(svmL1_cv, file="svmL1-cv.RData")



### Naive Bayes
parameters_grid <- seq(5, 200, length.out = 30)

naivebayes_cv <- cross_validation_function(X = Xnorm, Y = factor(Y), parameters_grid = parameters_grid, 
                                      methodname = "naivebayes",
                                      train_classifier = train_naivebayes, test_classifier = test_naivebayes,
                                      folds = 10, fold_list = folds_list, 
                                      parallel = T, num_clusters = 4, windows = T, 
                                      nested_cv = F, save_files = F, filename = "", sparsity_results=T) 
save(naivebayes_cv, file="naivebayes-cv4.RData")


### SVM L2
parameters_grid <- list(0.01,0.1, 1, 10, 100, 1000)
svme1071_cv <- cross_validation_function(X = Xnorm, Y = factor(Y), parameters_grid = parameters_grid, 
                                      methodname = "svm-e1071",algorithm_parameters = "linear",
                                      train_classifier = train_svml1, test_classifier = test_svml1,
                                      folds = 10, fold_list = folds_list, 
                                      parallel = T, num_clusters = 10, windows = T, 
                                      nested_cv = F, save_files = F, filename = "", sparsity_results=F) 

save(svme1071_cv, file="svme1071-cv.RData")

### Random forest
parameters_grid <- list(c(100, 500), c(200, 500), c(300, 500),
                        c(100, 1000), c(200, 1000), c(300, 1000), 
                        c(100, 2000), c(200, 2000), c(300, 2000))

rf1_cv <- cross_validation_function(X = Xnorm, Y = factor(Y), parameters_grid = parameters_grid, 
                                    methodname = "randomforest1",
                                    train_classifier = train_randomforest, test_classifier = test_randomforest,
                                    folds = 10, fold_list = folds_list, 
                                    parallel = T, num_clusters = 10, windows = T, 
                                    nested_cv = F, save_files = F, filename = "", sparsity_results=F) 

save(rf1_cv, file="rf1-cv.RData")

# With nested CV
source("CV/Cross-validation-function.R")
source("CV/train-test-functions-cv.R")

load("/data/jarroyor/data/COBRE/Data/COBRE_Xranks")

Xnorm <- apply(Xranks, 2, function(v) (v-mean(v))/sd(v))

load("/data/jarroyor/data/COBRE/Data/COBRE_Folds_1jun_1_10.RData")

#install.packages("penalizedSVM", repos='http://cran.us.r-project.org')
#library(penalizedSVM)



### GC
rho_seq <- c(10^(seq(2.5, -2, length.out = 31)))
grid = data.frame(rbind(1e-4, rho_seq, 1e-5))
parameters_grid <- lapply(grid, function(x) x)

NODES <- (1+sqrt(1+8*ncol(Xnorm)))/2
D <- construct_D(NODES)
gc1_cv <- cross_validation_function(X = Xnorm, Y = factor(Y), parameters_grid = parameters_grid, 
                                    methodname = "gc1",
                                    train_classifier = train_graphclass, test_classifier = test_graphclass,
                                    folds = 10, fold_list = folds_list, algorithm_parameters = list(D = D),
                                    parallel = T, num_clusters = 10, windows = T, 
                                    nested_cv = F, save_files = F, filename = "", sparsity_results=T) 

save(gc1_cv, file="gc1-cv.RData")

# GC2

rho_seq <- c(10^(seq(2.5, -2, length.out = 31)))
grid = data.frame(rbind(1e-7, rho_seq, 1e-5))
parameters_grid <- lapply(grid, function(x) x)

NODES <- (1+sqrt(1+8*ncol(Xnorm)))/2
D <- construct_D(NODES)
gc2_cv <- cross_validation_function(X = Xnorm, Y = factor(Y), parameters_grid = parameters_grid, 
                                    methodname = "gc2",
                                    train_classifier = train_graphclass, test_classifier = test_graphclass,
                                    folds = 10, fold_list = folds_list, algorithm_parameters = list(D = D),
                                    parallel = T, num_clusters = 10, windows = T, 
                                    nested_cv = F, save_files = F, filename = "", sparsity_results=T) 

gc2_cv = gc1_cv
save(gc2_cv, file="gc2-cv.RData")




### Relaxed GC
rho_seq <- c(10^(seq(2.5, -2, length.out = 31)))[1:15]
grid = data.frame(rbind(1e-4, rho_seq, 1e-5))
parameters_grid <- lapply(grid, function(x) x)

NODES <- (1+sqrt(1+8*ncol(Xnorm)))/2
D <- construct_D(NODES)
relaxed_gc1_cv <- cross_validation_function(X = Xnorm, Y = factor(Y), parameters_grid = parameters_grid, 
                                            methodname = "relaxed_gc1",
                                            train_classifier = train_relaxed_graphclass, test_classifier = test_relaxed_graphclass,
                                            folds = 10, fold_list = folds_list, algorithm_parameters = list(D = D),
                                            parallel = T, num_clusters = 10, windows = T, 
                                            nested_cv = F, save_files = F, filename = "", sparsity_results=T)

save(relaxed_gc1_cv = gc1_cv, file="relaxed_gc1-cv.RData")




### glmnet
lambda_seq <- 10^seq(from = 0, to = -7, length.out = 30)
lambda_seq <- lapply(lambda_seq, function (x) x)
glmnet2_cv <- cross_validation_function(X = Xnorm, Y = factor(Y), parameters_grid = lambda_seq, 
                                            methodname = "glmnet0.2",
                                            train_classifier = train_glmnet, test_classifier = test_glmnet,
                                            folds = 10, fold_list = folds_list, algorithm_parameters = list(alpha = 0.2),
                                            parallel = T, num_clusters = 2, windows = T, 
                                            nested_cv = F, save_files = F, filename = "", sparsity_results=T)

save(glmnet2_cv = glmnet2_cv, file="glmnet2_cv-cv.RData")
