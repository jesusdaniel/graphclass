# Generate population
source("Simulate_SBM.R")
communities = 6
NODES = 60
N1 = 200
N2 = 200
X1 = lapply(as.list(1:N1),function(k)sbm_q1(0.7,0.3,K = communities,nk = NODES/communities,0.3))
X2 = lapply(as.list(1:N2),function(k)sbm_q2(0.7,0.3,K = communities,nk = NODES/communities,0.5))
X_1 = sapply(X1,function(X) as.vector(X$Aobs[upper.tri(X$Aobs)]))
X_2 = sapply(X2,function(X) as.vector(X$Aobs[upper.tri(X$Aobs)]))
X = rbind(t(X_1),t(X_2))
Y = c(rep(-1,N1),rep(1,N2))


# Classifier
source("graphclass.R")
train_population = c(1:50, 201:250)
test_population = c(51:200, 251:400)

u = graphclass(X[train_population,], Y = Y[train_population],Xtest = X[test_population,],Ytest = Y[test_population],
               type = "union",
               lambda1=0.01, lambda2=0.05, params = NULL, id = "", verbose = T, D = NULL)

u$train_error
u$test
source("Plots.R")
plot_adjmatrix(u$beta,type = "union")

