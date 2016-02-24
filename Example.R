# Generate population

# Simulated example from stochastic blockmodel
# Networks in population 2 have a different connectivity between communities 3 and 4
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


u = graphclass(X[train_population,],    # Matrix with network samples by row, with upper triangle of the adjacency
               Y = Y[train_population], # Binary response
               Xtest = X[test_population,],  #Test matrix
               Ytest = Y[test_population],   # Test response
               type = "intersection",               # One between "union" and "intersection"
               lambda1=2, lambda2=5,  # Lambda1 (Lasso parameter) and Lambda2 (group lasso parameter)
               params = NULL,               # List of parameters for optimization algorithm
                                            # If not provided, default parameters: MAX_ITER =300,  CONV_CRIT = 1e-05,   MAX_TIME = Inf
               id = "", verbose = T,       # verbose = T, then prints progress on every iteration, id identifies job when multiple jobs are running
               D = NULL)                  # For big networks, D can be constructed with D = construct_D(number of nodes), and provided each time function is called

u$train_error
u$test_error
source("Plots.R")
plot_adjmatrix(u$beta,
               type = "intersection") # Plot Adjacency matrix, depending whether is union or intersection solution

