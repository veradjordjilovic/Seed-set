# Simulation study 1

require(igraph)
require(gRbase)
require(mvtnorm)
require(SourceSet)
source("SeedSetSource.R")

# generate a random graph
p<-100
seed<-123
appeal <- 20
power <- 2.5
directed <- TRUE

set.seed(seed)
ig<-igraph::barabasi.game(p,zero.appeal = appeal,power = power,directed = directed)
graph<-igraph::igraph.to.graphNEL(ig)
PaperGraph<-list(directed=graph,
                 triangulated= gRbase::triangulate(gRbase::moralize(graph)), 
                 parameters=list(seed=seed,zero.appeal=appeal,power=power,directed=directed,call="barabasi.game(100,zero.appeal = appeal,power = power,directed = directed)"))

max(sapply(rip(PaperGraph$triangulated)$cliques,length))
nodes <- nodes(PaperGraph$triangulated)
graph <- PaperGraph$triangulated

# generate parameters of the control condition
set.seed(1357)
mu1 <- rnorm(p,mean=0.5)
names(mu1) <- nodes
S1<- matrix(0.4, ncol=p,nrow=p)+diag(0.6, ncol=p,nrow=p)
colnames(S1) <- rownames(S1) <- nodes
S1<- simPATHy::fitSgraph(graph,S1)
# uncomment to save the parameters
#saveRDS(mu1, file = "SimResults/mu1.rds")
#saveRDS(S1, file = "SimResults/S1.rds")


# set the seed set
seedset <- c("2","5")
gseedset <- c("2", "5", "17")

# plot the graph
nAttrs<-list()
cn <- c(rep("red", 2), rep("grey", 98))
names(cn) <- c(seedset, setdiff(nodes(graph), seedset))
nAttrs$fillcolor <- cn
plot(graph, nodeAttrs = nAttrs)

# lambda determines the strength of the perturbation (uncomment the chosen line)

  lambda <- 1.1   # mild perturbation
# lambda <- 1.3   # medium perturbation
# lambda <- 1.7   # strong perturbation

# generate parameters of the perturbed condition
resp <- setdiff(nodes(graph),seedset)
param_t <- fromCanToMixed(mu1,S1, resp=resp, cond=seedset)
ordering <- order(as.numeric(c(seedset,resp)))
param_t2 <- mu1[seedset]
param_t2[seedset] <- lambda*param_t2[seedset]
param_t3 <- as.matrix(S1[seedset,seedset])
var_change <- rep(0.5,length(seedset)) #decrease variance by 50%
param_t3 <-diag(sqrt(var_change))%*%cov2cor(param_t3)%*%diag(sqrt(var_change))

new_param <- fromMixedToCan(param_t2, param_t3,
                            param_t$w1, param_t$w2, param_t$w3)
mu2 <- new_param$mean
mu2<-mu2[ordering]
names(mu2) <- names(mu1)


S2 <- new_param$Sigma
colnames(S2) <- rownames(S2) <- c(seedset,resp)
S2 <- S2[ordering,ordering]


# saveRDS(mu2, file = "SimResults/mu2MildInt.rds") 
# saveRDS(S2, file = "SimResults/S2MildInt.rds") 
# saveRDS(mu2, file ="SimResults/mu2MediumInt.rds") 
# saveRDS(S2, file ="SimResults/S2MediumInt.rds") 
# saveRDS(mu2, file ="SimResults/mu2StrongInt.rds") 
# saveRDS(S2, file ="SimResults/S2StrongInt.rds")


SimSeedSet <- function (n, mu1,S1,mu2,S2, B = 500, graph) {
  # Wrapper function to perfrom the simulation study
  # Input:
  #   n... sample size 
  #   mu1... mean of the first condition
  #   S1... variance matrix of the first condition
  #   mu2, S2... mean and variance of the second condition
  #   B... number of Monte Carlo runs
  #   graph... decomposable undirected graph associated to S1 and S2
  #
  # Output:
  #   B x length(mu1) zero one matrix, where element (i,j) equals 1 if variable/node j belongs
  #   to the graphical seed set estimator computed in the i-th simulated dataset
  
  gsset <- matrix(0, nrow = B, ncol = length(mu1))
  colnames(gsset) <- nodes(graph)
  for (i in (1:B)) {
    data.class1 <- rmvnorm(n, mean = mu1, sigma = S1)
    data.class2 <- rmvnorm(n, mean = mu2, sigma = S2)
    data<-rbind(data.class1,data.class2)
    colnames(data)<-graph::nodes(graph)
    classes<-c(rep(1,nrow(data.class1)),rep(2,nrow(data.class2)))
    ripped<-SourceSet::ripAllRootsClique(graph)
    res1<-SourceSet:::singleSourceSet(ordering = ripped, data = data, classes = classes,
                                      permute = FALSE, return.permutations = FALSE,
                                      seed = 123)
    gsset[i, res1$primarySet] <- 1
  }
  return(gsset)
} 


# Read in the parameters if previously saved

# mu1 <-readRDS("SimResults/mu1.rds")
# S1  <- readRDS("SimResults/S1.rds")
# mu2mild <-readRDS("SimResults/mu2MildInt.rds")
# mu2medium <-readRDS("SimResults/mu2MediumInt.rds")
# mu2strong <-readRDS("SimResults/mu2StrongInt.rds")
# S2  <- readRDS("SimResults/S2MildInt.rds")


# To perform the simulation study for different values of n, mu1, mu2, S1 and S2, 
# uncomment the desired line below.
# For a quick check, run the code with B = 10

# t<- proc.time()
# t1 <- SimSeedSet(n=50, mu1=mu1, S1=S1, mu2=mu2mild, S2=S2, B=1000, graph=graph)
# t2 <- SimSeedSet(n=75, mu1=mu1, S1=S1, mu2=mu2mild, S2=S2, B=1000, graph=graph)
# t3 <- SimSeedSet(n=100, mu1=mu1, S1=S1, mu2=mu2mild, S2=S2, B=1000, graph=graph)
# 
# t4 <- SimSeedSet(n=50, mu1=mu1, S1=S1, mu2=mu2medium, S2=S2, B=1000, graph=graph)
# t5 <- SimSeedSet(n=75, mu1=mu1, S1=S1, mu2=mu2medium, S2=S2, B=1000, graph=graph)
# t6 <- SimSeedSet(n=100, mu1=mu1, S1=S1, mu2=mu2medium, S2=S2, B=1000, graph=graph)
# 
# t7 <- SimSeedSet(n=50, mu1=mu1, S1=S1, mu2=mu2strong, S2=S2, B=1000, graph=graph)
# t8 <- SimSeedSet(n=75, mu1=mu1, S1=S1, mu2=mu2strong, S2=S2, B=1000, graph=graph)
# t9 <- SimSeedSet(n=100, mu1=mu1, S1=S1, mu2=mu2strong, S2=S2, B=1000, graph=graph)
# 
# t10 <- SimSeedSet(n=50, mu1=mu1, S1=S1, mu2=mu1, S2=S1, B=1000, graph=graph)
# t11 <- SimSeedSet(n=75, mu1=mu1, S1=S1, mu2=mu1, S2=S1, B=1000, graph=graph)
# t12<- SimSeedSet(n=100, mu1=mu1, S1=S1, mu2=mu1, S2=S1, B=1000, graph=graph)



# Analyze the results
tv <- rep(0,100)
names(tv) <- nodes(graph)
tv[gseedset] <- 1
tres <- list(t1,t2,t3,t4,t5,t6,t7,t8,t9)
sapply(tres, function(temp){ sum(apply(temp, 1, function(x)identical(x,tv)))})

# how many succesful recoveries of D_G under the null
tvundernull <-rep(0,100)
names(tvundernull) <- nodes(graph)
tundernull <- list(t10, t11, t12)
sapply(tundernull, function(temp){ sum(apply(temp, 1, function(x)identical(x,tvundernull)))})


#check FWER control in general
sapply(tres, function(temp){ sum(apply(temp, 1, function(x)sum(x[!(names(x)%in% gseedset)])>0) )})