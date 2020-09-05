# Simulation study 2

require(igraph)
require(gRbase)
require(mvtnorm)
require(SourceSet)
require(ggm)
require(graph)
require(simPATHy)


# Uncomment the following lines for comaprison with the method of Zhao et al.(2014)
# dyn.load("dpm.so");
# library(Rcpp)
# source("dpm.R");
# library(MASS);



# specify the graph
graph <- gRbase::ug(~ 1:2:3:4 + 4:5:6 + 4:6:7 + 7:8+6:9:10)
nodes <- graph::nodes(graph) 
p <- length(nodes)
plot(graph)
cliques <- gRbase::getCliques(graph) 

# specify the parameters of the first condition
set.seed(1357)
mu1 <- rnorm(p, mean = 0.5)
names(mu1) <- nodes
S1<- matrix(0.6, ncol=p,nrow=p) + diag(0.4, ncol = p, nrow = p)
colnames(S1) <- rownames(S1) <- nodes
S1<- simPATHy::fitSgraph(graph, S1)

# for the second condition, remove the link between 4 and 6
K1 <- solve(S1)
K2 <- K1
K2["4","6"] <- K2["6", "4"] <-  0
S2 <- solve(K2)



# simulation study
B <- 20 # number of Monte Carlo runs, keep low for illustration purposed
n <- 200 # number of observations


# prepare the objects for the results
gsset <- matrix(0, nrow = B, ncol = length(mu1))
colnames(gsset) <- nodes(graph)
gllist <- list()
dne <- list()
set.seed(137)

for (i in (1:B)) {
  data.class1 <- rmvnorm(n, mean = rep(0,p), sigma = S1)
  data.class2 <- rmvnorm(n, mean = rep(0,p), sigma = S2)
  data<-rbind(data.class1,data.class2)
  colnames(data)<-graph::nodes(graph)
  classes<-c(rep(1,nrow(data.class1)),rep(2,nrow(data.class2)))
  ripped<-SourceSet::ripAllRootsClique(graph)
  res1<-SourceSet:::singleSourceSet(ordering = ripped, data = data, classes = classes,
                                    permute = FALSE, return.permutations = FALSE,
                                    seed = 123)
  gsset[i, res1$primarySet] <- 1
  # uncomment the following line for the method of Zhao et al.(2014)
  # dne[[i]]<-  dpm(data.class1, data.class2, nlambda=10, tuning="cv")
}

# Proprtion of runs in which each of the nodes belongs to the graphical seed set estimate
colMeans(gsset)

# If dna has been comupted, here is the comparison
# optimalfindings <- sapply(dne, function(x) (x$dpm[x$opt[1]]))
# test <- sapply(optimalfindings, function(x) colSums(x) != 0)
# rbind(colMeans(gsset), rowMeans(test))
