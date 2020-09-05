# Chronic myeloid leukemia example

# load required packages
require(simPATHy)
require(gRbase)
require(graph)
require(graphite)
require(SourceSet)
source("SeedSetSource.R")


# prepare the pathway graph
graph <- graphite::pathwayGraph(pathways("hsapiens", "kegg")$'Chronic myeloid leukemia')
genes <- nodes(graph)
genes <- unlist(regmatches(genes, gregexpr("[[:digit:]]+", genes)))
plot(graph)
nodes(graph) <- genes

# This graph contains biderected and directed edges so we
# transform it into an undirected graph. 

# But first, prepare the data and the subgraph induced by the 
# measured genes pariticpating in the chosen pathway 
data(chimera)
position <- match(genes, rownames(chimera))
position <- na.omit(position)
data <- chimera[position, ]
data1 <- t(data[,colnames(data)=="1"])
data2 <- t(data[,colnames(data)=="2"])
genes <- rownames(data)
plot(graph)

# We focus on the largest connected component, so we eliminate the rest
genes <- genes[!is.element(genes, c("5595","7040","7043","7046","4088","5594",
                                    "1019","1021","1869","5925", "10912", "4616"))]
graph <- graph::subGraph(genes, graph)
plot(graph)

# Subset the data as well
data1 <- data1[, genes]
data2 <- data2[, genes]

# Moralize and triangulate to obtain a decomposable graph
# and find its junction tree
m <- gRbase::graphNEL2M(graph)
m <- gRbase::moralizeMAT(m)
graphU <-gRbase::coerceGraph(m, "graphNEL")
graphD <- gRbase::triangulate(graphU)
plot(graphD)  # this is the graph we will use for the analysis
jt<-gRbase::rip(gRbase::triangulate(graphD))

# Global test of equality 
n1 <- nrow(data1)
n2 <- nrow(data2)
Sigma <- (n1+n2-1)/(n1+n2)*var(rbind(data1,data2))
Sigma1 <- (n1-1)/n1*var(data1)
Sigma2 <- (n2-1)/n2*var(data2)
Sigmag <- fitSgraph(graphD, Sigma)
Sigma1g <- fitSgraph(graphD, Sigma1)
Sigma2g <- fitSgraph(graphD, Sigma2)
p <- nrow(Sigma)
df <- 2*p+sum(as(graphD, "matrix"))/2
TS <- n1*log(det(Sigmag)/det(Sigma1g))+n2*log(det(Sigmag)/det(Sigma2g))
(pval <- pchisq(TS, df, lower.tail =F))  # we reject the global hypothesis of equality. 

# Estimate the graphical seed set with the SourceSet package
# We use the function singleSourceSet sicne we consider a single graph. 
data <- rbind(data1, data2) # prepare the data
classes <-  c(rep("1", nrow(data1)), rep("2", nrow(data2))) # prepare the vector of classes
ripped <- SourceSet::ripAllRootsClique(graphD)  # the function computing all possible orderings of cliques

res <- SourceSet:::singleSourceSet(ordering = ripped, data = data, classes = classes,
                                  permute = TRUE, return.permutations = FALSE,
                                  seed = 123, shrink = FALSE)
# permute = TRUE indicates that  permutation p-values instead of the asymptotic p-values
#           should be computed
# shrink = FALSE computes maximum likelihood estimates of the parameters 

# The estimated graphical seed set
res$primarySet

# The set of variables that are marginally different in two conditions (in addition to the graphical seed set)
res$secondarySet

# p-values of all component hypotheses
res$Components

# Clique and separator identifiers are here.
res$Elements


# Results when asymptotic p-values are used (permute= FALSE)
res.alt<-SourceSet:::singleSourceSet(ordering = ripped, data = data, classes = classes,
                                 permute = FALSE, return.permutations = FALSE,
                                 seed = 123, shrink = FALSE)
# The graphical seed set now consists only of 25 and 613
res.alt$primarySet
# The third gene is now in the set of nodes that are marginally significant
res.alt$secondarySet



# Prepare the table with local tests results
listCl <- as.list(res$Components[,1])
l1 <- lapply(listCl, function(x)res$Elements[[x]])
listSep <- as.list(res$Components[,2])
l2 <- lapply(listSep, function(x)res$Elements[[x]])

ldiff <- mapply( function(x,y) setdiff(x,y), l1,l2)
components <- mapply(function(x,y) paste(paste(x,collapse=","),paste(y, collapse=","), sep="|"), ldiff,l2)
p.val <- cbind(components, round(res$Components$pvalue,4))
View(p.val)


