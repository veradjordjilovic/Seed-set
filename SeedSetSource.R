### functions for the estimation of the graphical seed set


EqDist.test <- function(Sigma, Sigma1, Sigma2, n1,  n2){
  # var1... estimate of the covariance matrix in the first sample
  # n1 ... size of the first sample
  # var2 ... estimate of the covariance matrix in the second sample
  # n2... size of the second sample
  Sigma <- as.matrix(Sigma)
  Sigma1 <- as.matrix(Sigma1)
  Sigma2 <- as.matrix(Sigma2)
  p <- nrow(Sigma)
  df <- p*(p+3)/2
  TS <- n1*log(det(Sigma)/det(Sigma1))+n2*log(det(Sigma)/det(Sigma2))
  p_value <- pchisq(TS, df = df, lower.tail = F )
  return(list("TS" = TS, "p-value" = p_value))
}

EqDist<-function(index, Sigma, Sigma1, Sigma2, n1, n2) 
  EqDist.test(Sigma[index, index],
              Sigma1[index,index],
              Sigma2[index,index], n1, n2)

EqDistGraph <- function(graph, sample1, sample2, root = NULL){
  # function for testing equality of two  normal distributions 
  # Markov to the graph 
  # sample1... a random sample from the first distribution
  # sample2... a random sample from the second distribution
  ripped<-gRbase::rip(graph, root = root)
  indexClique_list<- ripped$cliques
  indexSep_list<- (ripped$sep)[-1]
  
  Set_list <- append(indexClique_list, indexSep_list)
  
  n1 <- nrow(sample1)
  n2 <- nrow(sample2)
  
  Sigma <- (n1+n2-1)/(n1+n2)*var(rbind(sample1,sample2))
  Sigma1 <- (n1-1)/n1*var(sample1)
  Sigma2 <- (n2-1)/n2*var(sample2)
  
  
  TestClique<-lapply(indexClique_list, EqDist,
                          Sigma, Sigma1, Sigma2, n1, n2)
  
  TestSep<-lapply(indexSep_list, EqDist,
                          Sigma, Sigma1, Sigma2, n1, n2)
  Set_listR <- list(SetName=Set_list, 
                         TestST=c(sapply(TestClique, "[[",1),
                                  sapply(TestSep, "[[",1)),
                         pval= c(sapply(TestClique, "[[",2),
                                 sapply(TestSep, "[[",2)))
  return(Set_listR)
}
    


Components <- function (graph, root = NULL, results){
  ripped<-gRbase::rip(graph, root = root)
  indexClique_list<- ripped$clique
  indexSep_list<- (ripped$sep)[-1]
  comNames <- CliqSepNames(ripped)$namesD
   components <- list()
   ind <-  which(sapply(results$SetName, function(x)
   setequal(x,indexClique_list[[1]])))
   TS <- numeric()
   pval <- numeric()
   TS[1]<- results$TestST[ind]
   pval[1] <-results$pval[ind]
  for (i in (2:length(indexClique_list))){
        p1 <- length(indexClique_list[[i]])
        p2 <- length(indexSep_list[[i-1]])
        df <- p1*(p1+3)/2 - p2*(p2+3)/2
        t1 <- which(sapply(results$SetName, function(x)
                    setequal(x,indexClique_list[[i]])))
        t2 <- min( which(sapply(results$SetName, function(x)
                    setequal(x,indexSep_list[[i-1]]))))
        TS[i] <- results$TestST[t1] - results$TestST[t2]
        pval[i] <-  pchisq(TS[i], df, lower.tail = F) 
  }
   order_comp <- order_Components(graph =  graph, root1 = NULL, root2 = root)
   # clique numbering with respect to the default junction tree
   components  <- list("ComponentName"= comNames, 
                      "TestST" = TS,
                      "pval" = pval,
                      "clique_number"= order_comp
                      )
}

CliqSepNames <- function(jt){
  namesC <- list()
  namesC[[1]] <- paste(sort(jt$cliques[[1]]), collapse = ",")
  namesS <- list()
  namesS[[1]] <- " "
  namesD <- list()
  namesD[[1]] <- namesC[[1]]
  for (i in 2:length(jt$cliques)){
    namesC[i] <- paste(sort(jt$cliques[[i]]), collapse = ",")
    namesS[i] <- paste(sort(jt$separators[[i]]), collapse = ",")
    namesD[i] <- paste(sort(setdiff(jt$cliques[[i]],
                                    jt$separators[[i]])), collapse=",")
    namesD[i] <- paste("(",namesD[i],") | (",namesS[i],")")
  }
  list("namesC"=namesC, "namesS"=namesS, "namesD"=namesD)}


order_Components <- function(graph = graphD, root1 = NULL, root2=NULL){
  jt1 <- gRbase::rip(graph, root = root1)
  jt2 <- gRbase::rip(graph, root = root2)
  cliques1 <- jt1$cliques
  cliques2 <- jt2$cliques
  which_clique <- function(cl) which(sapply(cliques1, function(x)
    setequal(x,cl)))
  order <- sapply(cliques2, which_clique)
  order
}

treshold_Components <- function(results_components = components_chimera,
                                  treshold = 0.05){
  # function that returns indices of significant components 
  # with respect to the default junction tree <=> root = NULL
  # significant means with p_value not exceeding the treshold
ind <- which(results_components$pval<=treshold)
sig_cliques <- (results_components$clique_number)[ind]
  }


  
### function that passes from mixed parametrization to the canonical 

fromMixedToCan <- function(mu1, Sigma1, w1, w2, w3){
  # Y = (Y_1, Y_2) ... decomposition
  # mu1, Sigma1 ... paramteres of Y_1
  # w1, w2, w3 ... parameters of Y_2|Y_1
  p1 <- length(mu1)
  p2 <- length(w2)
  #w2 <- matrix(data = w2, ncol = p1, nrow = p2, byrow = T)
  Sigma21  <- w2 %*% Sigma1
  Sigma2 <- w3 + w2 %*% Sigma1 %*% t(w2)
  mu2 <- w1 + w2%*%mu1
  mu <- c(mu1, mu2)
  Sigma <- rbind(cbind(Sigma1, t(Sigma21)), cbind(Sigma21, Sigma2))
  return(list("mean" = mu, "Sigma"= Sigma))
}


fromCanToMixed <- function(mu, Sigma, resp, cond){
  #resp ... response variables index
  #cond ... conditioning variables index
  mu1 <- mu[cond]
  Sigma1 <- Sigma[cond, cond]
  mu2 <- mu[resp]
  Sigma2 <- Sigma[resp, resp]
  Sigma12 <- matrix(data=Sigma[cond, resp], nrow=length(cond),
                    ncol=length(resp))
  w1 <- mu2 - t(Sigma12) %*% solve(Sigma1) %*% 
    mu1
  w2 <- t(Sigma12) %*% solve(Sigma1)
  w3 <- Sigma2 - t(Sigma12) %*% solve(Sigma1) %*% Sigma12
  return(list ("w1"=w1, "w2" = w2, "w3"= w3))
}

simPar <- function(A, ss=2, rmean=c(-2,2), rbeta=c(0.2,1),
                   rsigma=c(0.5,1)){
  # A... adjacency matrix of a DAG
  # ss... the value of the seed for set.seed()
  # rmean... interval range for mean values
  # rbeta... interval range for paramteres of regressions
  # rsigma... interval range for the residual variances
  p <- ncol(A)
  nb <- sum(A==1)
  set.seed(ss)
  
  B <- matrix(0, nrow = p, ncol=p,
              dimnames=list(c(1:p), c(1:p))) # matrix of beta coefficients
  sign <- c(-1,1)
  indnz <- which(A==1, arr.ind=T)
  Rsign <- sample(sign, size=nb, replace=T)
  beta <- runif(nb, min=rbeta[1], max=rbeta[2])
  beta <- Rsign*beta
  B[indnz] <- beta
  
  sigma2 <- runif(p, min=rsigma[1], max=rsigma[2]) # the vector of residual variances
  
  m <- runif(p, min=rmean[1], max=rmean[2])  # intercepts (residual mean)
  
  L <- solve(diag(p)-B)
  mu <- m%*%L
  Sigma <- t(L)%*%diag(sigma2)%*%L # the mean and the variance
  return(list("mean"=mu, "Sigma"=Sigma, "Beta"=B, "L"=L))
}

EqDistGraphGlobal <- function(graph, sample1, sample2){
  # function for testing equality of two  normal distributions 
  # Markov to the graph 
  # sample1... a random sample from the first distribution
  # sample2... a random sample from the second distribution
  ripped<-gRbase::rip(graph)
  indexClique_list<- ripped$cliques
  indexSep_list<- (ripped$sep)[-1]
  
  Set_list <- append(indexClique_list, indexSep_list)
  
  n1 <- nrow(sample1)
  n2 <- nrow(sample2)
  
  Sigma <- (n1+n2-1)/(n1+n2)*var(rbind(sample1,sample2))
  Sigma1 <- (n1-1)/n1*var(sample1)
  Sigma2 <- (n2-1)/n2*var(sample2)
  
  
  TestClique<-lapply(indexClique_list, EqDist,
                     Sigma, Sigma1, Sigma2, n1, n2)
  
  TestSep<-lapply(indexSep_list, EqDist,
                  Sigma, Sigma1, Sigma2, n1, n2)
  TS <- sum(unlist(TestClique))-sum(unlist(TestSep))
  df <- 2*ncol(sample1)+0.5*sum(as(graph,"matrix"))
  pval <- pchisq(TS, df, lower.tail = F)                           
  return(pval=pval)
}





minP <- function(sample1,sample2, sig_lev = c(0.01,0.03,0.05,0.1),T=2000, graph){
  sig_q <- T*sig_lev
  results <- EqDistGraph(graph, sample1, sample2)
  Comp_complete <- lapply(cliques, tf, resultsG=results, graph=graph)
  Comp_complete_dataf <- cbind.data.frame("ComponentName" = unlist(lapply( Comp_complete, function(x)x[,1])),
                                          "TS"= unlist(lapply( Comp_complete, function(x)x[,2])))
  Comp_complete_dataf <- unique(Comp_complete_dataf)
  ts <- Comp_complete_dataf[,2]
  resST <-matrix(0,ncol=dim(Comp_complete_dataf)[1],nrow=T)
  resST[1,] <-ts
  sampleC <- rbind(sample1,sample2)
  n <- nrow(sample1)
  for (i in 2:T){
    print(i)
    ind1 <- sample(2*n,n, replace=F)
    sample1p <-sampleC[ind1,] 
    sample2p <- sampleC[-ind1,]  
    resultsp <- EqDistGraph(graph, sample1p, sample2p)
    Comp_complete <- lapply(cliques, tf, resultsG=resultsp, graph=graph)
    Comp_complete_dataf <- cbind.data.frame("ComponentName" = unlist(lapply( Comp_complete, function(x)x[,1])),
                                            "TS"= unlist(lapply( Comp_complete, function(x)x[,2])))
    
    Comp_complete_dataf <- unique(Comp_complete_dataf)
    resST[i,] <-  Comp_complete_dataf[,2]
  }
 
 
  resTm = apply(-resST, 2, rank, ties.method = "max", na.last = "keep")/T
  resTm = as.matrix(resTm)
  colnames(resTm)<- 1:dim(resTm)[2]
  mins <- apply(resTm,1,min)
  thresh <- sort(mins)[sig_q]
  rejected <- colnames(resTm)[which(resTm[1,]<thresh)]
  rt <- rejected
  resTmt <- resTm
  while ((length(rt)>0) && (length(rejected)<dim(resTm)[2])){
    resTmt <- as.matrix(resTmt[,-which(colnames(resTmt) %in% rt)])
    mins <- apply(resTmt,1,min)
    thresh <- sort(mins)[sig_q]
    rt <- colnames(resTmt)[which(resTmt[1,]<thresh)]
    rejected <- union(rejected,rt)
  }
  rejectedV <- rep(0, dim(resTm)[2])
  rejectedV[as.numeric(rejected)] <-1
  return(list(rejectedV = rejectedV, thresh=thresh, permM = resTm, resST = resST))
}


sf <- function(cl,  resultsG, graphG, threshold=0.05){
  t <-  Components( graphG, root = cl, results= resultsG)
  t_res <- treshold_Components(results_components=t, treshold = threshold)
  t_res
} 

vf <- function(cl, resultsG, graphG){
  t <-  Components(graphG, root = cl, results= resultsG)
  cbind("ComponentName"=t$ComponentName,"p-value"= t$pval)
}

tf <- function(cl, resultsG, graphG){
  t <-  Components(graphG, root = cl, results= resultsG)
  cbind("ComponentName"=t$ComponentName,"TestST"= t$TestST)
}

Components2 <- function (graph, root=NULL){
  ripped<-gRbase::rip(graph, root = root)
  indexClique_list<- ripped$clique
  indexSep_list<- (ripped$sep)[-1]
  comNames <- CliqSepNames(ripped)$namesD
  components <- list()
  # ind <-  which(sapply(results$SetName, function(x)
  #   setequal(x,indexClique_list[[1]])))
  order_comp <- order_Components(graph =  graph, root1 = NULL, root2 = root)
  # clique numbering with respect to the default junction tree
  components  <- list("ComponentNames"= comNames, 
                      "clique_number"= order_comp)
}

Seed_set <- function(rejectedComp,graphG){
  jt <- gRbase::rip(graphG)
  cliques <- jt$cliques
  Comp_complete <- lapply(cliques, Components2, graph=graphG)
  Comp_complete_dataf <- cbind.data.frame("ComponentName" = unlist(lapply( Comp_complete, function(x)x[[1]])),
                                          
                                          "clique"= unlist(lapply( Comp_complete, function(x)x[[2]])))
  # select unique tests
  UniqueComp <- subset(Comp_complete_dataf, !duplicated(Comp_complete_dataf$ComponentName))
  names(rejectedComp) <-  UniqueComp[,1]
  sig_comp <- lapply(Comp_complete, 
                     function(x) x$clique_number[which(rejectedComp[unlist(x$ComponentNames)]==1)])
  un_cl <-lapply(sig_comp, function(x) unique(unlist(cliques[x])))
  temp <- Reduce(intersect, un_cl)
  if (is.null(temp)| (length(temp)==0)) {
    seed_set <-list(NULL)}else{seed_set<-sort(temp)}
  return(seed_set)}

