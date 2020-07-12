calculate.performance.statistics <- function(true.graph, learned.graph){
  #' Takes in two graph objects of 'bn' class
  adj.true <- amat(true.graph)
  adj.learned <- amat(learned.graph)
  
  TP <- sum((adj.true == 1) * (adj.learned == 1))
  FP <- sum((adj.true == 0) * (adj.learned == 1))
  FN <- sum((adj.true == 1) * (adj.learned == 0))
  
  false.disc.rate <- FP/(FP + TP)
  sensitivity <- TP/(TP + FN)
  
  return(list(
    'fdr'=false.disc.rate,
    'sensitivity'=sensitivity
  ))
}


ebic <- function(bn.graph, data, gamma){
  #' Calculate EBIC criteria for a multivariate
  #' Gaussian graphical model. The `lambda` variable
  #' regularises the sparsity of the scoring criterion.
  ggm <- fitGGM(data=data, graph=amat(moral(bn.graph)), model='concentration')
  N <- dim(data)[1]; k <- dim(data)[2]
  
  # ggm$loglik returns postive values, probably defined as -log(L)
  bic <- 2 * ggm$loglik - log(N) * sum(amat(moral(bn.graph)))
  ebic <- bic - 4 * gamma * sum(amat(moral(bn.graph))) * log(k)
  return(ebic)
}


generate.dag <- function(num.data, num.nodes, sparsity, seed=0){
  #' Generates an observational dataset generated from
  #' a simulated multivariate Bayesian Network. Requires
  #' "DAGsim" function from "equSA" package
  #'
  #' param: num.data: Number of observations in dataset 
  #' param: num.nodes: Number of nodes in the generated DAG
  #'    (i.e. its dimension)
  #' param: sparsity: probability of a weight being active 
  #'    in the generated DAG
  #' 
  #' return: data
  set.seed(seed)
  sim.dag <- DAGsim(
    num.data,
    num.nodes,
    sparsity<-sparsity,
    p.binary <- 0,
    type<-"random"
  )
  
  adj.matrix <- sim.dag$edgematrix
  colnames(adj.matrix) <- paste("V", 1:num.nodes, sep='')
  rownames(adj.matrix) <- paste("V", 1:num.nodes, sep='')
  
  sim.data <- sim.dag$data
  colnames(sim.data) <- paste("V", 1:num.nodes, sep='')
  
  # Generate graph object using bnlearn package
  dag.graph <- empty.graph(paste("V", 1:num.nodes, sep=''))
  amat(dag.graph) <- adj.matrix
  
  return(list("graph" = dag.graph, "data" = sim.data))
}