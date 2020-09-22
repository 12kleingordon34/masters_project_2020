library(bnlearn)
library(dplyr)
library(equSA)
library(pcalg)
library(reticulate)

use_condaenv("r-reticulate")
source_python('./notears/notears/linear.py')
source_python('./notears/notears/utils.py')


calculate.performance.statistics <- function(learned.graph, true.graph){
  #' Takes in two graph objects of 'bn' class
  adj.true <- amat(true.graph)
  adj.learned <- amat(learned.graph)
  
  comparison <- compare(cpdag(true.graph), cpdag(learned.graph))
  
  TP <- comparison$tp
  FP <- comparison$fp
  FN <- comparison$fn
  
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


convert.amat.to.bn <- function(adj.matrix){
  bn.g <- empty.graph(colnames(adj.matrix))
  amat(bn.g) <- (adj.matrix != 0) * 1
  return(bn.g)
}


hybrid.notears <- function(data, params){
  l1 = params$l1; alpha = params$alpha; test = params$test; B = params$B 
  adj.mat <- amat(hpc(as.data.frame(data), alpha=alpha, test=test, B=B))
  vectorised.bounds <- as.vector(t(adj.mat))
  learned.dag.amat <- notears_linear(
    as.matrix(data),
    lambda1=l1,
    loss_type='l2',
    bounds=vectorised.bounds
  )

  # Check boundaries have been applied correctly
  stopifnot(sum(learned.dag.amat * adj.mat) == sum(learned.dag.amat))
  colnames(learned.dag.amat) <- colnames(data)
  return(learned.dag.amat)
}


lingam.amat <- function(data){
  #' Learns a DAG structure from input data using
  #' LinGAM function, and returns the adjacency matrix.
  bn.g <- empty.graph(colnames(data))
  learned.lingam <- lingam(data)
  amat(bn.g) <- (learned.lingam$Bpruned != 0) * 1

  return(amat(bn.g))
}


bootstrap.alt <- function(algorithm, algo.params, data, bootstrap.size, num.bootstraps, seed=0){
  i <- 1
  shift <- 0
  model.list <- vector(num.bootstraps, mode='list')
  while (i <= num.bootstraps+shift){
    set.seed(i+seed)
    bootstrap.sample <- sample_n(data, bootstrap.size)
    if (algo.params$hybrid == TRUE){
      amat.dag <- algorithm(as.matrix(bootstrap.sample), algo.params)
    } else {
      l1 <- algo.params$l1; loss.type <- algo.params$loss.type
      amat.dag <- algorithm(as.matrix(bootstrap.sample))
    }
    # Check if DAG is acyclic
    if (!is_dag(amat.dag)){
      shift <- shift + 1
      print(sprintf("Cycles detected. Shift set to %s", shift))
      break
    }
    colnames(amat.dag) <- colnames(data)
    learned.dag <- convert.amat.to.bn(amat.dag)
    model.list[[i]] <- learned.dag
    i <- i + 1
  }#FOR
  
  # Recompile model arcs
  arclist = list()
  for (i in seq_along(model.list)) {
    # extract the models
    run = model.list[[i]]
    arclist[[length(arclist) + 1]] = arcs(run)
  }#FOR
  
  # compute the arc strengths.
  nodes = unique(unlist(arclist))
  strength = custom.strength(arclist, nodes = nodes)
  return(strength)
}