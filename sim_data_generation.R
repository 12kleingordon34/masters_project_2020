library(bnlear)
library(equSA)

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


run.simulation.test <- function(N, S, D.range, sparsity, function.list, 
                                alpha.range, tabu.range){
  #' Function which 
  #'
  #'
  num.algos <- length(function.list)
  
  mean.shd.array <- array(numeric(), c(length(D.range), num.algos + 2))
  sd.shd.array <- array(numeric(), c(length(D.range), num.algos + 2))
  best.param.array <- array(numeric(), c(length(D.range), num.algos+1))  
  mean.run.time <- array(numeric(), c(length(D.range), num.algos+1))  
  
    
  colnames(mean.shd.array) <- c(c("N", "D"), names(function.list))
  colnames(sd.shd.array) <- c(c("N", "D"), names(function.list))
  colnames(best.param.array) <- c("D", names(function.list))
  colnames(mean.run.time) <- c("D", names(function.list))
    
  d.num <- 1
  for (d in D.range){
    mean.shd.array[d.num, 1] <- N; sd.shd.array[d.num, 1] <- N
    mean.shd.array[d.num, 2] <- d; sd.shd.array[d.num, 2] <- d
    best.param.array[d.num, 1] <- d; mean.run.time[d.num, 1] <- d

    for (algo in 1:num.algos){
      algo.name = names(function.list)[algo]
      if (substr(algo.name, 1, 4) == "tabu"){
        best.hyperparam <- 1e9; best.shd.mean <- 1e9; best.shd.std <- 1e9
        for (tabu.frac in tabu.range){
          # Multiply tabu frac by number of nodes
          tabu <- ceiling(d * tabu.frac)
          shd.vals <- c()
          run.times <- c()
          for (s in 1:S){
            sim.dag <- generate.dag(N, d, sparsity, s)
            sim.data.df <- as.data.frame(sim.dag$data)
            start.time <- Sys.time()
            print(sprintf(
              'D: %s. Sparsity: %s. Seed: %s. Algo: %s. Tabu Val: %s. Time: %s',
              d, sparsity,s,names(function.list)[algo], tabu, start.time
            ))
            tabu.dag <- function.list[[algo]](
              sim.data.df, 
              tabu=tabu, 
              max.tabu=tabu,
            )
            run.times <- c(run.times, Sys.time() - start.time)
            shd <- shd(cpdag(tabu.dag), cpdag(sim.dag$graph))
            shd.vals <- c(c(shd.vals), shd)
          }
          if (mean(shd.vals) < best.shd.mean) {
            best.shd.mean <- mean(shd.vals)
            best.shd.std <- sd(shd.vals)
            best.hyperparam <- tabu.frac
            av.run.time <- mean(run.times)
          }
        }
      } else {
        best.hyperparam <- 1e9; best.shd.mean <- 1e9; best.shd.std <- 1e9
        for (a in alpha.range) {
          shd.vals <- c()
          run.times <- c()
          for (s in 1:S){
            sim.dag <- generate.dag(N, d, sparsity, s)
            sim.data.df <- as.data.frame(sim.dag$data)    
            start.time <- Sys.time()
            print(sprintf(
              'D: %s. Sparsity: %s. Seed: %s. Algo: %s. Alpha Val: %s. Time: %s',
              d, sparsity,s,names(function.list)[algo], a, start.time
            ))                   
            alpha.dag <- function.list[[algo]](sim.data.df, alpha=a)  
            run.times <- c(run.times, Sys.time() - start.time)
            shd <- shd(cpdag(alpha.dag), cpdag(sim.dag$graph))
            shd.vals <- c(c(shd.vals), shd)
          }
          if (mean(shd.vals) < best.shd.mean) {
            best.shd.mean <- mean(shd.vals)
            best.shd.std <- sd(shd.vals)
            best.hyperparam <- a
            av.run.time <- mean(run.times)
          }
        }
      }
        mean.shd.array[d.num, (algo+2)] <- best.shd.mean
        sd.shd.array[d.num, (algo+2)] <- best.shd.std
        best.param.array[d.num, (algo+1)] <- best.hyperparam
        mean.run.time[d.num, (algo+1)] <- av.run.time
      }
    d.num <- d.num + 1
    } 
  return(list('mean'=mean.shd.array, 'std'=sd.shd.array, 'best.params'=best.param.array, "times"=mean.run.time))
}


tabu.xval.test <- function(N, S, D.range, true.sparsity, hyperparam.sparsity, constraint.algo){
  mean.shd.array <- array(numeric(), c(length(D.range), num.algos + 2))
  sd.shd.array <- array(numeric(), c(length(D.range), num.algos + 2))
  best.param.array <- array(numeric(), c(length(D.range), num.algos+1))  
  mean.run.time <- array(numeric(), c(length(D.range), num.algos+1))  
  
  
  colnames(mean.shd.array) <- c(c("N", "D"), c("tabu.null", "tabu.init", "tabu.constraint"))
  colnames(mean.shd.array) <- c(c("N", "D"), c("tabu.null", "tabu.init", "tabu.constraint"))
  colnames(best.param.array) <- c("D", c("tabu.null", "tabu.init.tabu", "tabu.init.sparsity", "tabu.constraint"))
  colnames(mean.run.time) <- c("D", c("tabu.null", "tabu.init", "tabu.constraint"))
  
  d.num <- 1
  for (d in D.range){
    mean.shd.array[d.num, 1] <- N; sd.shd.array[d.num, 1] <- N
    mean.shd.array[d.num, 2] <- d; sd.shd.array[d.num, 2] <- d
    best.param.array[d.num, 1] <- d; mean.run.time[d.num, 1] <- d
    
    # Run for tabu.null
    tabu.null <- tabu.hyperparam.iteration(N, d, S, sparsity, tabu.range, start=NULL, constraint.algo=NULL)
    mean.shd.array[d.num, 2] <- tabu.null$mean
    sd.shd.array[d.num, 2] <- tabu.null$std
    best.param.array[d.num, 2] <- tabu.null$hyperparam
    mean.run.time[d.num, 2] <- tabu.null$av.run.time  
    
    # Run for tabu.random
    tabu.random <- tabu.hyperparam.iteration(N, d, S, sparsity, tabu.range, start="random", constraint.algo=NULL)
    mean.shd.array[d.num, 3] <- tabu.random$mean
    sd.shd.array[d.num, 3] <- tabu.random$std
    best.param.array[d.num, 3] <- tabu.random$hyperparam
    best.param.array[d.num, 4] <- 0.2
    mean.run.time[d.num, 3] <- tabu.random$av.run.time  
    
    # Run for tabu.constraint
    tabu.constraint <- tabu.hyperparam.iteration(N, d, S, sparsity, tabu.range, start=NULL, constraint.algo=constraint.algo)
    mean.shd.array[d.num, 4] <- tabu.constraint$mean
    sd.shd.array[d.num, 4] <- tabu.constraint$std
    best.param.array[d.num, 5] <- tabu.constraint$hyperparam
    mean.run.time[d.num, 4] <- tabu.constraint$av.run.time  
    
    d.num <- d.num + 1
  } 
  return(list('mean'=mean.shd.array, 'std'=sd.shd.array, 'best.params'=best.param.array, "times"=mean.run.time))
}

tabu.learning.iteration <- function(N, d, s, sparsity, tabu, start=NULL, constraint.algo=NULL){
  sim.dag <- generate.dag(N, d, sparsity, s)
  sim.data.df <- as.data.frame(sim.dag$data)
  start.time <- Sys.time()
  print(sprintf(
    'D: %s. Sparsity: %s. Seed: %s. start: %s. constraint: %s Tabu Val: %s. Time: %s',
    d, sparsity, s, start, constraint.algo, tabu, start.time
  ))
  if (start='random'){
    set.seed(s)
    start <- randomDAG(d, 0.2, V=nodes(sim.dag$graph))
  } else {
    if (!is.null(constraint.algo)){
      start <- constraint.algo(sim.data.df)
    }
    tabu.dag <- function.list[[algo]](
      sim.data.df, 
      tabu=tabu, 
      max.tabu=tabu,
      start=start
    )
  }
  run.time <- Sys.time() - start.time
  shd <- shd(cpdag(tabu.dag), cpdag(sim.dag$graph))
  return("shd"=shd, "run.time"=run.time)
}


tabu.hyperparam.iteration <- function(N, d, S, tabu.range, start, constraint.algo=NULL){
  #'
  best.shd.mean <- 1e9
  for (tabu.frac in tabu.range){
    # Multiply tabu frac by number of nodes
    tabu <- ceiling(d * tabu.frac)
    shd.vals <- c()
    run.times <- c()
    for (s in 1:S){
      tabu.results <- tabu.learning.iteration(N, d, s, sparsity, start, constraint.algo)
      shd.vals <- c(shd.vals, tabu.results$shd)
      run.times <- c(run.times, tabu.results$run.time)
    }
    if (mean(shd.vals) < best.shd.mean) {
      best.shd.mean <- mean(shd.vals)
      best.shd.std <- sd(shd.vals)
      best.hyperparam <- tabu.frac
      av.run.time <- mean(run.times)
    }
  }
  return('mean'=best.shd.mean, 
         'std'=best.shd.std,
         'hyperparam'=best.hyperparam, 
         'run.time'=av.run.time)
}