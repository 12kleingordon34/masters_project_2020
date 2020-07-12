library(bnlearn)
library(equSA)
require(graph)
require(igraph)
library(parallel)
library(mixggm)

source("./utils.R")


run.simulation.test <- function(N, S, D.range, sparsity, function.list, 
                                alpha.range, tabu.range, gamma){
  #' Function which 

  num.algos <- length(function.list)
  
  # CPDAG Structural Hamming Distance
  mean.cpdag.shd.array <- array(numeric(), c(length(D.range), num.algos + 2))
  sd.cpdag.shd.array <- array(numeric(), c(length(D.range), num.algos + 2))
  # Moralized Structural Hamming Distance
  mean.moral.shd.array <- array(numeric(), c(length(D.range), num.algos + 2))
  sd.moral.shd.array <- array(numeric(), c(length(D.range), num.algos + 2))
  
  # BIC Criterion
  mean.bic.array <- array(numeric(), c(length(D.range), num.algos + 2))
  sd.bic.array <- array(numeric(), c(length(D.range), num.algos + 2))
  # False Discovery Rate
  mean.fdr.array <- array(numeric(), c(length(D.range), num.algos + 2))
  sd.fdr.array <- array(numeric(), c(length(D.range), num.algos + 2))
  # Sensitivity
  mean.sensitivity.array <- array(numeric(), c(length(D.range), num.algos + 2))
  sd.sensitivity.array <- array(numeric(), c(length(D.range), num.algos + 2))
  # Hyperparameters
  best.param.array <- array(numeric(), c(length(D.range), num.algos+1))  
  # Run time
  mean.run.time <- array(numeric(), c(length(D.range), num.algos+1))  
  
  colnames(mean.cpdag.shd.array) <- c(c("N", "D"), names(function.list))
  colnames(sd.cpdag.shd.array) <- c(c("N", "D"), names(function.list))
  colnames(mean.moral.shd.array) <- c(c("N", "D"), names(function.list))
  colnames(sd.moral.shd.array) <- c(c("N", "D"), names(function.list))
  colnames(mean.bic.array) <- c(c("N", "D"), names(function.list))
  colnames(sd.bic.array) <- c(c("N", "D"), names(function.list))  
  colnames(mean.fdr.array) <- c(c("N", "D"), names(function.list))
  colnames(sd.fdr.array) <- c(c("N", "D"), names(function.list))
  colnames(mean.sensitivity.array) <- c(c("N", "D"), names(function.list))
  colnames(sd.sensitivity.array) <- c(c("N", "D"), names(function.list))
  colnames(best.param.array) <- c("D", names(function.list))
  colnames(mean.run.time) <- c("D", names(function.list))
    
  d.num <- 1
  cl = makeCluster(2)
  for (d in D.range){
    mean.cpdag.shd.array[d.num, 1] <- N; sd.cpdag.shd.array[d.num, 1] <- N
    mean.cpdag.shd.array[d.num, 2] <- d; sd.cpdag.shd.array[d.num, 2] <- d
    mean.moral.shd.array[d.num, 1] <- N; sd.cpdag.shd.array[d.num, 1] <- N
    mean.moral.shd.array[d.num, 2] <- d; sd.cpdag.shd.array[d.num, 2] <- d
    best.param.array[d.num, 1] <- d; mean.run.time[d.num, 1] <- d

    for (algo in 1:num.algos){
      algo.name = names(function.list)[algo]
      if (substr(algo.name, 1, 4) == "tabu"){
        best.hyperparam <- 0; best.bic.mean <- -1e10
        for (tabu.frac in tabu.range){
          # Multiply tabu frac by number of nodes
          tabu <- ceiling(d * tabu.frac)
          cpdag.shd.vals <- c()
          moral.shd.vals <- c()
          bic.vals <- c()
          fdr.vals <- c()
          sensitivity.vals <- c()
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
            cpdag.shd <- shd(cpdag(tabu.dag), cpdag(sim.dag$graph))
            cpdag.shd.vals <- c(c(cpdag.shd.vals), cpdag.shd)
            moral.shd <- shd(moral(tabu.dag), moral(sim.dag$graph))
            moral.shd.vals <- c(c(moral.shd.vals), moral.shd)
            
            bic <- ebic(tabu.dag, sim.data.df, gamma=gamma)
            bic.vals <- c(c(bic.vals), bic)
            performance.metrics <- calculate.performance.statistics(
              cpdag(tabu.dag), cpdag(sim.dag$graph)
            )
            fdr.vals <- c(c(fdr.vals), performance.metrics$fdr)
            sensitivity.vals <- c(
              c(sensitivity.vals), performance.metrics$sensitivity
            )
          }
          if (mean(bic.vals) > best.bic.mean) {
            best.bic.mean <- mean(bic.vals)
            best.bic.std <- sd(bic.vals)
            best.cpdag.shd.mean <- mean(cpdag.shd.vals)
            best.cpdag.shd.std <- sd(cpdag.shd.vals)
            best.moral.shd.mean <- mean(moral.shd.vals)
            best.moral.shd.std <- sd(moral.shd.vals)
            best.fdr.mean <- mean(fdr.vals)
            best.std.mean <- sd(fdr.vals)
            best.sensitivity.mean <- mean(sensitivity.vals)
            best.sensitivity.std <- sd(sensitivity.vals)            
            best.hyperparam <- tabu.frac
            av.run.time <- mean(run.times)
          }
        }
      } else {
        best.hyperparam <- 0; best.bic.mean <- -1e10
        for (a in alpha.range) {
          cpdag.shd.vals <- c()
          moral.shd.vals <- c()
          bic.vals <- c()
          fdr.vals <- c()
          sensitivity.vals <- c()
          run.times <- c()
          for (s in 1:S){
            sim.dag <- generate.dag(N, d, sparsity, s)
            sim.data.df <- as.data.frame(sim.dag$data)    
            start.time <- Sys.time()
            print(sprintf(
              'D: %s. Sparsity: %s. Seed: %s. Algo: %s. Alpha Val: %s. Time: %s',
              d, sparsity,s,names(function.list)[algo], a, start.time
            ))       
            alpha.dag <- function.list[[algo]](sim.data.df, alpha=a, cluster=cl)  
            
            run.times <- c(run.times, Sys.time() - start.time)
            cpdag.shd <- shd(cpdag(alpha.dag), cpdag(sim.dag$graph))
            cpdag.shd.vals <- c(c(cpdag.shd.vals), cpdag.shd)
            moral.shd <- shd(moral(alpha.dag), moral(sim.dag$graph))
            moral.shd.vals <- c(c(moral.shd.vals), moral.shd)
            bic <- ebic(alpha.dag, sim.data.df, gamma=gamma)
            bic.vals <- c(c(bic.vals), bic)
            performance.metrics <- calculate.performance.statistics(
              cpdag(alpha.dag), cpdag(sim.dag$graph)
            )
            fdr.vals <- c(c(fdr.vals), performance.metrics$fdr)
            sensitivity.vals <- c(
              c(sensitivity.vals), performance.metrics$sensitivity
            )
          }
          if (mean(bic.vals) > best.bic.mean) {
            best.bic.mean <- mean(bic.vals)
            best.bic.std <- sd(bic.vals)
            best.cpdag.shd.mean <- mean(cpdag.shd.vals)
            best.cpdag.shd.std <- sd(cpdag.shd.vals)
            best.moral.shd.mean <- mean(moral.shd.vals)
            best.moral.shd.std <- sd(moral.shd.vals)
            best.fdr.mean <- mean(fdr.vals)
            best.fdr.std <- sd(fdr.vals)
            best.sensitivity.mean <- mean(sensitivity.vals)
            best.sensitivity.std <- sd(sensitivity.vals)            
            best.hyperparam <- a
            av.run.time <- mean(run.times)
          }
        }
      }
        mean.cpdag.shd.array[d.num, (algo+2)] <- best.cpdag.shd.mean
        sd.cpdag.shd.array[d.num, (algo+2)] <- best.cpdag.shd.std
        mean.moral.shd.array[d.num, (algo+2)] <- best.moral.shd.mean
        sd.moral.shd.array[d.num, (algo+2)] <- best.moral.shd.std
        mean.bic.array[d.num, (algo+2)] <- best.bic.mean
        sd.bic.array[d.num, (algo+2)] <- best.bic.std
        mean.fdr.array[d.num, (algo+2)] <- best.fdr.mean
        sd.fdr.array[d.num, (algo+2)] <- best.fdr.std
        mean.sensitivity.array[d.num, (algo+2)] <- best.sensitivity.mean
        sd.sensitivity.array[d.num, (algo+2)] <- best.sensitivity.std
        best.param.array[d.num, (algo+1)] <- best.hyperparam
        mean.run.time[d.num, (algo+1)] <- av.run.time
      }
    d.num <- d.num + 1
  }
  stopCluster(cl)
  return(list(
    'mean.bic'=mean.bic.array, 
    'std.bic'=sd.bic.array,    
    'mean.cpdag.shd'=mean.cpdag.shd.array, 
    'std.cpdag.shd'=sd.cpdag.shd.array,    
    'mean.moral.shd'=mean.moral.shd.array, 
    'std.moral.shd'=sd.moral.shd.array,
    'mean.fdr'=mean.fdr.array,
    'std.fdr'=sd.fdr.array,
    'mean.sensitivity'=mean.sensitivity.array,
    'std.sensitivity'=sd.sensitivity.array,
    'best.params'=best.param.array, 
    "times"=mean.run.time
  ))
}


tabu.xval.test <- function(N, S, D.range, true.sparsity, hyperparam.sparsity, constraint.algo){
  mean.cpdag.shd.array <- array(numeric(), c(length(D.range),5))
  sd.cpdag.shd.array <- array(numeric(), c(length(D.range), 5))
  best.param.array <- array(numeric(), c(length(D.range), 5))  
  mean.run.time <- array(numeric(), c(length(D.range), 4))  
  
  
  colnames(mean.cpdag.shd.array) <- c(c("N", "D"), c("tabu.null", "tabu.init", "tabu.constraint"))
  colnames(mean.cpdag.shd.array) <- c(c("N", "D"), c("tabu.null", "tabu.init", "tabu.constraint"))
  colnames(best.param.array) <- c("D", c("tabu.null", "tabu.init.tabu", "tabu.init.sparsity", "tabu.constraint"))
  colnames(mean.run.time) <- c("D", c("tabu.null", "tabu.init", "tabu.constraint"))
  
  d.num <- 1
  for (d in D.range){
    mean.cpdag.shd.array[d.num, 1] <- N; sd.cpdag.shd.array[d.num, 1] <- N
    mean.cpdag.shd.array[d.num, 2] <- d; sd.cpdag.shd.array[d.num, 2] <- d
    best.param.array[d.num, 1] <- d; mean.run.time[d.num, 1] <- d
    
    # Run for tabu.null
    tabu.null <- tabu.hyperparam.iteration(N, d, S, sparsity, tabu.range, start=NULL, constraint.algo=NULL)
    mean.cpdag.shd.array[d.num, 3] <- tabu.null$mean
    sd.cpdag.shd.array[d.num, 3] <- tabu.null$std
    best.param.array[d.num, 2] <- tabu.null$hyperparam
    mean.run.time[d.num, 2] <- tabu.null$av.run.time  
    
    # Run for tabu.random
    tabu.random <- tabu.hyperparam.iteration(N, d, S, sparsity, tabu.range, start="random", constraint.algo=NULL)
    mean.cpdag.shd.array[d.num, 4] <- tabu.random$mean
    sd.cpdag.shd.array[d.num, 4] <- tabu.random$std
    best.param.array[d.num, 3] <- tabu.random$hyperparam
    best.param.array[d.num, 4] <- 0.2
    mean.run.time[d.num, 3] <- tabu.random$av.run.time  
    
    # # Run for tabu.constraint
    # tabu.constraint <- tabu.hyperparam.iteration(N, d, S, sparsity, tabu.range, start=NULL, constraint.algo=constraint.algo)
    # mean.cpdag.shd.array[d.num, 5] <- tabu.constraint$mean
    # sd.cpdag.shd.array[d.num, 5] <- tabu.constraint$std
    # best.param.array[d.num, 5] <- tabu.constraint$hyperparam
    # mean.run.time[d.num, 4] <- tabu.constraint$av.run.time  
    
    d.num <- d.num + 1
  } 
  return(list('mean'=mean.cpdag.shd.array, 'std'=sd.cpdag.shd.array, 'best.params'=best.param.array, "times"=mean.run.time))
}

tabu.learning.iteration <- function(N, d, s, sparsity, tabu, start=NULL, constraint.algo=NULL){
  sim.dag <- generate.dag(N, d, sparsity, s)
  sim.data.df <- as.data.frame(sim.dag$data)
  start.time <- Sys.time()
  print(sprintf(
    'D: %s. Sparsity: %s. Seed: %s. Tabu Val: %s. Time: %s',
    d, sparsity, s, tabu, start.time
  ))
  if (is.null(start)){
    if (!is.null(constraint.algo)){
      start <- constraint.algo(sim.data.df)
    }
  } else if (start=='random'){
    set.seed(s)
    start <- random.graph(nodes(sim.dag$graph), 1, method='ordered')
  }
  tabu.dag <- tabu(
    sim.data.df, 
    tabu=tabu,
    score='bge',
    max.tabu=tabu,
    start=start
  )
  run.time <- Sys.time() - start.time
  cpdag.shd <- cpdag.shd(cpdag(tabu.dag), cpdag(sim.dag$graph))
  return(list("cpdag.shd"=cpdag.shd, "run.time"=run.time))
}


tabu.hyperparam.iteration <- function(N, d, S, sparsity, tabu.range, start=NULL, constraint.algo=NULL){
  #'
  best.cpdag.shd.mean <- 1e9
  for (tabu.frac in tabu.range){
    # Multiply tabu frac by number of nodes
    tabu <- ceiling(d * tabu.frac)
    cpdag.shd.vals <- c()
    run.times <- c()
    for (s in 1:S){
      tabu.results <- tabu.learning.iteration(N, d, s, sparsity, tabu, start, constraint.algo)
      cpdag.shd.vals <- c(cpdag.shd.vals, tabu.results$cpdag.shd)
      run.times <- c(run.times, tabu.results$run.time)
    }
    if (mean(cpdag.shd.vals) < best.cpdag.shd.mean) {
      best.cpdag.shd.mean <- mean(cpdag.shd.vals)
      best.cpdag.shd.std <- sd(cpdag.shd.vals)
      best.hyperparam <- tabu.frac
      av.run.time <- mean(run.times)
    }
  }
  return(list(
     'mean'=best.cpdag.shd.mean, 
     'std'=best.cpdag.shd.std,
     'hyperparam'=best.hyperparam, 
     'av.run.time'=av.run.time))
}

hc.xval.test <- function(N, S, D.range, true.sparsity, hyperparam.sparsity, constraint.algo){
  mean.cpdag.shd.array <- array(numeric(), c(length(D.range),5))
  sd.cpdag.shd.array <- array(numeric(), c(length(D.range), 5))
  best.param.array <- array(numeric(), c(length(D.range), 4))  
  mean.run.time <- array(numeric(), c(length(D.range), 4))  
  
  
  colnames(mean.cpdag.shd.array) <- c(c("N", "D"), c("hc.null", "hc.init", "hc.constraint"))
  colnames(mean.cpdag.shd.array) <- c(c("N", "D"), c("hc.null", "hc.init", "hc.constraint"))
  colnames(best.param.array) <- c("D", c("hc.null", "hc.init.sparsity", "hc.constraint"))
  colnames(mean.run.time) <- c("D", c("hc.null", "hc.init", "hc.constraint"))
  
  d.num <- 1
  for (d in D.range){
    mean.cpdag.shd.array[d.num, 1] <- N; sd.cpdag.shd.array[d.num, 1] <- N
    mean.cpdag.shd.array[d.num, 2] <- d; sd.cpdag.shd.array[d.num, 2] <- d
    best.param.array[d.num, 1] <- d; mean.run.time[d.num, 1] <- d
    
    # Run for hc.null
    hc.null <- hc.hyperparam.iteration(N, d, S, sparsity, start=NULL, constraint.algo=NULL)
    mean.cpdag.shd.array[d.num, 3] <- hc.null$mean
    sd.cpdag.shd.array[d.num, 3] <- hc.null$std
    best.param.array[d.num, 2] <- hc.null$hyperparam
    mean.run.time[d.num, 2] <- hc.null$av.run.time  
    
    # Run for hc.random
    hc.random <- hc.hyperparam.iteration(N, d, S, sparsity, start="random", constraint.algo=NULL)
    mean.cpdag.shd.array[d.num, 4] <- hc.random$mean
    sd.cpdag.shd.array[d.num, 4] <- hc.random$std
    best.param.array[d.num, 3] <- hc.random$hyperparam
    mean.run.time[d.num, 3] <- hc.random$av.run.time  
    
    # # Run for hc.constraint
    # hc.constraint <- hc.hyperparam.iteration(N, d, S, sparsity, start=NULL, constraint.algo=constraint.algo)
    # mean.cpdag.shd.array[d.num, 5] <- hc.constraint$mean
    # sd.cpdag.shd.array[d.num, 5] <- hc.constraint$std
    # best.param.array[d.num, 4] <- hc.constraint$hyperparam
    # mean.run.time[d.num, 4] <- hc.constraint$av.run.time  
    
    d.num <- d.num + 1
  } 
  return(list('mean'=mean.cpdag.shd.array, 'std'=sd.cpdag.shd.array, 'best.params'=best.param.array, "times"=mean.run.time))
}

hc.learning.iteration <- function(N, d, s, sparsity, start=NULL, constraint.algo=NULL){
  sim.dag <- generate.dag(N, d, sparsity, s)
  sim.data.df <- as.data.frame(sim.dag$data)
  start.time <- Sys.time()
  print(sprintf(
    'D: %s. Sparsity: %s. Seed: %s. Time: %s',
    d, sparsity, s, start.time
  ))
  if (is.null(start)){
    if (!is.null(constraint.algo)){
      start <- constraint.algo(sim.data.df)
    }
  } else if (start=='random'){
    set.seed(s)
    start <- random.graph(nodes(sim.dag$graph), 1, method='ordered')
  }
  hc.dag <- hc(
    sim.data.df, 
    score='bge',
    restart=10,
    peturb=5,
    start=start
  )
  run.time <- Sys.time() - start.time
  cpdag.shd <- cpdag.shd(cpdag(hc.dag), cpdag(sim.dag$graph))
  return(list("cpdag.shd"=cpdag.shd, "run.time"=run.time))
}


hc.hyperparam.iteration <- function(N, d, S, sparsity, start=NULL, constraint.algo=NULL){
  #'
  best.cpdag.shd.mean <- 1e9
  # Multiply hc frac by number of nodes
  hc <- ceiling(d * hc.frac)
  cpdag.shd.vals <- c()
  run.times <- c()
  for (s in 1:S){
    hc.results <- hc.learning.iteration(N, d, s, sparsity, start, constraint.algo)
    cpdag.shd.vals <- c(cpdag.shd.vals, hc.results$cpdag.shd)
    run.times <- c(run.times, hc.results$run.time)
  }
  if (mean(cpdag.shd.vals) < best.cpdag.shd.mean) {
    best.cpdag.shd.mean <- mean(cpdag.shd.vals)
    best.cpdag.shd.std <- sd(cpdag.shd.vals)
    av.run.time <- mean(run.times)
  }
  return(list(
    'mean'=best.cpdag.shd.mean, 
    'std'=best.cpdag.shd.std,
    'av.run.time'=av.run.time))
}

alpha.constraint.xval <- function(functions.bn, alpha.range, n.runs, d.range, max.n, gamma=0){
  #' For `n.runs` iterations, sample D (the number of graph nodes) uniformally
  #' between the two values in the list d.range, and sample 
  #' N (number of observations) uniformally between values in n.range. 
  #' For each algorithm in functions.bn, learn a graph for each alpha within
  #' alpha.range, and calculate the 
  #' * CPDAG SHD
  #' * moralised SHD 
  #' * BIC criterion for the moralised graph. 
  #' * CPDAG False Discovery Rate (FDR)
  #' * CPDAG Sensitivity
  num.algos <- length(functions.bn)
  
  # CPDAG Structural Hamming Distance
  cpdag.shd.array <- array(numeric(), c(0, num.algos + 4))
  # Moralized Structural Hamming Distance
  moral.shd.array <- array(numeric(), c(0, num.algos + 4))
  # BIC Criterion
  bic.array <- array(numeric(), c(0, num.algos + 4))
  # False Discovery Rate
  fdr.array <- array(numeric(), c(0, num.algos + 4))
  # Sensitivity
  sensitivity.array <- array(numeric(), c(0, num.algos + 4))

  colnames(cpdag.shd.array) <- c(c("N", "D", "sparsity", "alpha"), names(functions.bn))
  colnames(moral.shd.array) <- c(c("N", "D", "sparsity", "alpha"), names(functions.bn))
  colnames(bic.array) <- c(c("N", "D", "sparsity", "alpha"), names(functions.bn))
  colnames(fdr.array) <- c(c("N", "D", "sparsity", "alpha"), names(functions.bn))
  colnames(sensitivity.array) <- c(c("N", "D", "sparsity", "alpha"), names(functions.bn))

  for (n in seq(1, n.runs)){
    D <- sample(d.range[1]:d.range[2], 1)
    N <- sample(D:max.n, 1)
    sparsity = 3/(D-1)
    dag <- generate.dag(N, D, sparsity, seed=n)
    true.graph <- dag$graph; sim.data <- dag$data
    row <- c(N, D, sparsity)
    for (a in alpha.range){
      sim.data.df <- as.data.frame(sim.data)    
      start.time <- Sys.time()
      row <- c(row, alpha)
      row.cpdag.shd <- row
      row.moral.shd <- row
      bic.row <- row
      fdr.row <- row
      sensitivity.row <- row
      for (algo in 1:length(functions.bn)){
        print(sprintf(
          'Run number: %s/%s. Algo: %s Alpha Val: %s. Time: %s',
          n, n.runs, names(functions.bn)[algo], a, start.time
        ))       
        learned.dag <- functions.bn[[algo]](sim.data.df, alpha=a)  
        cpdag.shd <- shd(cpdag(learned.dag), cpdag(true.graph))
        row.cpdag.shd <- c(row.cpdag.shd, cpdag.shd)
        moral.shd <- shd(moral(learned.dag), moral(true.graph))
        row.moral.shd <- c(row.moral.shd, moral.shd)
        bic <- ebic(learned.dag, sim.data.df, gamma=gamma)
        bic.row <- c(bic.row, bic)
        performance.metrics <- calculate.performance.statistics(
          cpdag(learned.dag), cpdag(true.graph)
        )
        fdr.row <- c(fdr.row, performance.metrics$fdr)
        sensitivity.row <- c(sensitivity.row, performance.metrics$sensitivity)        
      } # For algo loop
      cpdag.shd.array <- rbind(cpdag.shd.array, row.cpdag.shd)
      moral.shd.array <- rbind(moral.shd.array, row.moral.shd)
      bic.array <- rbind(bic.array, bic.row)
      fdr.array <- rbind(fdr.array, fdr.row)
      sensitivity.array <- rbind(sensitivity.array, sensitivity.row)
    } # FOR alpha.range
  } # FOR n.run
  return(list(
    'cpdag'=cpdag.shd.array,
    'moral'=moral.shd.array,
    'bic'=bic.array,
    'fdr'=fdr.array,
    'sensitivity'=sensitivity.array
  ))
}
