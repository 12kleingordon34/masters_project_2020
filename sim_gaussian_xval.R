library(bnlearn)
library(equSA)
require(graph)
require(igraph)
library(parallel)
library(mixggm)
library(reticulate)
import('scipy')
source_python('./notears/notears/linear.py')
source_python('./notears/notears/utils.py')

source("./utils.R")


data.size.study.gaussian <- function(algorithms, data.size.ranges, d.range, n.runs){
  #' Take a fitted BN and generate `n.data.size` samples, and learn
  #' a network across all algorithms in function.list. Repeat
  #' over n.runs. 
  #' 
  #' Returns arrays of the SHD for CPDAG and Moralised graphs.
  num.algos <- length(algorithms)
  # CPDAG Structural Hamming Distance
  shd.array <- array(numeric(), c(0, num.algos + 3))
  # Moralized Structural Hamming Distance
  hd.array <- array(numeric(), c(0, num.algos + 3))
  # BIC Criterion
  bic.array <- array(numeric(), c(0, num.algos + 3))
  # False Discovery Rate
  fdr.array <- array(numeric(), c(0, num.algos + 3))
  # Sensitivity
  sensitivity.array <- array(numeric(), c(0, num.algos + 3))
  
  colnames(shd.array) <- c(c("run", "N", "D"), names(algorithms))
  colnames(hd.array) <- c(c("run", "N", "D"), names(algorithms))
  colnames(fdr.array) <- c(c("run", "N", "D"), names(algorithms))
  colnames(sensitivity.array) <- c(c("run", "N", "D"), names(algorithms))
  for (D in d.range){
    for (data.size in data.size.ranges){
      for (r in seq(1, n.runs, 1)){
        algo.graphs <- vector('list', length(algorithms))
        set.seed(r)
        sparsity <- 3/(D-1)
        sim.dag <- generate.dag(data.size, D, sparsity, seed=D+r)
        true.graph <- sim.dag$graph; sim.data <- sim.dag$data; sim.data.df <- as.data.frame(sim.data)
        
        row <- c(r, data.size, D)
        row.shd <- row;  row.hd <- row;  bic.row <- row; fdr.row <- row
        sensitivity.row <- row
        for (algo in seq(1, length(algorithms), 1)){
          start.time <- Sys.time()
          print(sprintf(
            'Run number: %s/%s. Algo: %s. Data size: %s/%s. Network size: %s/%s Time: %s',
            r, n.runs, names(algorithms)[algo],
            data.size, data.size.ranges[length(data.size.ranges)], 
            D, d.range[length(d.range)],
            start.time
          ))
          if (grepl('notears', names(algorithms)[algo], fixed=TRUE)){
            amat.dag <- algorithms[[algo]](sim.data)
            colnames(amat.dag) <- colnames(sim.data)
            fitted.network <- convert.amat.to.bn(amat.dag)
          } else if (grepl('tabu', names(algorithms)[algo], fixed=TRUE)){
            fitted.network <- algorithms[[algo]](
              sim.data.df,
              tabu=D,
              score='bge'
            )
          } else fitted.network <- algorithms[[algo]](sim.data.df)
          shd <- bnlearn::shd(cpdag(fitted.network), cpdag(true.graph))
          row.shd <- c(row.shd, shd)
          hd <- bnlearn::hamming(moral(fitted.network), moral(true.graph))
          row.hd <- c(row.hd, hd)
          performance.metrics <- calculate.performance.statistics(
            cpdag(fitted.network), cpdag(true.graph)
          )
          fdr.row <- c(fdr.row, performance.metrics$fdr)
          sensitivity.row <- c(sensitivity.row, performance.metrics$sensitivity)        
        } # For algo loop
        shd.array <- rbind(shd.array, row.shd)
        hd.array <- rbind(hd.array, row.hd)
        fdr.array <- rbind(fdr.array, fdr.row)
        sensitivity.array <- rbind(sensitivity.array, sensitivity.row)
      } # FOR runs
    } # For Data size
  } # FOR d.range
  return(list(
    'shd'=data.frame(shd.array),
    'hd'=data.frame(hd.array),
    'fdr'=fdr.array,
    'sensitivity'=sensitivity.array
  ))
}


run.simulation.test <- function(N, S, D.range, sparsity, function.list, 
                                alpha.range, tabu.range, restart.range,
                                l1.range, gamma){
  #' Function which 

  num.algos <- length(function.list)
  
  # CPDAG Structural Hamming Distance
  mean.shd.array <- array(numeric(), c(length(D.range), num.algos + 2))
  sd.shd.array <- array(numeric(), c(length(D.range), num.algos + 2))
  # Moralized Structural Hamming Distance
  mean.hd.array <- array(numeric(), c(length(D.range), num.algos + 2))
  sd.hd.array <- array(numeric(), c(length(D.range), num.algos + 2))
  
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
  
  colnames(mean.shd.array) <- c(c("N", "D"), names(function.list))
  colnames(sd.shd.array) <- c(c("N", "D"), names(function.list))
  colnames(mean.hd.array) <- c(c("N", "D"), names(function.list))
  colnames(sd.hd.array) <- c(c("N", "D"), names(function.list))
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
    mean.shd.array[d.num, 1] <- N; sd.shd.array[d.num, 1] <- N
    mean.shd.array[d.num, 2] <- d; sd.shd.array[d.num, 2] <- d
    mean.hd.array[d.num, 1] <- N; sd.shd.array[d.num, 1] <- N
    mean.hd.array[d.num, 2] <- d; sd.shd.array[d.num, 2] <- d
    best.param.array[d.num, 1] <- d; mean.run.time[d.num, 1] <- d
    mean.bic.array[d.num, 1] <- N; sd.bic.array[d.num, 1] <- N
    mean.bic.array[d.num, 2] <- d; sd.bic.array[d.num, 2] <- d
    

    for (algo in 1:num.algos){
      algo.name = names(function.list)[algo]
      if (substr(algo.name, 1, 4) == "hc"){
        best.hyperparam <- 0; best.bic.mean <- -1e10
        for (restart.no in restart.range){
          # Multiply tabu frac by number of nodes
          tabu <- ceiling(d * tabu.frac)
          shd.vals <- c()
          hd.vals <- c()
          bic.vals <- c()
          fdr.vals <- c()
          sensitivity.vals <- c()
          run.times <- c()
          for (s in 1:S){
            sim.dag <- generate.dag(N, d, sparsity, s)
            sim.data.df <- as.data.frame(sim.dag$data)
            start.time <- Sys.time()
            print(sprintf(
              'D: %s. Sparsity: %s. Seed: %s. Algo: %s. Restart no: %s. Time: %s',
              d, sparsity,s,names(function.list)[algo], restart.no, start.time
            ))
            tabu.dag <- function.list[[algo]](
              sim.data.df, 
              restart=restart.no
            )
            run.times <- c(run.times, Sys.time() - start.time)
            shd <- bnlearn::shd(cpdag(tabu.dag), cpdag(sim.dag$graph))
            shd.vals <- c(c(shd.vals), shd)
            hd <- bnlearn::shd(moral(tabu.dag), moral(sim.dag$graph))
            hd.vals <- c(c(hd.vals), hd)
            
            bic <- ebic(tabu.dag, sim.data.df, gamma=gamma)
            bic.vals <- c(c(bic.vals), bic)
            performance.metrics <- calculate.performance.statistics(
              cpdag(tabu.dag), cpdag(sim.dag$graph)
            )
            fdr.vals <- c(c(fdr.vals), performance.metrics$fdr)
            sensitivity.vals <- c(
              c(sensitivity.vals), performance.metrics$sensitivity
            )
          } # for s in 1:s
          if (mean(bic.vals) > best.bic.mean) {
            best.bic.mean <- mean(bic.vals)
            best.bic.std <- sd(bic.vals)
            best.shd.mean <- mean(shd.vals)
            best.shd.std <- sd(shd.vals)
            best.hd.mean <- mean(hd.vals)
            best.hd.std <- sd(hd.vals)
            best.fdr.mean <- mean(fdr.vals)
            best.std.mean <- sd(fdr.vals)
            best.sensitivity.mean <- mean(sensitivity.vals)
            best.sensitivity.std <- sd(sensitivity.vals)            
            best.hyperparam <- restart.no
            av.run.time <- mean(run.times)
          } # if mean(bic) > best.bic
        } # for tabu ftac
      
      } else if (substr(algo.name, 1, 4) == "tabu"){
        best.hyperparam <- 0; best.bic.mean <- -1e10
        for (tabu.frac in tabu.range){
          # Multiply tabu frac by number of nodes
          tabu <- ceiling(d * tabu.frac)
          shd.vals <- c()
          hd.vals <- c()
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
            shd <- bnlearn::shd(cpdag(tabu.dag), cpdag(sim.dag$graph))
            shd.vals <- c(c(shd.vals), shd)
            hd <- bnlearn::shd(moral(tabu.dag), moral(sim.dag$graph))
            hd.vals <- c(c(hd.vals), hd)
            
            bic <- ebic(tabu.dag, sim.data.df, gamma=gamma)
            bic.vals <- c(c(bic.vals), bic)
            performance.metrics <- calculate.performance.statistics(
              cpdag(tabu.dag), cpdag(sim.dag$graph)
            )
            fdr.vals <- c(c(fdr.vals), performance.metrics$fdr)
            sensitivity.vals <- c(
              c(sensitivity.vals), performance.metrics$sensitivity
            )
          } # for s in 1:s
          if (mean(bic.vals) > best.bic.mean) {
            best.bic.mean <- mean(bic.vals)
            best.bic.std <- sd(bic.vals)
            best.shd.mean <- mean(shd.vals)
            best.shd.std <- sd(shd.vals)
            best.hd.mean <- mean(hd.vals)
            best.hd.std <- sd(hd.vals)
            best.fdr.mean <- mean(fdr.vals)
            best.std.mean <- sd(fdr.vals)
            best.sensitivity.mean <- mean(sensitivity.vals)
            best.sensitivity.std <- sd(sensitivity.vals)            
            best.hyperparam <- tabu.frac
            av.run.time <- mean(run.times)
          } # if mean(bic) > best.bic
        } # for tabu ftac
      } else if (substr(algo.name, 1, 7) == 'notears'){
        best.hyperparam <- 0; best.bic.mean <- -1e10
        for (l1 in l1.range) {
          shd.vals <- c()
          hd.vals <- c()
          bic.vals <- c()
          fdr.vals <- c()
          sensitivity.vals <- c()
          run.times <- c()
          for (s in 1:S){
            sim.dag <- generate.dag(N, d, sparsity, s)
            sim.data <- sim.dag$data
            sim.data.df <- as.data.frame(sim.dag$data)
            start.time <- Sys.time()
            print(sprintf(
              'D: %s. Sparsity: %s. Seed: %s. Algo: %s. L1 Val: %s. Time: %s',
              d, sparsity,s,names(function.list)[algo], l1, start.time
            ))       
            amat.dag <- notears_linear(sim.data, lambda1=l1, loss_type='l2')  
            colnames(amat.dag) <- colnames(sim.data)
            alpha.dag <- convert.amat.to.bn(amat.dag)
            run.times <- c(run.times, Sys.time() - start.time)
            shd <- bnlearn::shd(cpdag(alpha.dag), cpdag(sim.dag$graph))
            shd.vals <- c(c(shd.vals), shd)
            hd <- bnlearn::shd(moral(alpha.dag), moral(sim.dag$graph))
            hd.vals <- c(c(hd.vals), hd)
            bic <- ebic(alpha.dag, sim.data.df, gamma=gamma)
            bic.vals <- c(c(bic.vals), bic)
            performance.metrics <- calculate.performance.statistics(
              cpdag(alpha.dag), cpdag(sim.dag$graph)
            )
            fdr.vals <- c(c(fdr.vals), performance.metrics$fdr)
            sensitivity.vals <- c(
              c(sensitivity.vals), performance.metrics$sensitivity
            )
          } # for s in 1:S
          if (mean(bic.vals) > best.bic.mean) {
            best.bic.mean <- mean(bic.vals)
            best.bic.std <- sd(bic.vals)
            best.shd.mean <- mean(shd.vals)
            best.shd.std <- sd(shd.vals)
            best.hd.mean <- mean(hd.vals)
            best.hd.std <- sd(hd.vals)
            best.fdr.mean <- mean(fdr.vals)
            best.fdr.std <- sd(fdr.vals)
            best.sensitivity.mean <- mean(sensitivity.vals)
            best.sensitivity.std <- sd(sensitivity.vals)            
            best.hyperparam <- l1
            av.run.time <- mean(run.times)
          } # if mean(bic) > best.bic
        }
      } else if (substr(algo.name, 1, 7) == 'lingam'){
        best.hyperparam <- 0; best.bic.mean <- -1e10
        shd.vals <- c()
        hd.vals <- c()
        bic.vals <- c()
        fdr.vals <- c()
        sensitivity.vals <- c()
        run.times <- c()
        for (s in 1:S){
          sim.dag <- generate.dag(N, d, sparsity, s)
          sim.data <- sim.dag$data
          sim.data.df <- as.data.frame(sim.dag$data)    
          start.time <- Sys.time()
          print(sprintf(
            'D: %s. Sparsity: %s. Seed: %s. Algo: %s. Time: %s',
            d, sparsity, s, names(function.list)[algo], start.time
          ))       
          alpha.dag <- function.list[[algo]](sim.data.df)  
          
          run.times <- c(run.times, Sys.time() - start.time)
          shd <- bnlearn::shd(cpdag(alpha.dag), cpdag(sim.dag$graph))
          shd.vals <- c(c(shd.vals), shd)
          hd <- bnlearn::shd(moral(alpha.dag), moral(sim.dag$graph))
          hd.vals <- c(c(hd.vals), hd)
          bic <- ebic(alpha.dag, sim.data.df, gamma=gamma)
          bic.vals <- c(c(bic.vals), bic)
          performance.metrics <- calculate.performance.statistics(
            cpdag(alpha.dag), cpdag(sim.dag$graph)
          )
          fdr.vals <- c(c(fdr.vals), performance.metrics$fdr)
          sensitivity.vals <- c(
            c(sensitivity.vals), performance.metrics$sensitivity
          )
        } # for s in 1:S
        if (mean(bic.vals) > best.bic.mean) {
          best.bic.mean <- mean(bic.vals)
          best.bic.std <- sd(bic.vals)
          best.shd.mean <- mean(shd.vals)
          best.shd.std <- sd(shd.vals)
          best.hd.mean <- mean(hd.vals)
          best.hd.std <- sd(hd.vals)
          best.fdr.mean <- mean(fdr.vals)
          best.fdr.std <- sd(fdr.vals)
          best.sensitivity.mean <- mean(sensitivity.vals)
          best.sensitivity.std <- sd(sensitivity.vals)            
          best.hyperparam <- -999999
          av.run.time <- mean(run.times)
        } # if mean(bic) > best.bic
      } else {
        best.hyperparam <- 0; best.bic.mean <- -1e10
        for (a in alpha.range) {
          shd.vals <- c()
          hd.vals <- c()
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
            shd <- bnlearn::shd(cpdag(alpha.dag), cpdag(sim.dag$graph))
            shd.vals <- c(c(shd.vals), shd)
            hd <- bnlearn::shd(moral(alpha.dag), moral(sim.dag$graph))
            hd.vals <- c(c(hd.vals), hd)
            bic <- ebic(alpha.dag, sim.data.df, gamma=gamma)
            bic.vals <- c(c(bic.vals), bic)
            performance.metrics <- calculate.performance.statistics(
              cpdag(alpha.dag), cpdag(sim.dag$graph)
            )
            fdr.vals <- c(c(fdr.vals), performance.metrics$fdr)
            sensitivity.vals <- c(
              c(sensitivity.vals), performance.metrics$sensitivity
            )
          } # for s in 1:S
          if (mean(bic.vals) > best.bic.mean) {
            best.bic.mean <- mean(bic.vals)
            best.bic.std <- sd(bic.vals)
            best.shd.mean <- mean(shd.vals)
            best.shd.std <- sd(shd.vals)
            best.hd.mean <- mean(hd.vals)
            best.hd.std <- sd(hd.vals)
            best.fdr.mean <- mean(fdr.vals)
            best.fdr.std <- sd(fdr.vals)
            best.sensitivity.mean <- mean(sensitivity.vals)
            best.sensitivity.std <- sd(sensitivity.vals)            
            best.hyperparam <- a
            av.run.time <- mean(run.times)
          } # if mean(bic) > best.bic
        } # for a in alpha.range
      
        }
        mean.shd.array[d.num, (algo+2)] <- best.shd.mean
        sd.shd.array[d.num, (algo+2)] <- best.shd.std
        mean.hd.array[d.num, (algo+2)] <- best.hd.mean
        sd.hd.array[d.num, (algo+2)] <- best.hd.std
        mean.bic.array[d.num, (algo+2)] <- best.bic.mean
        sd.bic.array[d.num, (algo+2)] <- best.bic.std
        mean.fdr.array[d.num, (algo+2)] <- best.fdr.mean
        sd.fdr.array[d.num, (algo+2)] <- best.fdr.std
        mean.sensitivity.array[d.num, (algo+2)] <- best.sensitivity.mean
        sd.sensitivity.array[d.num, (algo+2)] <- best.sensitivity.std
        best.param.array[d.num, (algo+1)] <- best.hyperparam
        mean.run.time[d.num, (algo+1)] <- av.run.time
      } #
    d.num <- d.num + 1
  } # for d in d.range
  stopCluster(cl)
  return(list(
    'mean.bic'=mean.bic.array, 
    'std.bic'=sd.bic.array,    
    'mean.shd'=mean.shd.array, 
    'std.shd'=sd.shd.array,    
    'mean.hd'=mean.hd.array, 
    'std.hd'=sd.hd.array,
    'mean.fdr'=mean.fdr.array,
    'std.fdr'=sd.fdr.array,
    'mean.sensitivity'=mean.sensitivity.array,
    'std.sensitivity'=sd.sensitivity.array,
    'best.params'=best.param.array, 
    "times"=mean.run.time
  ))
}


tabu.xval.test <- function(N, S, D.range, sparsity, tabu.range){
  #' Learn graph structureof randomly generated DAG (with simulated
  #' data size N and number of nodes in `D.range` using TABU search with
  #' a range of hyperparameters `tabu.range`.
  mean.shd.array <- array(numeric(), c(length(D.range),4))
  sd.shd.array <- array(numeric(), c(length(D.range), 4))
  best.param.array <- array(numeric(), c(length(D.range), 3))  
  mean.run.time <- array(numeric(), c(length(D.range), 3))  
  
  
  colnames(mean.shd.array) <- c(c("N", "D"), c("tabu.null", "tabu.init"))
  colnames(mean.shd.array) <- c(c("N", "D"), c("tabu.null", "tabu.init"))
  colnames(best.param.array) <- c("D", c("tabu.null", "tabu.init.tabu"))
  colnames(mean.run.time) <- c("D", c("tabu.null", "tabu.init"))
  
  d.num <- 1
  for (d in D.range){
    mean.shd.array[d.num, 1] <- N; sd.shd.array[d.num, 1] <- N
    mean.shd.array[d.num, 2] <- d; sd.shd.array[d.num, 2] <- d
    best.param.array[d.num, 1] <- d; mean.run.time[d.num, 1] <- d
    
    # Run for tabu.null
    tabu.null <- tabu.hyperparam.iteration(N, d, S, sparsity, tabu.range, start=NULL)
    mean.shd.array[d.num, 3] <- tabu.null$mean
    sd.shd.array[d.num, 3] <- tabu.null$std
    best.param.array[d.num, 2] <- tabu.null$hyperparam
    mean.run.time[d.num, 2] <- tabu.null$av.run.time  
    
    # Run for tabu.random
    tabu.random <- tabu.hyperparam.iteration(N, d, S, sparsity, tabu.range, start="random")
    mean.shd.array[d.num, 4] <- tabu.random$mean
    sd.shd.array[d.num, 4] <- tabu.random$std
    best.param.array[d.num, 3] <- tabu.random$hyperparam
    mean.run.time[d.num, 3] <- tabu.random$av.run.time  
    d.num <- d.num + 1
  } 
  return(list('mean'=mean.shd.array, 'std'=sd.shd.array, 'best.params'=best.param.array, "times"=mean.run.time))
}

tabu.learning.iteration <- function(N, d, s, sparsity, tabu.length, start=NULL){
  #' Generate random DAG parameterised by number of nodes `d`
  #' and `sparsity`. Learn graph structure using TABU
  #' search paramterised by `tabu.length` list length and `tabu.length`
  #' as the max.tabu size.
  sim.dag <- generate.dag(N, d, sparsity, s)
  sim.data.df <- as.data.frame(sim.dag$data)
  start.time <- Sys.time()
  print(sprintf(
    'D: %s. Sparsity: %s. Seed: %s. Tabu Length Val: %s. Time: %s',
    d, sparsity, s, tabu.length, start.time
  ))
  if (!is.null(start)){
    set.seed(s)
    initial.graph <- random.graph(nodes(sim.dag$graph), 1, method='ordered')
  } else initial.graph <- random.graph(nodes(sim.dag$graph), 1, method='empty')
  tabu.dag <- tabu(
    sim.data.df, 
    tabu=tabu.length,
    score='bic-g',
    max.tabu=tabu.length,
    start=initial.graph
  )

  run.time <- Sys.time() - start.time
  shd <- shd(cpdag(tabu.dag), cpdag(sim.dag$graph))
  return(list("shd"=shd, "run.time"=run.time))
}

tabu.hyperparam.iteration <- function(N, d, S, sparsity, tabu.range, start=NULL){
  #' Iterate over tabu list size fractions in `tabu.range`, find the SHD
  #' for tabu-learned graphs learnt from simulated data of size `N`x`d`.
  #' Repeat over `S` different randomly generated graphs. If `start`='random',
  #' the Tabu initialisation will be a random DAG.
  best.shd.mean <- 1e9
  for (tabu.frac in tabu.range){
    # Multiply tabu frac by number of nodes
    tabu.length <- ceiling(d * tabu.frac)
    shd.vals <- c()
    run.times <- c()
    for (s in 1:S){
      tabu.results <- tabu.learning.iteration(N, d, s, sparsity, tabu.length, start)
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
  return(list(
     'mean'=best.shd.mean, 
     'std'=best.shd.std,
     'hyperparam'=best.hyperparam, 
     'av.run.time'=av.run.time))
}

hc.xval.test <- function(N, S, D.range, sparsity, restart.hyperparams){
  mean.shd.array <- array(numeric(), c(length(D.range), 4))
  sd.shd.array <- array(numeric(), c(length(D.range), 4))
  best.param.array <- array(numeric(), c(length(D.range), 3))  
  mean.run.time <- array(numeric(), c(length(D.range), 3))
  
  
  colnames(mean.shd.array) <- c(c("N", "D"), c("hc.null", "hc.init"))
  colnames(mean.shd.array) <- c(c("N", "D"), c("hc.null", "hc.init"))
  colnames(best.param.array) <- c("D", c("hc.null.restart", "hc.init.restart"))
  colnames(mean.run.time) <- c("D", c("hc.null", "hc.init"))
  
  d.num <- 1
  for (d in D.range){
    mean.shd.array[d.num, 1] <- N; sd.shd.array[d.num, 1] <- N
    mean.shd.array[d.num, 2] <- d; sd.shd.array[d.num, 2] <- d
    best.param.array[d.num, 1] <- d; mean.run.time[d.num, 1] <- d
    
    # Run for hc.null
    hc.null <- hc.hyperparam.iteration(N, d, S, sparsity, restart.hyperparams, start=NULL)
    mean.shd.array[d.num, 3] <- hc.null$mean
    sd.shd.array[d.num, 3] <- hc.null$std
    best.param.array[d.num, 2] <- hc.null$best.hyperparam
    mean.run.time[d.num, 2] <- hc.null$av.run.time  
    
    # Run for hc.random
    hc.random <- hc.hyperparam.iteration(N, d, S, sparsity, restart.hyperparams, start="random")
    mean.shd.array[d.num, 4] <- hc.random$mean
    sd.shd.array[d.num, 4] <- hc.random$std
    best.param.array[d.num, 3] <- hc.random$best.hyperparam
    mean.run.time[d.num, 3] <- hc.random$av.run.time  
    
    d.num <- d.num + 1
  } 
  return(list('mean'=mean.shd.array, 'std'=sd.shd.array, 'best.params'=best.param.array, "times"=mean.run.time))
}

hc.learning.iteration <- function(N, d, s, sparsity, restart.hyperparam, start=NULL){
  sim.dag <- generate.dag(N, d, sparsity, s)
  sim.data.df <- as.data.frame(sim.dag$data)
  start.time <- Sys.time()
  print(sprintf(
    'D: %s. Sparsity: %s. Seed: %s. Time: %s',
    d, sparsity, s, start.time
  ))
  if (!is.null(start)){
    set.seed(s)
    initial.graph <- random.graph(nodes(sim.dag$graph), 1, method='ordered')
  } else initial.graph <- random.graph(nodes(sim.dag$graph), 1, method='empty')
  hc.dag <- hc(
    sim.data.df, 
    score='bic-g',
    restart=restart.hyperparam,
    perturb=d,
    start=initial.graph
  )
  run.time <- Sys.time() - start.time
  shd <- shd(cpdag(hc.dag), cpdag(sim.dag$graph))
  return(list("shd"=shd, "run.time"=run.time))
}

hc.hyperparam.iteration <- function(N, d, S, sparsity, restart.hyperparams, start=NULL){
  #'
  best.shd.mean <- 1e9
  # Multiply hc frac by number of nodes
  shd.vals <- c()
  run.times <- c()
  for (restart.hyperparam in restart.hyperparams){
    for (s in 1:S){
      hc.results <- hc.learning.iteration(N, d, s, sparsity, restart.hyperparam, start)
      shd.vals <- c(shd.vals, hc.results$shd)
      run.times <- c(run.times, hc.results$run.time)
    }
    if (mean(shd.vals) < best.shd.mean) {
      best.shd.mean <- mean(shd.vals)
      best.shd.std <- sd(shd.vals)
      av.run.time <- mean(run.times)
      best.hyperparam <- restart.hyperparam
    }
  }
  return(list(
    'mean'=best.shd.mean, 
    'std'=best.shd.std,
    'av.run.time'=av.run.time,
    'best.hyperparam'=best.hyperparam
  ))
}

alpha.constraint.xval.ranDAG <- function(
  functions.bn, alpha.range, n.runs, d.range, n.range, gamma=0
){
  #' For `n.runs` iterations, sample D (the number of graph nodes) 
  #' uniformally between the two values in the list d.range, and sample 
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
  shd.array <- array(numeric(), c(0, num.algos + 4))
  # Moralized Structural Hamming Distance
  hd.array <- array(numeric(), c(0, num.algos + 4))
  # BIC Criterion
  bic.array <- array(numeric(), c(0, num.algos + 4))
  # False Discovery Rate
  fdr.array <- array(numeric(), c(0, num.algos + 4))
  # Sensitivity
  sensitivity.array <- array(numeric(), c(0, num.algos + 4))

  colnames(shd.array) <- c(c("N", "D", "sparsity", "alpha"), names(functions.bn))
  colnames(hd.array) <- c(c("N", "D", "sparsity", "alpha"), names(functions.bn))
  colnames(bic.array) <- c(c("N", "D", "sparsity", "alpha"), names(functions.bn))
  colnames(fdr.array) <- c(c("N", "D", "sparsity", "alpha"), names(functions.bn))
  colnames(sensitivity.array) <- c(c("N", "D", "sparsity", "alpha"), names(functions.bn))

  for (n in seq(1, n.runs)){
    set.seed(n)
    D <- sample(d.range[1]:d.range[2], 1)
    set.seed(n)
    N <- sample(n.range[1]:n.range[2], 1)
    sparsity = 3/(D-1)
    dag <- generate.dag(N, D, sparsity, seed=n)
    true.graph <- dag$graph; sim.data <- dag$data
    for (a in alpha.range){
      row <- c(N, D, sparsity)
      sim.data.df <- as.data.frame(sim.data)    
      row <- c(row, a)
      row.shd <- row
      row.hd <- row
      bic.row <- row
      fdr.row <- row
      sensitivity.row <- row
      for (algo in 1:length(functions.bn)){
        algo.name <- names(functions.bn)[algo]
        start.time <- Sys.time()
        print(sprintf(
          'Run number: %s/%s. Algo: %s Alpha Val: %s. Time: %s',
          n, n.runs, algo.name, a, start.time
        ))       
        if (substr(algo.name, 1, 6) == "hybrid"){
          if (grepl('tabu', algo.name, fixed=TRUE)){
            learned.dag <- functions.bn[[algo]](
              sim.data.df,
              restrict.args=list(alpha=a),
              maximize.args=list(tabu=D)
            )  
          } else{
            learned.dag <- functions.bn[[algo]](
              sim.data.df,
              restrict.args=list(alpha=a)
            )  
          }
        } else learned.dag <- functions.bn[[algo]](sim.data.df, alpha=a)  
        shd <- bnlearn::shd(cpdag(learned.dag), cpdag(true.graph))
        row.shd <- c(row.shd, shd)
        hd <- bnlearn::hamming(moral(learned.dag), moral(true.graph))
        row.hd <- c(row.hd, hd)
        
        # bic <- ebic(learned.dag, sim.data.df, gamma=gamma) # Commenting out since requires full rank
        bic <- 999
        bic.row <- c(bic.row, bic)
        performance.metrics <- calculate.performance.statistics(
          cpdag(learned.dag), cpdag(true.graph)
        )
        fdr.row <- c(fdr.row, performance.metrics$fdr)
        sensitivity.row <- c(sensitivity.row, performance.metrics$sensitivity)        
      } # For algo loop
      shd.array <- rbind(shd.array, row.shd)
      hd.array <- rbind(hd.array, row.hd)
      bic.array <- rbind(bic.array, bic.row)
      fdr.array <- rbind(fdr.array, fdr.row)
      sensitivity.array <- rbind(sensitivity.array, sensitivity.row)
    } # FOR alpha.range
  } # FOR n.run
  return(list(
    'shd'=shd.array,
    'hd'=hd.array,
    'bic'=bic.array,
    'fdr'=fdr.array,
    'sensitivity'=sensitivity.array
  ))
}

tabu.xval.ranDAG <- function(tabu.range, n.runs, d.range, n.range, gamma=0){
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

  # CPDAG Structural Hamming Distance
  shd.array <- array(numeric(), c(0, 7))
  # Moralized Structural Hamming Distance
  hd.array <- array(numeric(), c(0, 7))
  # BIC Criterion
  bic.array <- array(numeric(), c(0, 7))
  # False Discovery Rate
  fdr.array <- array(numeric(), c(0, 7))
  # Sensitivity
  sensitivity.array <- array(numeric(), c(0, 7))
  
  colnames(shd.array) <- c(c("N", "D", "sparsity", "graph.init", "tabu.frac", "score", 'shd'))
  colnames(hd.array) <- c(c("N", "D", "sparsity", "graph.init", "tabu.frac", "score", 'hd'))
  colnames(bic.array) <- c(c("N", "D", "sparsity", "graph.init", "tabu.frac", "score", 'bic'))
  colnames(fdr.array) <- c(c("N", "D", "sparsity", "graph.init", "tabu.frac", "score", "fdr"))
  colnames(sensitivity.array) <- c(c("N", "D", "sparsity", "graph.init", "tabu.frac", "score", "sensitivity"))
  
  for (n in seq(1, n.runs)){
    D <- sample(d.range[1]:d.range[2], 1)
    N <- sample(n.range[1]:n.range[2], 1)
    sparsity <- round(3/(D-1), 4)
    dag <- generate.dag(N, D, sparsity, seed=n)
    true.graph <- dag$graph; sim.data <- dag$data
    for (tabu.frac in tabu.range){
      row <- c(N, D, sparsity)
      tabu.length <- ceiling(D * tabu.frac)
      sim.data.df <- as.data.frame(sim.data)    
      for (i in c('random', 'empty')){
        for (s in c('bic-g', 'bge')){
          row <- c(N, D, sparsity, i, tabu.frac, score)
          row.shd <- row
          row.hd <- row
          bic.row <- row
          fdr.row <- row
          sensitivity.row <- row
          start.time <- Sys.time()
          print(sprintf(
            'Run number: %s/%s. Graph Type: %s Tabu length: %s. Score: %s. Time: %s',
            n, n.runs, i, tabu.length, s, start.time
          ))
          if (!is.null(start)){
            initial.graph <- random.graph(nodes(true.graph), 1, method='ordered')
          } else initial.graph <- random.graph(nodes(true.graph), 1, method='empty')
          learned.dag <- tabu(
            sim.data.df,
            start=initial.graph,
            tabu=tabu.length,
            score=s
          )
          shd <- bnlearn::shd(cpdag(learned.dag), cpdag(true.graph))
          row.shd <- c(row.shd, shd)
          hd <- bnlearn::hamming(moral(learned.dag), moral(true.graph))
          row.hd <- c(row.hd, hd)
          # bic <- ebic(learned.dag, sim.data.df, gamma=0)
          bic <- 999
          bic.row <- c(bic.row, bic)
          performance.metrics <- calculate.performance.statistics(
            cpdag(learned.dag), cpdag(true.graph)
          )
          fdr.row <- c(fdr.row, performance.metrics$fdr)
          sensitivity.row <- c(sensitivity.row, performance.metrics$sensitivity)        
  
          shd.array <- rbind(shd.array, row.shd)
          hd.array <- rbind(hd.array, row.hd)
          bic.array <- rbind(bic.array, bic.row)
          fdr.array <- rbind(fdr.array, fdr.row)
          sensitivity.array <- rbind(sensitivity.array, sensitivity.row)
        } # for bic/bge
      } # for empty/random
    } # FOR alpha.range
  } # FOR n.run
  return(list(
    'shd'=shd.array,
    'hd'=hd.array,
    'bic'=bic.array,
    'fdr'=fdr.array,
    'sensitivity'=sensitivity.array
  ))
}

hc.xval.ranDAG <- function(restart.range, perturb.range, n.runs, d.range, n.range, gamma=0){
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
  
  # CPDAG Structural Hamming Distance
  shd.array <- array(numeric(), c(0, 8))
  # Moralized Structural Hamming Distance
  hd.array <- array(numeric(), c(0, 8))
  # BIC Criterion
  bic.array <- array(numeric(), c(0, 8))
  # False Discovery Rate
  fdr.array <- array(numeric(), c(0, 8))
  # Sensitivity
  sensitivity.array <- array(numeric(), c(0, 8))
  
  colnames(shd.array) <- c(c("N", "D", "sparsity", "graph.init", "restart.no", "perturb", "score", 'shd'))
  colnames(hd.array) <- c(c("N", "D", "sparsity", "graph.init", "restart.no", "perturb", "score", 'hd'))
  colnames(bic.array) <- c(c("N", "D", "sparsity", "graph.init", "restart.no", "perturb", "score", 'bic'))
  colnames(fdr.array) <- c(c("N", "D", "sparsity", "graph.init", "restart.no", "perturb", "score", "fdr"))
  colnames(sensitivity.array) <- c(c("N", "D", "sparsity", "graph.init", "restart.no", "perturb", "score", "sensitivity"))
  
  for (n in seq(1, n.runs)){
    D <- sample(d.range[1]:d.range[2], 1)
    N <- sample(n.range[1]:n.range[2], 1)
    sparsity <- round(3/(D-1), 4)
    dag <- generate.dag(N, D, sparsity, seed=n)
    true.graph <- dag$graph; sim.data <- dag$data
    sim.data.df <- as.data.frame(sim.data)    
    for (i in c('random', 'empty')){
      for (s in c('bic-g', 'bge')){
        if (i == 'random'){
          initial.graph <- random.graph(nodes(true.graph), 1, method='ordered')
        } else initial.graph <- random.graph(nodes(true.graph), 1, method='empty')
        for (p in perturb.range){
          for (r in restart.range){
            row <- c(N, D, sparsity, i, r, p, s)
            row.shd <- row
            row.hd <- row
            bic.row <- row
            fdr.row <- row
            sensitivity.row <- row
            start.time <- Sys.time()
            print(sprintf(
              'Run number: %s/%s. Graph Type: %s Restart no: %s. Perturb: %s Score: %s. Time: %s',
              n, n.runs, i, r, p,  s, start.time
            ))
            learned.dag <- hc(
              sim.data.df, 
              start=initial.graph, 
              restart=r, 
              perturb=p,
              score=s
            )
            shd <- bnlearn::shd(cpdag(learned.dag), cpdag(true.graph))
            row.shd <- c(row.shd, shd)
            hd <- bnlearn::hamming(moral(learned.dag), moral(true.graph))
            row.hd <- c(row.hd, hd)
            # bic <- ebic(learned.dag, sim.data.df, gamma=0)
            bic <- 999
            bic.row <- c(bic.row, bic)
            performance.metrics <- calculate.performance.statistics(
              cpdag(learned.dag), cpdag(true.graph)
            )
            fdr.row <- c(fdr.row, performance.metrics$fdr)
            sensitivity.row <- c(sensitivity.row, performance.metrics$sensitivity)        
            
            shd.array <- rbind(shd.array, row.shd)
            hd.array <- rbind(hd.array, row.hd)
            bic.array <- rbind(bic.array, bic.row)
            fdr.array <- rbind(fdr.array, fdr.row)
            sensitivity.array <- rbind(sensitivity.array, sensitivity.row)
          } # FOR
        } # FOR
      } # for empty/random
    } # FOR alpha.range
  } # FOR n.run
  return(list(
    'shd'=shd.array,
    'hd'=hd.array,
    'bic'=bic.array,
    'fdr'=fdr.array,
    'sensitivity'=sensitivity.array
  ))
}

notears.xval.ranDAG <- function(l1.range, n.runs, d.range, n.range, gamma=0){
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
  
  # CPDAG Structural Hamming Distance
  shd.array <- array(numeric(), c(0, 6))
  # Moralized Structural Hamming Distance
  hd.array <- array(numeric(), c(0, 6))
  # BIC Criterion
  bic.array <- array(numeric(), c(0, 6))
  # False Discovery Rate
  fdr.array <- array(numeric(), c(0, 6))
  # Sensitivity
  sensitivity.array <- array(numeric(), c(0, 6))
  
  colnames(shd.array) <- c(c("N", "D", "sparsity", "seed", 'l1', 'shd'))
  colnames(hd.array) <- c(c("N", "D", "sparsity", "seed", 'l1', 'hd'))
  colnames(bic.array) <- c(c("N", "D", "sparsity", "seed", "l1", 'bic'))
  colnames(fdr.array) <- c(c("N", "D", "sparsity", "seed", "l1", "fdr"))
  colnames(sensitivity.array) <- c(c("N", "D", "sparsity", "seed", "l1", "sensitivity"))
 
  n <- 1
  shift <- 0  
  while (n < n.runs){
    for (l1 in l1.range){
      set.seed(n+shift)
      D <- sample(d.range[1]:d.range[2], 1)
      set.seed(n+shift)
      N <- sample(n.range[1]:n.range[2], 1)
      sparsity <- round(3/(D-1), 4)
      set.seed(n+shift)
      dag <- generate.dag(N, D, sparsity, seed=n)
      true.graph <- dag$graph; sim.data <- dag$data
      
      sim.data.df <- as.data.frame(sim.data)    
      start.time <- Sys.time()
      row <- c(N, D, sparsity, n+shift, l1)
      row.shd <- row; row.hd <- row; bic.row <- row; fdr.row <- row; sensitivity.row <- row
      print(sprintf(
        'Run number: %s/%s. L1: %s. Time: %s',
        n, n.runs, l1, start.time
      ))
      amat.dag <- notears_linear(
        sim.data, 
        lambda1=l1, 
        loss_type='l2'
      )
      # Check if DAG is acyclic
      if (!is_dag(amat.dag)){
        shift = shift + 1
        print(sprintf("Cycles detected. Shift set to %s", shift))
        break
      } else print('Valid DAG')
      
      colnames(amat.dag) <- colnames(sim.data)
      learned.dag <- convert.amat.to.bn(amat.dag)
      shd <- bnlearn::shd(cpdag(learned.dag), cpdag(true.graph))
      row.shd <- c(row.shd, shd)
      hd <- bnlearn::hamming(moral(learned.dag), moral(true.graph))
      row.hd <- c(row.hd, hd)
      # bic <- ebic(learned.dag, sim.data.df, gamma=0)
      bic <- 999
      bic.row <- c(bic.row, bic)
      performance.metrics <- calculate.performance.statistics(
        cpdag(learned.dag), cpdag(true.graph)
      )
      fdr.row <- c(fdr.row, performance.metrics$fdr)
      sensitivity.row <- c(sensitivity.row, performance.metrics$sensitivity)        
      
      shd.array <- rbind(shd.array, row.shd)
      hd.array <- rbind(hd.array, row.hd)
      bic.array <- rbind(bic.array, bic.row)
      fdr.array <- rbind(fdr.array, fdr.row)
      sensitivity.array <- rbind(sensitivity.array, sensitivity.row)
    } # FOR alpha.range
    n = n + 1
  } # FOR n.run
  return(list(
    'shd'=shd.array,
    'hd'=hd.array,
    'bic'=bic.array,
    'fdr'=fdr.array,
    'sensitivity'=sensitivity.array
  ))
}

