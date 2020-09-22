library(bnlearn)
library(mixggm)
library(reticulate)
library(R.utils)
# import('scipy')
use_condaenv("r-reticulate")
source_python('./notears/notears/linear.py')
source_python('./notears/notears/utils.py')

source("./utils.R")


confusion.shd.simulation <- function(algorithms, data.size, run.number, bn.fitted){
  #' Take a fitted BN and generate `n.data.size` samples, and learn
  #' a network across all algorithms in function.list. Repeat
  #' over n.runs. 
  #' 
  #' Returns arrays of the SHD for CPDAG and Moralised graphs.
  shd.array.true.cpdag <- array(
    numeric(),
    c(run.number, length(algorithms)),
  )  
  shd.array.true.moral <- array(
    numeric(),
    c(run.number, length(algorithms)),
  )
  colnames(shd.array.true.cpdag) <- names(algorithms)
  colnames(shd.array.true.moral) <- names(algorithms)
  
  shd.array.rel.cpdag <- array(
    numeric(),
    c(length(algorithms), length(algorithms), run.number),
    dimnames = list(
      names(algorithms),
      names(algorithms),
      paste("run_", 1:run.number, sep='')
    )
  )
  shd.array.rel.moral <- array(
    numeric(),
    c(length(algorithms), length(algorithms), run.number),
    dimnames = list(
      names(algorithms),
      names(algorithms),
      paste("run_", 1:run.number, sep='')
    )
  )
  r <- 1
  shift <- 0
  while (r <= run.number){
    algo.graphs <- vector('list', length(algorithms))
    set.seed(r+shift)
    sim.data <- rbn(bn.fitted, n=data.size)
    for (algo in seq(1, length(algorithms), 1)){
      amat.dag <- NULL
      if (grepl('notears', names(algorithms)[algo], fixed=TRUE)){
        amat.dag <- algorithms[[algo]](as.matrix(sim.data))
        colnames(amat.dag) <- colnames(sim.data)
        if (!is_dag(amat.dag)) {
          shift <- shift + 1
          break
        }
        fitted.network <- convert.amat.to.bn(amat.dag)
      } else {
        fitted.network <- algorithms[[algo]](sim.data)  
      }
      if (!is.null(amat.dag)){
        if (!is_dag(amat.dag)) break
      }
      shd.array.true.cpdag[r, algo] <- bnlearn::shd(fitted.network, bn.net(bn.fitted))
      shd.array.true.cpdag[r, algo] <- bnlearn::shd(cpdag(fitted.network), cpdag(bn.net(bn.fitted)))
      shd.array.true.moral[r, algo] <- bnlearn::shd(moral(fitted.network), moral(bn.net(bn.fitted)))
      algo.graphs[[algo]] <- fitted.network
      print(sprintf(
        '%s - %s - Run: %s/%s - CPDAG SHD: %s - Moralised SHD: %s',
        Sys.time(), 
        names(algorithms)[algo], 
        r, 
        run.number,
        shd.array.true.cpdag[r, algo],
        shd.array.true.moral[r, algo]
      ))    
    }
    if (is_dag(amat.dag)){
      for (i in seq(1, length(algo.graphs), 1)){
        for (j in seq(1, length(algo.graphs), 1)){
          shd.array.rel.cpdag[i, j, r] = bnlearn::shd(cpdag(algo.graphs[[i]]), cpdag(algo.graphs[[j]]))
          shd.array.rel.moral[i, j, r] = bnlearn::shd(moral(algo.graphs[[i]]), moral(algo.graphs[[j]]))
        }
      }
      r <- r + 1
    }
  }
  return(list(
    'shd.true.cpdag'=data.frame(shd.array.true.cpdag),
    'shd.true.moral'=data.frame(shd.array.true.moral),
    'shd.relative.cpdag'=shd.array.rel.cpdag,
    'shd.relative.moral'=shd.array.rel.moral
  ))
}

data.size.study.network <- function(algorithms, data.size.ranges, n.runs, network){
  #' Take a fitted BN and generate `n.data.size` samples, and learn
  #' a network across all algorithms in function.list. Repeat
  #' over n.runs. 
  #' 
  #' Returns arrays of the SHD for CPDAG and Moralised graphs.
  num.algos <- length(algorithms)
  D <- length(nodes(network))
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
  for (data.size in data.size.ranges){
    shift <- 10
    for (r in seq(1, n.runs, 1)){
      algo.graphs <- vector('list', length(algorithms))
      set.seed(r+shift)
      sim.data <- rbn(network, n=data.size)
      row <- c(r, data.size, D)
      row.shd <- row;  row.hd <- row;  bic.row <- row; fdr.row <- row
      sensitivity.row <- row
      for (algo in seq(1, length(algorithms), 1)){
        start.time <- Sys.time()
        if (grepl('notears|lingam', names(algorithms)[algo], fixed=FALSE)){
          amat.dag <- algorithms[[algo]](as.matrix(sim.data))
          while (!is_dag(amat.dag)){
            shift <- shift + 1
            print(sprintf('Shift changed to %s. Time: %s', shift, Sys.time()))
            sim.data <- rbn(network, n=data.size)
            amat.dag <- algorithms[[algo]](as.matrix(sim.data))
          }
          colnames(amat.dag) <- colnames(as.matrix(sim.data))
          fitted.network <- convert.amat.to.bn(amat.dag)
        } else if (grepl('tabu', names(algorithms)[algo], fixed=TRUE)){
          fitted.network <- algorithms[[algo]](
            sim.data,
            maximize.args=list('tabu'=D, 'score'='bge')
          )
          while (!acyclic(fitted.network, directed=TRUE)){
            shift <- shift + 1
            print(sprintf('Shift changed to %s. Time: %s', shift, Sys.time()))
            sim.data <- rbn(network, n=data.size)
            fitted.network <- algorithms[[algo]](sim.data)
          }
        } else {
          fitted.network <- algorithms[[algo]](sim.data)
          while (!acyclic(fitted.network, directed=TRUE)){
            shift <- shift + 1
            print(sprintf('Shift changed to %s. Time: %s', shift, Sys.time()))
            sim.data <- rbn(network, n=data.size)
            fitted.network <- algorithms[[algo]](sim.data)
          }
        }
        
        shd <- bnlearn::shd(cpdag(fitted.network), cpdag(network))
        row.shd <- c(row.shd, shd)
        print(sprintf(
          'Run number: %s/%s. Algo: %s. Data size: %s. SHD: %s. Time: %s',
          r, n.runs, names(algorithms)[algo], data.size, shd, start.time
        ))       
        
        hd <- bnlearn::hamming(moral(fitted.network), moral(network))
        row.hd <- c(row.hd, hd)
        performance.metrics <- calculate.performance.statistics(
          cpdag(fitted.network), cpdag(network)
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
  return(list(
    'shd'=data.frame(shd.array),
    'hd'=data.frame(hd.array),
    'fdr'=fdr.array,
    'sensitivity'=sensitivity.array
  ))
}

alpha.xval.network <- function(functions.bn, alpha.range, n.runs, n.range, network, gamma=0){
  #' For `n.runs` iterations, sample `N` (the number of datapoints) points uniformally
  #' from a bn.fit object `network` . 
  #' For each algorithm in functions.bn, learn a graph for each alpha within
  #' alpha.range, and calculate the 
  #' * CPDAG SHD
  #' * moralised SHD 
  #' * BIC criterion for the moralised graph. 
  #' * CPDAG False Discovery Rate (FDR)
  #' * CPDAG Sensitivity
  
  num.algos <- length(functions.bn)
  D <- length(nodes(network))
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
  
  colnames(shd.array) <- c(c("seed", "N", "D", "alpha"), names(functions.bn))
  colnames(hd.array) <- c(c("seed", "N", "D", "alpha"), names(functions.bn))
  colnames(bic.array) <- c(c("seed", "N", "D", "alpha"), names(functions.bn))
  colnames(fdr.array) <- c(c("seed", "N", "D", "alpha"), names(functions.bn))
  colnames(sensitivity.array) <- c(c("seed", "N", "D", "alpha"), names(functions.bn))
  
  n <- 1
  shift <- 1
  while (n <= n.runs){
    set.seed(n+shift)
    N <- sample(n.range[1]:n.range[2], 1)
    set.seed(n)
    sim.data <- rbn(network, n=N)
    for (a in alpha.range){
      sim.data.df <- as.data.frame(sim.data)    
      row <- c(n+shift, N, D, a)
      row.shd <- row
      row.hd <- row
      bic.row <- row
      fdr.row <- row
      sensitivity.row <- row
      for (algo in 1:length(functions.bn)){
        start.time <- Sys.time()
        algo.name <-  names(functions.bn)[algo]
        print(sprintf(
          'Run number: %s/%s. Shift: %s. Algo: %s Alpha Val: %s. Time: %s',
          n, n.runs, shift, algo.name, a, start.time
        ))       
        learned.dag <- NULL
        if (substr(algo.name, 1, 6) == "hybrid"){
          if (grepl('tabu', algo.name, fixed=TRUE)){
            learned.dag <- withTimeout(functions.bn[[algo]](
              sim.data.df,
              restrict.args=list(alpha=a, test='cor'),
              maximize.args=list(tabu=D)
            ), timeout=900, onTimeout=c('silent'))
            } else {
              learned.dag <- withTimeout(functions.bn[[algo]](
                sim.data.df,
                restrict.args=list(alpha=a, test='cor')
              ), timeout=900, onTimeout=c('silent'))
            }
        } else {
          learned.dag <- withTimeout(
            functions.bn[[algo]](sim.data.df, alpha=a),
            timeout=900, onTimeout=c('silent')
          )
        }
        if (is.null(learned.dag)) {
          shift <- shift + 1
          break
        }
        if (!acyclic(learned.dag, directed=TRUE)) {
          shift <- shift + 1
          break
        }
        shd <- bnlearn::shd(cpdag(learned.dag), cpdag(network))
        row.shd <- c(row.shd, shd)
        hd <- bnlearn::hamming(moral(learned.dag), moral(network))
        row.hd <- c(row.hd, hd)
        # bic <- ebic(learned.dag, sim.data.df, gamma=gamma)
        bic <- 9
        bic.row <- c(bic.row, bic)
        performance.metrics <- calculate.performance.statistics(
          cpdag(learned.dag), cpdag(network)
        )
        fdr.row <- c(fdr.row, performance.metrics$fdr)
        sensitivity.row <- c(sensitivity.row, performance.metrics$sensitivity)        
      } # For algo loop
      if (is.null(learned.dag)) break
      if (!acyclic(learned.dag, directed=TRUE)) break
      shd.array <- rbind(shd.array, row.shd)
      hd.array <- rbind(hd.array, row.hd)
      bic.array <- rbind(bic.array, bic.row)
      fdr.array <- rbind(fdr.array, fdr.row)
      sensitivity.array <- rbind(sensitivity.array, sensitivity.row)
    } # FOR alpha.range
    if (!is.null(learned.dag)) {
      if (acyclic(learned.dag, directed=TRUE)){
        n <- n + 1
      }
    }
  } # FOR n.run
  return(list(
    'shd'=shd.array,
    'hd'=hd.array,
    'bic'=bic.array,
    'fdr'=fdr.array,
    'sensitivity'=sensitivity.array
  ))
}

hc.xval.network <- function(restart.range, perturb.range, n.runs, n.range, network, gamma=0.){
  #' For `n.runs` iterations, sample `N` (the number of datapoints) points uniformally
  #' from a bn.fit object `network` . 
  #' For each algorithm in functions.bn, learn a graph for each alpha within
  #' alpha.range, and calculate the 
  #' * CPDAG SHD
  #' * moralised SHD 
  #' * BIC criterion for the moralised graph. 
  #' * CPDAG False Discovery Rate (FDR)
  #' * CPDAG Sensitivity

  D <- length(nodes(network))
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
  
  colnames(shd.array) <- c(c("seed", "N", "D", "graph.init", "restart.no", "perturb", "score", 'shd'))
  colnames(hd.array) <- c(c("seed", "N", "D", "graph.init", "restart.no", "perturb", "score", 'hd'))
  colnames(bic.array) <- c(c("seed", "N", "D", "graph.init", "restart.no", "perturb", "score", 'bic'))
  colnames(fdr.array) <- c(c("seed", "N", "D", "graph.init", "restart.no", "perturb", "score", "fdr"))
  colnames(sensitivity.array) <- c(c("seed", "N", "D", "graph.init", "restart.no", "perturb", "score", "sensitivity"))
  
  n <- 1
  shift <- 0
  while (n <= n.runs){
    N <- sample(n.range[1]:n.range[2], 1)
    set.seed(n+shift)
    sim.data <- rbn(network, n=N)
    for (r in restart.range){
      for (p in perturb.range){
        sim.data.df <- as.data.frame(sim.data)    
        start.time <- Sys.time()
        for (i in c('random', 'empty')){
          for (s in c('bic-g', 'bge')){
            row <- c(n+shift, N, D, i, r, p, s)
            row.shd <- row; row.hd <- row; bic.row <- row; fdr.row <- row; sensitivity.row <- row
            print(sprintf(
              'Run number: %s/%s. Shift: %s. Graph Type: %s Restart number: %s. Perturbations: %s Time: %s',
              n, n.runs, shift, i, r, p, start.time
            ))
            if (i == "random"){
              initial.graph <- random.graph(nodes(network), 1, method='ordered')
            } else initial.graph <- random.graph(nodes(network), 1, method='empty')
            learned.dag <- NULL
            learned.dag <- withTimeout(hc(
              sim.data.df, 
              start=initial.graph, 
              restart=r, 
              perturb=p,
              score=s
            ), timeout=300, onTimeout=c('silent'))
            if (is.null(learned.dag)) {
              shift <- shift + 1
              break
            }
            shd <- bnlearn::shd(cpdag(learned.dag), cpdag(network))
            row.shd <- c(row.shd, shd)
            hd <- bnlearn::hamming(moral(learned.dag), moral(network))
            row.hd <- c(row.hd, hd)
            # bic <- ebic(learned.dag, sim.data.df, gamma=gamma)
            bic <- 9
            bic.row <- c(bic.row, bic)
            performance.metrics <- calculate.performance.statistics(
              cpdag(learned.dag), cpdag(network)
            )
            fdr.row <- c(fdr.row, performance.metrics$fdr)
            sensitivity.row <- c(sensitivity.row, performance.metrics$sensitivity)        
            
            shd.array <- rbind(shd.array, row.shd)
            hd.array <- rbind(hd.array, row.hd)
            bic.array <- rbind(bic.array, bic.row)
            fdr.array <- rbind(fdr.array, fdr.row)
            sensitivity.array <- rbind(sensitivity.array, sensitivity.row)
          }
          if (is.null(learned.dag)) break
        } # FOR score
        if (is.null(learned.dag)) break
      }
      if (is.null(learned.dag)) break
    } # FOR alpha.range
    if (!is.null(learned.dag)) n <- n + 1
  } # FOR n.run
  return(list(
    'shd'=shd.array,
    'hd'=hd.array,
    'bic'=bic.array,
    'fdr'=fdr.array,
    'sensitivity'=sensitivity.array
  ))
}

tabu.xval.network <- function(tabu.range, n.runs, n.range, network, gamma=0.){
  #' For `n.runs` iterations, sample `N` (the number of datapoints) points uniformally
  #' from a bn.fit object `network` . 
  #' For each algorithm in functions.bn, learn a graph for each alpha within
  #' alpha.range, and calculate the 
  #' * CPDAG SHD
  #' * moralised SHD 
  #' * BIC criterion for the moralised graph. 
  #' * CPDAG False Discovery Rate (FDR)
  #' * CPDAG Sensitivity

  D <- length(nodes(network))
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
  
  colnames(shd.array) <- c(c("N", "D", "graph.init", "tabu.frac", "score", 'shd'))
  colnames(hd.array) <- c(c("N", "D", "graph.init", "tabu.frac", "score", 'hd'))
  colnames(bic.array) <- c(c("N", "D", "graph.init", "tabu.frac", "score", 'bic'))
  colnames(fdr.array) <- c(c("N", "D", "graph.init", "tabu.frac", "score", "fdr"))
  colnames(sensitivity.array) <- c(c("N", "D", "graph.init", "tabu.frac", "score", "sensitivity"))
  
  for (n in seq(1, n.runs)){
    N <- sample(n.range[1]:n.range[2], 1)
    set.seed(n)
    sim.data <- rbn(network, n=N)
    for (tabu.frac in tabu.range){
      tabu.length <- ceiling(D * tabu.frac)
      sim.data.df <- as.data.frame(sim.data)    
      for (i in c('random', 'empty')){
        for (s in c('bic-g', 'bge')){
          row <- c(N, D, i, tabu.frac, s)
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
            initial.graph <- random.graph(nodes(network), 1, method='ordered')
          } else initial.graph <- random.graph(nodes(network), 1, method='empty')
          learned.dag <- tabu(
            sim.data.df, 
            start=initial.graph, 
            tabu=tabu.length,
            score=s
          )
          shd <- bnlearn::shd(cpdag(learned.dag), cpdag(network))
          row.shd <- c(row.shd, shd)
          hd <- bnlearn::hamming(moral(learned.dag), moral(network))
          row.hd <- c(row.hd, hd)
          # bic <- ebic(learned.dag, sim.data.df, gamma=gamma)
          bic <- 9
          bic.row <- c(bic.row, bic)
          performance.metrics <- calculate.performance.statistics(
            cpdag(learned.dag), cpdag(network)
          )
          fdr.row <- c(fdr.row, performance.metrics$fdr)
          sensitivity.row <- c(sensitivity.row, performance.metrics$sensitivity)        
          
          shd.array <- rbind(shd.array, row.shd)
          hd.array <- rbind(hd.array, row.hd)
          bic.array <- rbind(bic.array, bic.row)
          fdr.array <- rbind(fdr.array, fdr.row)
          sensitivity.array <- rbind(sensitivity.array, sensitivity.row)
        } # FOR score
      }
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

notears.xval.network <- function(l1.range, n.runs, n.range, network, gamma=0.){
  #' For `n.runs` iterations, sample `N` (the number of datapoints) points uniformally
  #' from a bn.fit object `network` . 
  #' For each algorithm in functions.bn, learn a graph for each alpha within
  #' alpha.range, and calculate the 
  #' * CPDAG SHD
  #' * moralised SHD 
  #' * BIC criterion for the moralised graph. 
  #' * CPDAG False Discovery Rate (FDR)
  #' * CPDAG Sensitivity

  D <- length(nodes(network))
  # CPDAG Structural Hamming Distance
  shd.array <- array(numeric(), c(0, 5))
  # Moralized Structural Hamming Distance
  hd.array <- array(numeric(), c(0, 5))
  # BIC Criterion
  bic.array <- array(numeric(), c(0, 5))
  # False Discovery Rate
  fdr.array <- array(numeric(), c(0, 5))
  # Sensitivity
  sensitivity.array <- array(numeric(), c(0, 5))
  
  colnames(shd.array) <- c(c("N", "D", 'seed', "l1", 'shd'))
  colnames(hd.array) <- c(c("N", "D", 'seed', "l1", 'hd'))
  colnames(bic.array) <- c(c("N", "D", 'seed', "l1", 'bic'))
  colnames(fdr.array) <- c(c("N", "D", 'seed', "l1", "fdr"))
  colnames(sensitivity.array) <- c(c("N", "D", 'seed', "l1", "sensitivity"))
  n <- 1
  shift <- 130
  while (n <= n.runs){
    for (l1 in l1.range){
      set.seed(n + shift)
      N <- sample(n.range[1]:n.range[2], 1)
      set.seed(n + shift)
      sim.data.df <- rbn(network, n=N)
      sim.data <- as.matrix(sim.data.df)    
  
      row <- c(N, D, n+shift, l1)
      row.shd <- row
      row.hd <- row
      bic.row <- row
      fdr.row <- row
      sensitivity.row <- row
      start.time <- Sys.time()
      print(sprintf(
        'Run number: %s/%s. Lambda1: %s. Time: %s',
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
      
      shd <- bnlearn::shd(cpdag(learned.dag), cpdag(network))
      row.shd <- c(row.shd, shd)
      hd <- bnlearn::hamming(moral(learned.dag), moral(network))
      row.hd <- c(row.hd, hd)
      # bic <- ebic(bn.fit(learned.dag, sim.data.df), sim.data.df, gamma=gamma)
      bic <- 9
      bic.row <- c(bic.row, bic)
      performance.metrics <- calculate.performance.statistics(
        cpdag(learned.dag), cpdag(network)
      )
      fdr.row <- c(fdr.row, performance.metrics$fdr)
      sensitivity.row <- c(sensitivity.row, performance.metrics$sensitivity)        
      
      shd.array <- rbind(shd.array, row.shd)
      hd.array <- rbind(hd.array, row.hd)
      bic.array <- rbind(bic.array, bic.row)
      fdr.array <- rbind(fdr.array, fdr.row)
      sensitivity.array <- rbind(sensitivity.array, sensitivity.row)
    } # FOR alpha.range
    if (is_dag(amat.dag)){
      n = n + 1
    }
  } # while n < n.runs
  return(list(
    'shd'=shd.array,
    'hd'=hd.array,
    'bic'=bic.array,
    'fdr'=fdr.array,
    'sensitivity'=sensitivity.array
  ))
}

alpha.network.search <- function(functions.bn, alpha.range, data.size.range, n.runs, bn.fitted){
    #' For datasizes N in the vector `data.size.range`, generate 
    #' `n.runs` datasets from the network `bn.fitted`. Calculate the 
    #' SHD, FDR and sensitivity of learned graphs to the true network.
  num.algos <- length(functions.bn)
  
  # CPDAG Structural Hamming Distance
  shd.array <- array(numeric(), c(0, num.algos + 2))
  # Moralized Structural Hamming Distance
  hd.array <- array(numeric(), c(0, num.algos + 2))
  # BIC Criterion
  bic.array <- array(numeric(), c(0, num.algos + 2))
  # False Discovery Rate
  fdr.array <- array(numeric(), c(0, num.algos + 2))
  # Sensitivity
  sensitivity.array <- array(numeric(), c(0, num.algos + 2))
  
  colnames(shd.array) <- c(c("N", "alpha"), names(functions.bn))
  colnames(hd.array) <- c(c("N", "alpha"), names(functions.bn))
  colnames(bic.array) <- c(c("N", "alpha"), names(functions.bn))
  colnames(fdr.array) <- c(c("N", "alpha"), names(functions.bn))
  colnames(sensitivity.array) <- c(c("N", "alpha"), names(functions.bn))
  for (data.size in data.size.range){
    for (run in seq(1, n.run)){
      set.seed(run)
      sim.data <- rbn(bn.fitted, n=data.size)
      for (alpha in alpha.range){
        row <- c(data.size, alpha)
        row.shd <- row
        row.hd <- row
        bic.row <- row
        fdr.row <- row
        sensitivity.row <- row
        for (algo in 1:num.algos){
          learned.dag <- functions.bn[[algo]](sim.data.df, alpha=a)  
          shd <- bnlearn:shd(cpdag(learned.dag), cpdag(bn.fitted))
          row.shd <- c(row.shd, shd)
          hd <- bnlearn::hamming(cpdag(learned.dag), cpdag(bn.fitted))
          row.hd <- c(row.hd, hd)
          bic <- ebic(learned.dag, sim.data.df, gamma=gamma)
          bic.row <- c(bic.row, bic)
          performance.metrics <- calculate.performance.statistics(
            cpdag(learned.dag), cpdag(bn.fitted)
          )
          fdr.row <- c(fdr.row, performance.metrics$fdr)
          sensitivity.row <- c(sensitivity.row, performance.metrics$sensitivity)
        } # for algo
        shd.array <- rbind(shd.array, row.shd)
        hd.array <- rbind(hd.array, row.hd)
        bic.array <- rbind(bic.array, bic.row)
        fdr.array <- rbind(fdr.array, fdr.row)
        sensitivity.array <- rbind(sensitivity.array, sensitivity.row)
      } # for alpha 
    } # for run
  } # for data.size
  return(list(
    'cpdag'=shd.array,
    'moral'=hd.array,
    'bic'=bic.array,
    'fdr'=fdr.array,
    'sensitivity'=sensitivity.array
  ))
}
  
network.search <- function(functions.bn, alpha.range, data.size.range, n.runs, bn.fitted){
  #' For datasizes N in the vector `data.size.range`, generate 
  #' `n.runs` datasets from the network `bn.fitted`. Calculate the 
  #' SHD, FDR and sensitivity of learned graphs to the true network.
  num.algos <- length(functions.bn)
  
  # CPDAG Structural Hamming Distance
  shd.array <- array(numeric(), c(0, num.algos + 2))
  # Moralized Structural Hamming Distance
  hd.array <- array(numeric(), c(0, num.algos + 2))
  # BIC Criterion
  bic.array <- array(numeric(), c(0, num.algos + 2))
  # False Discovery Rate
  fdr.array <- array(numeric(), c(0, num.algos + 2))
  # Sensitivity
  sensitivity.array <- array(numeric(), c(0, num.algos + 2))
  
  colnames(shd.array) <- c(c("N", "alpha"), names(functions.bn))
  colnames(hd.array) <- c(c("N", "alpha"), names(functions.bn))
  colnames(bic.array) <- c(c("N", "alpha"), names(functions.bn))
  colnames(fdr.array) <- c(c("N", "alpha"), names(functions.bn))
  colnames(sensitivity.array) <- c(c("N", "alpha"), names(functions.bn))
  for (data.size in data.size.range){
    for (run in seq(1, n.run)){
      set.seed(run)
      sim.data <- rbn(bn.fitted, n=data.size)
      row <- c(data.size)
      row.shd <- row
      row.hd <- row
      bic.row <- row
      fdr.row <- row
      sensitivity.row <- row
      for (algo in 1:num.algos){
        if (substr(algo.name, 1, 7) == "notears"){
          amat.dag <- functions.bn[[algo]](sim.data)
          # Check if DAG is acyclic
          if (!is_dag(amat.dag)){
            shift = shift + 1
            print(sprintf("Cycles detected. Shift set to %s", shift))
            break
          } else print('Valid DAG')
          colnames(amat.dag) <- colnames(sim.data)
          learned.dag <- convert.amat.to.bn(amat.dag)
        } else {
          learned.dag <- functions.bn[[algo]](sim.data.df)  
        }
        if (!is_dag(amat.dag)) break
        
        shd <- shd(cpdag(learned.dag), cpdag(bn.fitted))
        row.shd <- c(row.shd, shd)
        hd <- shd(moral(learned.dag), moral(bn.fitted))
        row.hd <- c(row.hd, hd)
        bic <- ebic(learned.dag, sim.data.df, gamma=gamma)
        bic.row <- c(bic.row, bic)
        performance.metrics <- calculate.performance.statistics(
          cpdag(learned.dag), cpdag(bn.fitted)
        )
        fdr.row <- c(fdr.row, performance.metrics$fdr)
        sensitivity.row <- c(sensitivity.row, performance.metrics$sensitivity)
      } # for algo
      shd.array <- rbind(shd.array, row.shd)
      hd.array <- rbind(hd.array, row.hd)
      bic.array <- rbind(bic.array, bic.row)
      fdr.array <- rbind(fdr.array, fdr.row)
      sensitivity.array <- rbind(sensitivity.array, sensitivity.row)
    } # for run
  } # for data.size
  return(list(
    'cpdag'=shd.array,
    'moral'=hd.array,
    'bic'=bic.array,
    'fdr'=fdr.array,
    'sensitivity'=sensitivity.array
  ))
}


# pc.f <- partial(pc.stable, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
# gs.f <- partial(gs, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
# iamb.f <- partial(iamb, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
# fast.iamb.f <- partial(fast.iamb, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
# inter.iamb.f <- partial(inter.iamb, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
# iamb.fdr.f <- partial(iamb.fdr, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
## mmpc.f <- partial(mmpc, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
## si.hiton.pc.f <- partial(si.hiton.pc, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
## hpc.f <- partial(hpc, test = 'mi-g-sh', undirected = FALSE, debug=FALSE, alpha=alpha)
# tabu.f <- partial(tabu, score = 'bic-g', tabu=20)
# hc.f <- partial(hc, score = 'bic-g', restart=20, perturb=10)
# 

# 
# algorithms <- list(
#   "pc"=pc.f,
#   "gs"=gs.f,
#   "iamb"=iamb.f,
#   "fast.iamb"=fast.iamb.f,
#   "inter.iamb"= inter.iamb.f,
#   "iamb.fdr"=inter.iamb.f,
#   # "mmpc"=mmpc.f,
#   "si.hiton.pc"=si.hiton.pc.f,
#   "hpc"=hpc.f
# )
# 
# N <- 10000
# n.runs <- 10
# ecoli.run = confusion.shd.simulation(algorithms, N, n.runs, ecoli)
# confusion.save.template <- "data/results/simulated_networks/ecoli/confusion_matrices/%s/%s_N_%s_run_%s.csv"
# write.csv(
#   apply(ecoli.run$shd.relative.cpdag, c(1,2), mean), 
#   sprintf(confusion.save.template, 'cpdag', "mean", N, n.runs), 
#   row.names = FALSE
# )
# write.csv(
#   apply(ecoli.run$shd.relative.moral, c(1,2), mean), 
#   sprintf(confusion.save.template, 'moral', "mean", N, n.runs), 
#   row.names = FALSE
# )
# 
# shd.runs.template <- "data/results/simulated_networks/ecoli/runs/%s/%s_N_%s_run_%s.csv"
# write.csv(
#   apply(ecoli.run$shd.true.cpdag, c(1,2), mean), 
#   sprintf(shd.runs.template, 'cpdag', "mean", N, n.runs), 
#   row.names = FALSE
# )
# write.csv(
#   apply(ecoli.run$shd.true.moral, c(1,2), mean), 
#   sprintf(shd.runs.template, 'moral', "mean", N, n.runs), 
#   row.names = FALSE
# )
# 
# # # z = melt(z, value.name='shd')
# # ggplot(z, aes(x=variable, y=shd)) + geom_violin()