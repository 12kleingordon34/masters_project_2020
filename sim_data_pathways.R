library(bnlearn)
library(mixggm)
library(reticulate)
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
  
  for (r in seq(1, run.number, 1)){
    algo.graphs <- vector('list', length(algorithms))
    set.seed(r)
    sim.data <- rbn(bn.fitted, n=data.size)
    for (algo in seq(1, length(algorithms), 1)){
      fitted.network <- algorithms[[algo]](sim.data)
      shd.array.true.cpdag[r, algo] <- shd(cpdag(fitted.network), cpdag(bn.fitted))
      shd.array.true.moral[r, algo] <- shd(moral(fitted.network), moral(bn.fitted))
      algo.graphs[[algo]] <- fitted.network
      print(sprintf(
        '%s - %s - Run: %s - CPDAG SHD: %s - Moralised SHD: %s',
        Sys.time(), 
        names(algorithms)[algo], 
        r, 
        shd.array.true.cpdag[r, algo],
        shd.array.true.moral[r, algo]
      ))    
    }

    for (i in seq(1, length(algo.graphs), 1)){
      for (j in seq(1, length(algo.graphs), 1)){
        shd.array.rel.cpdag[i, j, r] = shd(cpdag(algo.graphs[[i]]), cpdag(algo.graphs[[j]]))
        shd.array.rel.moral[i, j, r] = shd(moral(algo.graphs[[i]]), moral(algo.graphs[[j]]))
      }
    }
  }
  return(list(
    'shd.true.cpdag'=data.frame(shd.array.true.cpdag),
    'shd.true.moral'=data.frame(shd.array.true.moral),
    'shd.relative.cpdag'=shd.array.rel.cpdag,
    'shd.relative.moral'=shd.array.rel.moral
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
  cpdag.shd.array <- array(numeric(), c(0, num.algos + 3))
  # Moralized Structural Hamming Distance
  moral.shd.array <- array(numeric(), c(0, num.algos + 3))
  # BIC Criterion
  bic.array <- array(numeric(), c(0, num.algos + 3))
  # False Discovery Rate
  fdr.array <- array(numeric(), c(0, num.algos + 3))
  # Sensitivity
  sensitivity.array <- array(numeric(), c(0, num.algos + 3))
  
  colnames(cpdag.shd.array) <- c(c("N", "D", "alpha"), names(functions.bn))
  colnames(moral.shd.array) <- c(c("N", "D", "alpha"), names(functions.bn))
  colnames(bic.array) <- c(c("N", "D", "alpha"), names(functions.bn))
  colnames(fdr.array) <- c(c("N", "D", "alpha"), names(functions.bn))
  colnames(sensitivity.array) <- c(c("N", "D", "alpha"), names(functions.bn))
  
  for (n in seq(1, n.runs)){
    set.seed(n)
    N <- sample(n.range[1]:n.range[2], 1)
    set.seed(n)
    sim.data <- rbn(network, n=N)
    for (a in alpha.range){
      sim.data.df <- as.data.frame(sim.data)    
      start.time <- Sys.time()
      row <- c(N, D, a)
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
        cpdag.shd <- bnlearn::shd(cpdag(learned.dag), cpdag(network))
        row.cpdag.shd <- c(row.cpdag.shd, cpdag.shd)
        moral.shd <- bnlearn::shd(moral(learned.dag), moral(network))
        row.moral.shd <- c(row.moral.shd, moral.shd)
        bic <- ebic(learned.dag, sim.data.df, gamma=gamma)
        bic.row <- c(bic.row, bic)
        performance.metrics <- calculate.performance.statistics(
          cpdag(learned.dag), cpdag(network)
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

hc.xval.network <- function(restart.range, n.runs, n.range, network, gamma=0.){
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
  cpdag.shd.array <- array(numeric(), c(0, 5))
  # Moralized Structural Hamming Distance
  moral.shd.array <- array(numeric(), c(0, 5))
  # BIC Criterion
  bic.array <- array(numeric(), c(0, 5))
  # False Discovery Rate
  fdr.array <- array(numeric(), c(0, 5))
  # Sensitivity
  sensitivity.array <- array(numeric(), c(0, 5))
  
  colnames(cpdag.shd.array) <- c(c("N", "D", "graph.init", "restart.no", 'cpdag.shd'))
  colnames(moral.shd.array) <- c(c("N", "D", "graph.init", "restart.no", 'moral.shd'))
  colnames(bic.array) <- c(c("N", "D", "graph.init", "restart.no", 'bic'))
  colnames(fdr.array) <- c(c("N", "D", "graph.init", "restart.no", "fdr"))
  colnames(sensitivity.array) <- c(c("N", "D", "graph.init", "restart.no", "sensitivity"))
  
  for (n in seq(1, n.runs)){
    N <- sample(n.range[1]:n.range[2], 1)
    set.seed(n)
    sim.data <- rbn(network, n=N)
    for (restart.no in restart.range){
      sim.data.df <- as.data.frame(sim.data)    
      start.time <- Sys.time()
      for (i in c('random', 'empty')){
        row <- c(N, D, restart.no, i)
        row.cpdag.shd <- row
        row.moral.shd <- row
        bic.row <- row
        fdr.row <- row
        sensitivity.row <- row
        print(sprintf(
          'Run number: %s/%s. Graph Type: %s Restart number: %s. Time: %s',
          n, n.runs, i, restart.no, start.time
        ))
        if (!is.null(start)){
          initial.graph <- random.graph(nodes(network), 1, method='ordered')
        } else initial.graph <- random.graph(nodes(network), 1, method='empty')
        learned.dag <- hc(
          sim.data.df, 
          start=initial.graph, 
          restart=restart.no, 
          perturb=20
        )
        cpdag.shd <- bnlearn::shd(cpdag(learned.dag), cpdag(network))
        row.cpdag.shd <- c(row.cpdag.shd, cpdag.shd)
        moral.shd <- bnlearn::shd(moral(learned.dag), moral(network))
        row.moral.shd <- c(row.moral.shd, moral.shd)
        bic <- ebic(learned.dag, sim.data.df, gamma=gamma)
        bic.row <- c(bic.row, bic)
        performance.metrics <- calculate.performance.statistics(
          cpdag(learned.dag), cpdag(network)
        )
        fdr.row <- c(fdr.row, performance.metrics$fdr)
        sensitivity.row <- c(sensitivity.row, performance.metrics$sensitivity)        
        
        cpdag.shd.array <- rbind(cpdag.shd.array, row.cpdag.shd)
        moral.shd.array <- rbind(moral.shd.array, row.moral.shd)
        bic.array <- rbind(bic.array, bic.row)
        fdr.array <- rbind(fdr.array, fdr.row)
        sensitivity.array <- rbind(sensitivity.array, sensitivity.row)
      }
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
  cpdag.shd.array <- array(numeric(), c(0, 5))
  # Moralized Structural Hamming Distance
  moral.shd.array <- array(numeric(), c(0, 5))
  # BIC Criterion
  bic.array <- array(numeric(), c(0, 5))
  # False Discovery Rate
  fdr.array <- array(numeric(), c(0, 5))
  # Sensitivity
  sensitivity.array <- array(numeric(), c(0, 5))
  
  colnames(cpdag.shd.array) <- c(c("N", "D", "graph.init", "tabu.frac", 'cpdag.shd'))
  colnames(moral.shd.array) <- c(c("N", "D", "graph.init", "tabu.frac", 'moral.shd'))
  colnames(bic.array) <- c(c("N", "D", "graph.init", "tabu.frac", 'bic'))
  colnames(fdr.array) <- c(c("N", "D", "graph.init", "tabu.frac", "fdr"))
  colnames(sensitivity.array) <- c(c("N", "D", "graph.init", "tabu.frac", "sensitivity"))
  
  for (n in seq(1, n.runs)){
    N <- sample(n.range[1]:n.range[2], 1)
    set.seed(n)
    sim.data <- rbn(network, n=N)
    for (tabu.frac in tabu.range){
      tabu.length <- ceiling(D * tabu.frac)
      sim.data.df <- as.data.frame(sim.data)    
      start.time <- Sys.time()
      for (i in c('random', 'empty')){
        row <- c(N, D, i, tabu.frac)
        row.cpdag.shd <- row
        row.moral.shd <- row
        bic.row <- row
        fdr.row <- row
        sensitivity.row <- row
        print(sprintf(
          'Run number: %s/%s. Graph Type: %s Tabu length: %s. Time: %s',
          n, n.runs, i, tabu.length, start.time
        ))
        if (!is.null(start)){
          initial.graph <- random.graph(nodes(network), 1, method='ordered')
        } else initial.graph <- random.graph(nodes(network), 1, method='empty')
        learned.dag <- tabu(
          sim.data.df, 
          start=initial.graph, 
          tabu=tabu.length
        )
        cpdag.shd <- bnlearn::shd(cpdag(learned.dag), cpdag(network))
        row.cpdag.shd <- c(row.cpdag.shd, cpdag.shd)
        moral.shd <- bnlearn::shd(moral(learned.dag), moral(network))
        row.moral.shd <- c(row.moral.shd, moral.shd)
        bic <- ebic(learned.dag, sim.data.df, gamma=gamma)
        bic.row <- c(bic.row, bic)
        performance.metrics <- calculate.performance.statistics(
          cpdag(learned.dag), cpdag(network)
        )
        fdr.row <- c(fdr.row, performance.metrics$fdr)
        sensitivity.row <- c(sensitivity.row, performance.metrics$sensitivity)        
        
        cpdag.shd.array <- rbind(cpdag.shd.array, row.cpdag.shd)
        moral.shd.array <- rbind(moral.shd.array, row.moral.shd)
        bic.array <- rbind(bic.array, bic.row)
        fdr.array <- rbind(fdr.array, fdr.row)
        sensitivity.array <- rbind(sensitivity.array, sensitivity.row)
      }
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
  cpdag.shd.array <- array(numeric(), c(0, 5))
  # Moralized Structural Hamming Distance
  moral.shd.array <- array(numeric(), c(0, 5))
  # BIC Criterion
  bic.array <- array(numeric(), c(0, 5))
  # False Discovery Rate
  fdr.array <- array(numeric(), c(0, 5))
  # Sensitivity
  sensitivity.array <- array(numeric(), c(0, 5))
  
  colnames(cpdag.shd.array) <- c(c("N", "D", 'seed', "lamdba1", 'cpdag.shd'))
  colnames(moral.shd.array) <- c(c("N", "D", 'seed', "lamdba1", 'moral.shd'))
  colnames(bic.array) <- c(c("N", "D", 'seed', "lambda1", 'bic'))
  colnames(fdr.array) <- c(c("N", "D", 'seed', "lambda1", "fdr"))
  colnames(sensitivity.array) <- c(c("N", "D", 'seed', "lambda1", "sensitivity"))
  n <- 1
  shift <- 0  
  while (n <= n.runs){
    for (l1 in l1.range){
      set.seed(n + shift)
      N <- sample(n.range[1]:n.range[2], 1)
      set.seed(n + shift)
      sim.data.df <- rbn(network, n=N)
      sim.data <- as.matrix(sim.data.df)    
  
      start.time <- Sys.time()
      row <- c(N, D, n+shift, l1)
      row.cpdag.shd <- row
      row.moral.shd <- row
      bic.row <- row
      fdr.row <- row
      sensitivity.row <- row
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
      
      cpdag.shd <- bnlearn::shd(cpdag(learned.dag), cpdag(network))
      row.cpdag.shd <- c(row.cpdag.shd, cpdag.shd)
      moral.shd <- bnlearn::shd(moral(learned.dag), moral(network))
      row.moral.shd <- c(row.moral.shd, moral.shd)
      bic <- ebic(bn.fit(learned.dag, sim.data.df), sim.data.df, gamma=gamma)
      bic.row <- c(bic.row, bic)
      performance.metrics <- calculate.performance.statistics(
        cpdag(learned.dag), cpdag(network)
      )
      fdr.row <- c(fdr.row, performance.metrics$fdr)
      sensitivity.row <- c(sensitivity.row, performance.metrics$sensitivity)        
      
      cpdag.shd.array <- rbind(cpdag.shd.array, row.cpdag.shd)
      moral.shd.array <- rbind(moral.shd.array, row.moral.shd)
      bic.array <- rbind(bic.array, bic.row)
      fdr.array <- rbind(fdr.array, fdr.row)
      sensitivity.array <- rbind(sensitivity.array, sensitivity.row)
    } # FOR alpha.range
    if (is_dag(amat.dag)){
      n = n + 1
    }
  } # while n < n.runs
  return(list(
    'cpdag'=cpdag.shd.array,
    'moral'=moral.shd.array,
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
  cpdag.shd.array <- array(numeric(), c(0, num.algos + 2))
  # Moralized Structural Hamming Distance
  moral.shd.array <- array(numeric(), c(0, num.algos + 2))
  # BIC Criterion
  bic.array <- array(numeric(), c(0, num.algos + 2))
  # False Discovery Rate
  fdr.array <- array(numeric(), c(0, num.algos + 2))
  # Sensitivity
  sensitivity.array <- array(numeric(), c(0, num.algos + 2))
  
  colnames(cpdag.shd.array) <- c(c("N", "alpha"), names(functions.bn))
  colnames(moral.shd.array) <- c(c("N", "alpha"), names(functions.bn))
  colnames(bic.array) <- c(c("N", "alpha"), names(functions.bn))
  colnames(fdr.array) <- c(c("N", "alpha"), names(functions.bn))
  colnames(sensitivity.array) <- c(c("N", "alpha"), names(functions.bn))
  for (data.size in data.size.range){
    for (run in seq(1, n.run)){
      set.seed(run)
      sim.data <- rbn(bn.fitted, n=data.size)
      for (alpha in alpha.range){
        row <- c(data.size, alpha)
        row.cpdag.shd <- row
        row.moral.shd <- row
        bic.row <- row
        fdr.row <- row
        sensitivity.row <- row
        for (algo in 1:num.algos){
          learned.dag <- functions.bn[[algo]](sim.data.df, alpha=a)  
          cpdag.shd <- shd(cpdag(learned.dag), cpdag(bn.fitted))
          row.cpdag.shd <- c(row.cpdag.shd, cpdag.shd)
          moral.shd <- shd(moral(learned.dag), moral(bn.fitted))
          row.moral.shd <- c(row.moral.shd, moral.shd)
          bic <- ebic(learned.dag, sim.data.df, gamma=gamma)
          bic.row <- c(bic.row, bic)
          performance.metrics <- calculate.performance.statistics(
            cpdag(learned.dag), cpdag(bn.fitted)
          )
          fdr.row <- c(fdr.row, performance.metrics$fdr)
          sensitivity.row <- c(sensitivity.row, performance.metrics$sensitivity)
        } # for algo
        cpdag.shd.array <- rbind(cpdag.shd.array, row.cpdag.shd)
        moral.shd.array <- rbind(moral.shd.array, row.moral.shd)
        bic.array <- rbind(bic.array, bic.row)
        fdr.array <- rbind(fdr.array, fdr.row)
        sensitivity.array <- rbind(sensitivity.array, sensitivity.row)
      } # for alpha 
    } # for run
  } # for data.size
  return(list(
    'cpdag'=cpdag.shd.array,
    'moral'=moral.shd.array,
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
  cpdag.shd.array <- array(numeric(), c(0, num.algos + 2))
  # Moralized Structural Hamming Distance
  moral.shd.array <- array(numeric(), c(0, num.algos + 2))
  # BIC Criterion
  bic.array <- array(numeric(), c(0, num.algos + 2))
  # False Discovery Rate
  fdr.array <- array(numeric(), c(0, num.algos + 2))
  # Sensitivity
  sensitivity.array <- array(numeric(), c(0, num.algos + 2))
  
  colnames(cpdag.shd.array) <- c(c("N", "alpha"), names(functions.bn))
  colnames(moral.shd.array) <- c(c("N", "alpha"), names(functions.bn))
  colnames(bic.array) <- c(c("N", "alpha"), names(functions.bn))
  colnames(fdr.array) <- c(c("N", "alpha"), names(functions.bn))
  colnames(sensitivity.array) <- c(c("N", "alpha"), names(functions.bn))
  for (data.size in data.size.range){
    for (run in seq(1, n.run)){
      set.seed(run)
      sim.data <- rbn(bn.fitted, n=data.size)
      row <- c(data.size)
      row.cpdag.shd <- row
      row.moral.shd <- row
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
        
        cpdag.shd <- shd(cpdag(learned.dag), cpdag(bn.fitted))
        row.cpdag.shd <- c(row.cpdag.shd, cpdag.shd)
        moral.shd <- shd(moral(learned.dag), moral(bn.fitted))
        row.moral.shd <- c(row.moral.shd, moral.shd)
        bic <- ebic(learned.dag, sim.data.df, gamma=gamma)
        bic.row <- c(bic.row, bic)
        performance.metrics <- calculate.performance.statistics(
          cpdag(learned.dag), cpdag(bn.fitted)
        )
        fdr.row <- c(fdr.row, performance.metrics$fdr)
        sensitivity.row <- c(sensitivity.row, performance.metrics$sensitivity)
      } # for algo
      cpdag.shd.array <- rbind(cpdag.shd.array, row.cpdag.shd)
      moral.shd.array <- rbind(moral.shd.array, row.moral.shd)
      bic.array <- rbind(bic.array, bic.row)
      fdr.array <- rbind(fdr.array, fdr.row)
      sensitivity.array <- rbind(sensitivity.array, sensitivity.row)
    } # for run
  } # for data.size
  return(list(
    'cpdag'=cpdag.shd.array,
    'moral'=moral.shd.array,
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