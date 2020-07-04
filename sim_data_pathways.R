library(bnlearn)

generate.network.data <- function(N, network.path, seed=0){
  #' Generates simulated data from a well documented disease
  #' pathway. The network file must be in RDS format. The 
  #' networks can be downloaded from: 
  #' https://www.bnlearn.com/bnrepository/
  
  network.bn <- readRDS(network.path)
  set.seed(seed)
  sim.data <- rbn(network.bn, N)
  return (list('data'=sim.data, 'network'=network.bn))
}


fit.network.to.data <- function(algorithm, N, network.path, seed=0){
  #' 
  sim.network <- generate.network.data(N, network.path, seed)
  network.data <- sim.network$data; bn.network <- sim.network$network.bn
  fitted.bn <- algorithm(network.data)
  
  shd <- shd(cpdag(fitted.bn), cpdag(sim.network$network))
  return (list('shd'=shd, 'fitted.bn'=fitted.bn))
}

run.shd.simulation <- function(algorithms, data.size, run.number, network.path){
  #' 
  shd.array.true <- array(
    numeric(),
    c(run.number, length(algorithms)),
  )
  colnames(shd.array.true) <- names(algorithms)
  
  shd.array.rel <- array(
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
    for (algo in seq(1, length(algorithms), 1)){
      fitted.network <- fit.network.to.data(
        algorithms[[algo]], data.size, network.path, seed=r
      )
      shd.array.true[r, algo] <- fitted.network$shd
      algo.graphs[[algo]] <- fitted.network$fitted.bn
      print(sprintf(
        '%s - %s - Run: %s - SHD: %s',
        Sys.time(), names(algorithms)[algo], r, fitted.network$shd
      ))    
    }

    for (i in seq(1, length(algo.graphs), 1)){
      for (j in seq(1, length(algo.graphs), 1)){
        shd.array.rel[i, j, r] = shd(cpdag(algo.graphs[[i]]), cpdag(algo.graphs[[j]]))
      }
    }
  }
  return(list(
    'shd.true'=data.frame(shd.array.true),
    'shd.relative'=shd.array.rel
  ))
}



# z = melt(z, value.name='shd')
# ggplot(z, aes(x=variable, y=shd)) + geom_violin()