library(bnlearn)


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

pc.f <- partial(pc.stable, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
gs.f <- partial(gs, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
iamb.f <- partial(iamb, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
fast.iamb.f <- partial(fast.iamb, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
inter.iamb.f <- partial(inter.iamb, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
iamb.fdr.f <- partial(iamb.fdr, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
mmpc.f <- partial(mmpc, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
si.hiton.pc.f <- partial(si.hiton.pc, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
hpc.f <- partial(hpc, test = 'mi-g-sh', undirected = FALSE, debug=FALSE, alpha=alpha)

algorithms <- list(
  "pc"=pc.f,
  "gs"=gs.f,
  "iamb"=iamb.f,
  "fast.iamb"=fast.iamb.f,
  "inter.iamb"= inter.iamb.f,
  "iamb.fdr"=inter.iamb.f,
  "mmpc"=mmpc.f,
  "si.hiton.pc"=si.hiton.pc.f,
  "hpc"=hpc.f
)
ecoli <- readRDS('bnlearn_networks/ecoli70.rds')
N <- 15000
n.runs <- 10
ecoli.run = confusion.shd.simulation(algorithms, N, n.runs, ecoli)
save.template <- "data/results/simulated_networks/ecoli/%s/%_N_%s_run_%s.csv"
write.csv(
  apply(ecoli.run$shd.relative.cpdag, c(1,2), mean), 
  sprintf(save.template, 'cpdag', "mean", N, n.runs), 
  row.names = FALSE
)
write.csv(
  apply(ecoli.run$shd.relative.moral, c(1,2), mean), 
  sprintf(save.template, 'moral', "mean", N, n.runs), 
  row.names = FALSE
)

# # z = melt(z, value.name='shd')
# ggplot(z, aes(x=variable, y=shd)) + geom_violin()