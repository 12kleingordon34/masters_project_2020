library(bnlearn)
library(equSA)
library(Rgraphviz)
library(purrr)

source("sim_data_generation.R")

## Define functions
# Constraint Methods
pc.f <- partial(pc.stable, test = 'mi-g-sh', undirected = FALSE)
gs.f <- partial(gs, test = 'mi-g-sh', undirected = FALSE)
iamb.f <- partial(iamb, test = 'mi-g-sh', undirected = FALSE)
fast.iamb.f <- partial(fast.iamb, test = 'mi-g-sh', undirected = FALSE)
inter.iamb.f <- partial(inter.iamb, test = 'mi-g-sh', undirected = FALSE)
iamb.fdr.f <- partial(iamb.fdr, test = 'mi-g-sh', undirected = FALSE)
mmpc.f <- partial(mmpc, test = 'mi-g-sh', undirected = FALSE)
si.hiton.pc.f <- partial(si.hiton.pc, test = 'mi-g-sh', undirected = FALSE)
hpc.f <- partial(hpc, test = 'mi-g-sh', undirected = FALSE, debug=FALSE)

# Score based methods
tabu.f <- partial(tabu, score = 'bge')
hc.f <- partial(hc, score = 'bge')

# Hybrid methods


constraint.bn.func <- list(
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

# Define hyperparameter ranges
alpha.range <- c(0.1, 0.05, 0.01); tabu.range<- c(0.1, 0.5, 1, 1.5, 3, 7); 

# Define graph parameter ranges
D.range <- seq(15, 155, 10)
sparsity.range <- seq(0.1, 0.3, 0.1)
N <- 100 # Num of simulated observations
S <- 5 # Number of seeded runs

for (sparsity in sparsity.range){
  start.time <- Sys.time()
  sim.shd.results <- run.simulation.test(N, S, D.range, sparsity, constraint.bn.func, alpha.range, tabu.range)
  end.time <- Sys.time()
  print(sprintf("Start: %s. Finish: %s. Duration: %s", start.time, end.time, end.time-start.time))
  
  mean.shd <- sim.shd.results$mean
  std.shd <- sim.shd.results$std
  best.hyperparams <- sim.shd.results$best.params
  run.times <- sim.shd.results$times
  
  save.template <- "data/results/simulated_gaussian/constraint_%s_N_%s_shd_sparsity_%s.csv"
  write.csv(mean.shd, sprintf(save.template, "mean",  N, sparsity))
  write.csv(std.shd, sprintf(save.template, "std", N, sparsity))
  write.csv(best.hyperparams, sprintf(save.template, "hyperparams", N, sparsity))
  write.csv(run.times, sprintf(save.template,"runtimes", N,  sparsity))
}
