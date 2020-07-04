library(bnlearn)
library(equSA)
library(ggplot2)
library(reshape2)
library(Rgraphviz)
library(purrr)
library(plyr)

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

# # Hybrid methods
# mmpc.f <- partial(rsmax2, restrict='mmpc', maximise='hc', restrict.args=list(alpha=0.01), 
#                   maximize.args=list(restart=10, perturb=5, max.iter=100000))

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
alpha.range <- c(0.1, 0.05, 0.01, 0.001); tabu.range<- c(0.1, 0.5, 1, 1.5, 3, 7); 

# Define graph parameter ranges
D.range <- seq(15, 120, 15)
sparsity.range <- c(0.1)
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
  write.csv(mean.shd, sprintf(save.template, "mean",  N, sparsity), row.names = FALSE)
  write.csv(std.shd, sprintf(save.template, "std", N, sparsity), row.names = FALSE)
  write.csv(best.hyperparams, sprintf(save.template, "hyperparams", N, sparsity), row.names = FALSE)
  write.csv(run.times, sprintf(save.template,"runtimes", N,  sparsity), row.names = FALSE)
}

# Plot graphs
for (sparsity in sparsity.range) {
  mean <- read.table(
    sprintf('data/results/simulated_gaussian/constraint_mean_N_100_shd_sparsity_%s.csv', sparsity),
    header=TRUE, row.names=NULL, sep=','
  )
  mean$n.over.d <- mean$N / mean$D
  mean <- subset(mean, select=c(-1,-2,-3))
  std <- read.table(
    sprintf('data/results/simulated_gaussian/constraint_std_N_100_shd_sparsity_%s.csv', sparsity),
    header=TRUE, row.names=NULL, sep=','
  )
  std$n.over.d <- std$N / std$D
  std <- subset(std, select=c(-1,-2,-3))
  times <- read.table(
    sprintf('data/results/simulated_gaussian/constraint_runtimes_N_100_shd_sparsity_%s.csv', sparsity),
    header=TRUE, row.names=NULL, sep=','
  )
  times$n.over.d <- 100 / times$D
  times <- subset(times, select=c(-1,-2))

  mean.melted <- melt(mean, id="n.over.d")
  std.melted <- melt(std, id="n.over.d")
  times.melted <- melt(times, id='n.over.d')
  mean.melted <- rename(mean.melted, c('value'='shd'))
  std.melted <- rename(std.melted, c('value'='error'))
  times.melted <- rename(times.melted, c('value'='time'))

  total <- merge(mean.melted, std.melted, by=c('n.over.d', 'variable'))

  plot.shd <- ggplot(subset(total, n.over.d < 3), aes(x=n.over.d, y=shd, group=variable, ymin=shd-error, ymax=shd+error, colour=variable)) +
    geom_line() + geom_point() + scale_y_log10() +
    geom_errorbar(width=.2, position=position_dodge(0.05))
  plot.shd + labs(title=sprintf('Structural Hamming distance vs N/D - Sparsity=%s', sparsity), x='N/D', y='SHD')
  
  plot.times <- ggplot(times.melted, aes(x=n.over.d, y=time, group=variable, colour=variable)) + 
    geom_line() + geom_point() + scale_y_log10()
  plot.times + labs(title=sprintf('Algorithm run times - Sparsity=%s', sparsity), x='N/D', y='Time/s')
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
