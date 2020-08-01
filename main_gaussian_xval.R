library(reticulate)
# import('scipy')
use_condaenv("r-reticulate")
source_python('./notears/notears/linear.py')

library(bnlearn)
library(equSA)
# library(ggplot2)
library(reshape2)
library(Rgraphviz)
library(purrr)
library(plyr)

source("./sim_gaussian_xval.R")

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
hpc.f <- partial(hpc, test = 'mi-g-sh', undirected = FALSE)

# Score based methods
tabu.f <- partial(tabu, score = 'bic-g')
hc.f <- partial(hc, perturb=20, score = 'bic-g')

# Regression
notears.f <- partial(notears_linear, loss_type='l2')

constraints.bn.func <- list(
  "pc"=pc.f,
  "gs"=gs.f,
  "iamb"=iamb.f,
  "fast.iamb"=fast.iamb.f,
  "inter.iamb"= inter.iamb.f,
  "iamb.fdr"=inter.iamb.f
)
n.runs <- 100
d.range <- c(20, 120)
max.n <- 300
alpha.range <- c(0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001)
alpha.const.results <- alpha.constraint.xval.ranDAG(
  constraints.bn.func,
  alpha.range,
  n.runs,
  d.range = d.range,
  max.n
)

save.template <- "data/results/simulated_gaussian/alpha_validation/%s_runs_%s.csv"
write.csv(alpha.const.results$shd, sprintf(save.template, 'shd', n.runs), row.names = FALSE)
write.csv(alpha.const.results$hd, sprintf(save.template, 'hd', n.runs), row.names = FALSE)
write.csv(alpha.const.results$bic, sprintf(save.template, 'bic', n.runs), row.names = FALSE)
write.csv(alpha.const.results$fdr, sprintf(save.template, 'fdr', n.runs), row.names = FALSE)
write.csv(alpha.const.results$sensitivity, sprintf(save.template, 'sensitivity', n.runs), row.names = FALSE)

tabu.range <- c(0.1, 0.25, 0.5, 1, 5)
tabu.const.results <- tabu.xval.ranDAG(
  tabu.range, n.runs, d.range, max.n
)
save.template <- "data/results/simulated_gaussian/tabu_validation/%s_runs_%s.csv"
write.csv(tabu.const.results$shd, sprintf(save.template, 'shd', n.runs), row.names = FALSE)
write.csv(tabu.const.results$hd, sprintf(save.template, 'hd', n.runs), row.names = FALSE)
write.csv(tabu.const.results$bic, sprintf(save.template, 'bic', n.runs), row.names = FALSE)
write.csv(tabu.const.results$fdr, sprintf(save.template, 'fdr', n.runs), row.names = FALSE)
write.csv(tabu.const.results$sensitivity, sprintf(save.template, 'sensitivity', n.runs), row.names = FALSE)

# HC Restart
restart.range <- c(0, 1, 5, 10, 25, 50)
hc.const.results <- hc.xval.ranDAG(restart.range, NULL, n.runs, d.range, max.n)
save.template <- "data/results/simulated_gaussian/hc_validation/%s_runs_%s.csv"
write.csv(hc.const.results$shd, sprintf(save.template, 'shd', n.runs), row.names = FALSE)
write.csv(hc.const.results$hd, sprintf(save.template, 'hd', n.runs), row.names = FALSE)
write.csv(hc.const.results$bic, sprintf(save.template, 'bic', n.runs), row.names = FALSE)
write.csv(hc.const.results$fdr, sprintf(save.template, 'fdr', n.runs), row.names = FALSE)
write.csv(hc.const.results$sensitivity, sprintf(save.template, 'sensitivity', n.runs), row.names = FALSE)

l1.range <- c(0.0001, 0.001, 0.01, 0.1, 1)
notears.const.results <- notears.xval.ranDAG(l1.range, n.runs, d.range, max.n)
save.template <- "data/results/simulated_gaussian/notears_validation/%s_runs_%s.csv"
write.csv(notears.const.results$shd, sprintf(save.template, 'shd', n.runs), row.names = FALSE)
write.csv(notears.const.results$hd, sprintf(save.template, 'hd', n.runs), row.names = FALSE)
write.csv(notears.const.results$bic, sprintf(save.template, 'bic', n.runs), row.names = FALSE)
write.csv(notears.const.results$fdr, sprintf(save.template, 'fdr', n.runs), row.names = FALSE)
write.csv(notears.const.results$sensitivity, sprintf(save.template, 'sensitivity', n.runs), row.names = FALSE)
