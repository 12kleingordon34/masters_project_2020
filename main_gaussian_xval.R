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
hc.f <- partial(hc, score = 'bic-g')

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
n.runs <- 75
d.range <- c(20, 120)
n.range <- c(40, 200)
alpha.range <- c(0.2, 0.1, 0.05, 0.01, 0.005, 0.001)
# alpha.const.results <- alpha.constraint.xval.ranDAG(
#   constraints.bn.func,
#   alpha.range,
#   n.runs,
#   d.range = d.range,
#   n.range
# )
# 
# save.template <- "data/results/simulated_gaussian/alpha_validation/%s_runs_%s.csv"
# write.csv(alpha.const.results$shd, sprintf(save.template, 'shd', n.runs), row.names = FALSE)
# write.csv(alpha.const.results$hd, sprintf(save.template, 'hd', n.runs), row.names = FALSE)
# write.csv(alpha.const.results$bic, sprintf(save.template, 'bic', n.runs), row.names = FALSE)
# write.csv(alpha.const.results$fdr, sprintf(save.template, 'fdr', n.runs), row.names = FALSE)
# write.csv(alpha.const.results$sensitivity, sprintf(save.template, 'sensitivity', n.runs), row.names = FALSE)
# 
# tabu.range <- c(0.25, 0.5, 1, 2.5, 5)
# tabu.const.results <- tabu.xval.ranDAG(
#   tabu.range, n.runs, d.range, n.range
# )
# save.template <- "data/results/simulated_gaussian/tabu_validation/%s_runs_%s.csv"
# write.csv(tabu.const.results$shd, sprintf(save.template, 'shd', n.runs), row.names = FALSE)
# write.csv(tabu.const.results$hd, sprintf(save.template, 'hd', n.runs), row.names = FALSE)
# write.csv(tabu.const.results$bic, sprintf(save.template, 'bic', n.runs), row.names = FALSE)
# write.csv(tabu.const.results$fdr, sprintf(save.template, 'fdr', n.runs), row.names = FALSE)
# write.csv(tabu.const.results$sensitivity, sprintf(save.template, 'sensitivity', n.runs), row.names = FALSE)
# 
# HC Restart
restart.range <- c(1, 10, 25, 50)
perturb.range <- c(1, 10, 25, 50)
hc.const.results <- hc.xval.ranDAG(restart.range, perturb.range, n.runs, d.range, n.range)
save.template <- "data/results/simulated_gaussian/hc_validation/%s_runs_%s.csv"
write.csv(hc.const.results$shd, sprintf(save.template, 'shd', n.runs), row.names = FALSE)
write.csv(hc.const.results$hd, sprintf(save.template, 'hd', n.runs), row.names = FALSE)
write.csv(hc.const.results$bic, sprintf(save.template, 'bic', n.runs), row.names = FALSE)
write.csv(hc.const.results$fdr, sprintf(save.template, 'fdr', n.runs), row.names = FALSE)
write.csv(hc.const.results$sensitivity, sprintf(save.template, 'sensitivity', n.runs), row.names = FALSE)

# # NOTEARS
# l1.range <- c(0.0001, 0.001, 0.01, 0.1, 1)
# notears.const.results <- notears.xval.ranDAG(l1.range, n.runs, d.range, n.range)
# save.template <- "data/results/simulated_gaussian/notears_validation/%s_runs_%s.csv"
# write.csv(notears.const.results$shd, sprintf(save.template, 'shd', n.runs), row.names = FALSE)
# write.csv(notears.const.results$hd, sprintf(save.template, 'hd', n.runs), row.names = FALSE)
# write.csv(notears.const.results$bic, sprintf(save.template, 'bic', n.runs), row.names = FALSE)
# write.csv(notears.const.results$fdr, sprintf(save.template, 'fdr', n.runs), row.names = FALSE)
# write.csv(notears.const.results$sensitivity, sprintf(save.template, 'sensitivity', n.runs), row.names = FALSE)


# # Hybrid Approaches
# h2pc.f <- partial(h2pc, maximize.args=list(perturb=20, restart=20))
# mmhc.f <- partial(mmhc, maximize.args=list(perturb=20, restart=20))
# hiton.f <- partial(rsmax2, restrict='si.hiton.pc', maximize='hc')
# iamb.tabu.f <- partial(
#   rsmax2,
#   restrict='inter.iamb',
#   maximize='tabu',
# )
# hybrid.algos <- list(
#   'hybrid.h2pc.hc'=h2pc.f,
#   'hybrid.mmhc.hc'=mmhc.f,
#   'hybrid.hiton.hc'=hiton.f,
#   'hybrid.iamb.tabu'=iamb.tabu.f
# )
# hybrid.const.results = alpha.constraint.xval.ranDAG(
#   hybrid.algos,
#   alpha.range,
#   n.runs,
#   d.range = d.range,
#   n.range
# )
# 
# save.template <- "data/results/simulated_gaussian/alpha_validation/%s_runs_%s.csv"
# write.csv(alpha.const.results$shd, sprintf(save.template, 'shd', n.runs), row.names = FALSE)
# write.csv(alpha.const.results$hd, sprintf(save.template, 'hd', n.runs), row.names = FALSE)
# write.csv(alpha.const.results$bic, sprintf(save.template, 'bic', n.runs), row.names = FALSE)
# write.csv(alpha.const.results$fdr, sprintf(save.template, 'fdr', n.runs), row.names = FALSE)
# write.csv(alpha.const.results$sensitivity, sprintf(save.template, 'sensitivity', n.runs), row.names = FALSE)

# # Compare impact of varying network size/data size
# pc.f <- partial(pc.stable, test = 'mi-g-sh', undirected = FALSE, alpha=0.05)
# gs.f <- partial(gs, test = 'mi-g-sh', undirected = FALSE, alpha=0.05)
# iamb.f <- partial(iamb, test = 'mi-g-sh', undirected = FALSE, alpha=0.05)
# fast.iamb.f <- partial(fast.iamb, test = 'mi-g-sh', undirected = FALSE, alpha=0.05)
# inter.iamb.f <- partial(inter.iamb, test = 'mi-g-sh', undirected = FALSE, alpha=0.05)
# iamb.fdr.f <- partial(iamb.fdr, test = 'mi-g-sh', undirected = FALSE, alpha=0.05)
# h2pc.f <- partial(h2pc, maximize.args=list(perturb=20, restart=20), restrict.args=list(alpha=0.05))
# mmhc.f <- partial(mmhc, maximize.args=list(perturb=20, restart=20), restrict.args=list(alpha=0.05))
# hiton.hc.f <- partial(rsmax2, restrict='si.hiton.pc', maximize='hc', restrict.args=list(alpha=0.05))
# iamb.tabu.f <- partial(
#   rsmax2,
#   restrict='inter.iamb',
#   maximize='tabu',
# )
# 
# # Regression
# notears.f <- partial(notears_linear, loss_type='l2', lambda1=0.01)
# 
# 
# learning.functions <- list(
#   "pc"=pc.f,
#   "gs"=gs.f,
#   "iamb"=iamb.f,
#   "fast.iamb"=fast.iamb.f,
#   "inter.iamb"= inter.iamb.f,
#   "iamb.fdr"=inter.iamb.f,
#   'h2pc'=h2pc.f,
#   'mmhc'=mmhc.f,
#   'hiton.hc'=hiton.hc.f,
#   'iamb.tabu'=iamb.tabu.f,
#   'notears'=notears.f
# )
# 
# n.range <- c(25, 50, 75, 100, 125, 150)
# d.range <- c(50, 100, 150)
# n.runs <- 20
# size.data <-data.size.study.gaussian(learning.functions, n.range, d.range, n.runs)
# 
# save.template <- "data/results/simulated_gaussian/data_size/%s.csv"
# write.csv(size.data$shd, sprintf(save.template, 'shd'), row.names = FALSE)
# write.csv(size.data$hd, sprintf(save.template, 'hd'), row.names = FALSE)
# write.csv(size.data$fdr, sprintf(save.template, 'fdr'), row.names = FALSE)
# write.csv(size.data$sensitivity, sprintf(save.template, 'sensitivity'), row.names = FALSE)
