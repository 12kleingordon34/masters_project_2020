library(reticulate)
# import('scipy')
use_condaenv("r-reticulate")
source_python('./notears/notears/linear.py')
source_python('./notears/notears/utils.py')

library(bnlearn)
library(equSA)
# library(ggplot2)
library(reshape2)
library(Rgraphviz)
library(purrr)
library(plyr)

source("./sim_data_pathways.R")

## Define functions
# Constraint Methods
pc.f <- partial(pc.stable, test = 'mi-g-sh', undirected = FALSE,)
gs.f <- partial(gs, test = 'mi-g-sh', undirected = FALSE)
iamb.f <- partial(iamb, test = 'mi-g-sh', undirected = FALSE)
fast.iamb.f <- partial(fast.iamb, test = 'mi-g-sh', undirected = FALSE)
inter.iamb.f <- partial(inter.iamb, test = 'mi-g-sh', undirected = FALSE)
iamb.fdr.f <- partial(iamb.fdr, test = 'mi-g-sh', undirected = FALSE)

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
n.runs <- 50
n.range <- c(100, 300)
network.name <- 'ecoli'
ecoli <- readRDS('bnlearn_networks/ecoli70.rds')
alpha.range <- c(0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001)
alpha.const.results <- alpha.xval.network(
  constraints.bn.func,
  alpha.range,
  n.runs,
  n.range,
  ecoli
)

save.template <- "data/results/simulated_networks/%s/alpha_validation/%s_runs_%s.csv"
write.csv(alpha.const.results$cpdag, sprintf(save.template, network.name, 'shd_cpdag', n.runs), row.names = FALSE)
write.csv(alpha.const.results$moral, sprintf(save.template, network.name, 'shd_moral', n.runs), row.names = FALSE)
write.csv(alpha.const.results$bic, sprintf(save.template, network.name, 'bic', n.runs), row.names = FALSE)
write.csv(alpha.const.results$fdr, sprintf(save.template, network.name, 'fdr', n.runs), row.names = FALSE)
write.csv(alpha.const.results$sensitivity, sprintf(save.template, network.name, 'sensitivity', n.runs), row.names = FALSE)

tabu.range <- c(0.1, 0.25, 0.5, 1, 5)
tabu.const.results <- tabu.xval.network(
  tabu.range,
  n.runs,
  n.range,
  ecoli
)
save.template <- "data/results/simulated_networks/%s/tabu_validation/%s_runs_%s.csv"
write.csv(tabu.const.results$cpdag, sprintf(save.template, network.name, 'shd_cpdag', n.runs), row.names = FALSE)
write.csv(tabu.const.results$moral, sprintf(save.template, network.name, 'shd_moral', n.runs), row.names = FALSE)
write.csv(tabu.const.results$bic, sprintf(save.template, network.name, 'bic', n.runs), row.names = FALSE)
write.csv(tabu.const.results$fdr, sprintf(save.template, network.name, 'fdr', n.runs), row.names = FALSE)
write.csv(tabu.const.results$sensitivity, sprintf(save.template, network.name, 'sensitivity', n.runs), row.names = FALSE)

restart.range <- c(0, 1, 5, 10, 25, 50)
hc.const.results <- hc.xval.network(
  restart.range,
  n.runs,
  n.range,
  ecoli
)
save.template <- "data/results/simulated_networks/%s/hc_validation/%s_runs_%s.csv"
write.csv(tabu.const.results$cpdag, sprintf(save.template, network.name, 'shd_cpdag', n.runs), row.names = FALSE)
write.csv(tabu.const.results$moral, sprintf(save.template, network.name, 'shd_moral', n.runs), row.names = FALSE)
write.csv(tabu.const.results$bic, sprintf(save.template, network.name, 'bic', n.runs), row.names = FALSE)
write.csv(tabu.const.results$fdr, sprintf(save.template, network.name, 'fdr', n.runs), row.names = FALSE)
write.csv(tabu.const.results$sensitivity, sprintf(save.template, network.name, 'sensitivity', n.runs), row.names = FALSE)

l1.range <- c(0.001, 0.01, 0.1, 1, 10)
notears.const.results <- notears.xval.network(
  l1.range,
  n.runs,
  n.range,
  ecoli
)
save.template <- "data/results/simulated_networks/%s/notears_validation/%s_runs_%s.csv"
write.csv(notears.const.results$cpdag, sprintf(save.template, network.name, 'shd_cpdag', n.runs), row.names = FALSE)
write.csv(notears.const.results$moral, sprintf(save.template, network.name, 'shd_moral', n.runs), row.names = FALSE)
write.csv(notears.const.results$bic, sprintf(save.template, network.name, 'bic', n.runs), row.names = FALSE)
write.csv(notears.const.results$fdr, sprintf(save.template, network.name, 'fdr', n.runs), row.names = FALSE)
write.csv(notears.const.results$sensitivity, sprintf(save.template, network.name, 'sensitivity', n.runs), row.names = FALSE)
