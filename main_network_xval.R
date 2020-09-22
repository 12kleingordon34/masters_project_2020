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
pc.f <- partial(pc.stable, test = 'cor', undirected = FALSE,)
gs.f <- partial(gs, test = 'cor', undirected = FALSE)
iamb.f <- partial(iamb, test = 'cor', undirected = FALSE)
fast.iamb.f <- partial(fast.iamb, test = 'cor', undirected = FALSE)
inter.iamb.f <- partial(inter.iamb, test = 'cor', undirected = FALSE)
iamb.fdr.f <- partial(iamb.fdr, test = 'cor', undirected = FALSE)

# Score based methods
tabu.f <- partial(tabu, score = 'bic-g')
hc.f <- partial(hc, perturb=20, score = 'bic-g')

# Regression
notears.f <- partial(notears_linear, loss_type='l2')

constraints.bn.func <- list(
  "pc"=pc.f,
  # "gs"=gs.f,
  "iamb"=iamb.f,
  "fast.iamb"=fast.iamb.f,
  "inter.iamb"= inter.iamb.f,
  "iamb.fdr"=iamb.fdr.f
)
n.runs <- 35
n.range <- c(25, 300)
network.name <- 'arth'
ecoli <- readRDS('bnlearn_networks/arth150.rds')
alpha.range <- c(0.15, 0.1, 0.05, 0.01, 0.005, 0.001)
# alpha.const.results <- alpha.xval.network(
#   constraints.bn.func,
#   alpha.range,
#   n.runs,
#   n.range,
#   ecoli
# )
# 
# save.template <- "data/results/simulated_networks/%s/alpha_validation/%s_runs_%s.csv"
# write.csv(alpha.const.results$shd, sprintf(save.template, network.name, 'shd_cpdag', n.runs), row.names = FALSE)
# write.csv(alpha.const.results$hd, sprintf(save.template, network.name, 'shd_moral', n.runs), row.names = FALSE)
# write.csv(alpha.const.results$fdr, sprintf(save.template, network.name, 'fdr', n.runs), row.names = FALSE)
# write.csv(alpha.const.results$sensitivity, sprintf(save.template, network.name, 'sensitivity', n.runs), row.names = FALSE)
# 
#tabu.range <- c(0.25, 0.5, 1, 2.5, 5, 10)
#tabu.const.results <- tabu.xval.network(
#  tabu.range,
#  n.runs,
#  n.range,
#  ecoli
#)
#save.template <- "data/results/simulated_networks/%s/tabu_validation/%s_runs_%s.csv"
#write.csv(tabu.const.results$shd, sprintf(save.template, network.name, 'shd_cpdag', n.runs), row.names = FALSE)
#write.csv(tabu.const.results$hd, sprintf(save.template, network.name, 'shd_moral', n.runs), row.names = FALSE)
#write.csv(tabu.const.results$fdr, sprintf(save.template, network.name, 'fdr', n.runs), row.names = FALSE)
#write.csv(tabu.const.results$sensitivity, sprintf(save.template, network.name, 'sensitivity', n.runs), row.names = FALSE)

#restart.range <- c(1, 5, 10, 25, 50)
#perturb.range <- c(1, 5, 10, 25, 50)
#hc.const.results <- hc.xval.network(
#  restart.range,
#  perturb.range,
#  n.runs,
#  n.range,
#  ecoli
#)
#save.template <- "data/results/simulated_networks/%s/hc_validation/%s_runs_%s.csv"
#write.csv(hc.const.results$shd, sprintf(save.template, network.name, 'shd_cpdag', n.runs), row.names = FALSE)
#write.csv(hc.const.results$hd, sprintf(save.template, network.name, 'shd_moral', n.runs), row.names = FALSE)
#write.csv(hc.const.results$bic, sprintf(save.template, network.name, 'bic', n.runs), row.names = FALSE)
#write.csv(hc.const.results$fdr, sprintf(save.template, network.name, 'fdr', n.runs), row.names = FALSE)
#write.csv(hc.const.results$sensitivity, sprintf(save.template, network.name, 'sensitivity', n.runs), row.names = FALSE)

# l1.range <- c(0.001, 0.01, 0.1, 1)
# notears.const.results <- notears.xval.network(
#   l1.range,
#   n.runs,
#   n.range,
#   ecoli
# )
# save.template <- "data/results/simulated_networks/%s/notears_validation/%s_runs_%s_v4.csv"
# write.csv(notears.const.results$shd, sprintf(save.template, network.name, 'shd', n.runs), row.names = FALSE)
# write.csv(notears.const.results$hd, sprintf(save.template, network.name, 'shd_moral', n.runs), row.names = FALSE)
# write.csv(notears.const.results$bic, sprintf(save.template, network.name, 'bic', n.runs), row.names = FALSE)
# write.csv(notears.const.results$fdr, sprintf(save.template, network.name, 'fdr', n.runs), row.names = FALSE)
# write.csv(notears.const.results$sensitivity, sprintf(save.template, network.name, 'sensitivity', n.runs), row.names = FALSE)

# # Hybrid Approaches
# h2pc.f <- partial(h2pc, maximize.args=list(perturb=50, restarts=50))
# mmhc.f <- partial(mmhc, maximize.args=list(perturb=50, restarts=50))
# hiton.f <- partial(rsmax2, restrict='si.hiton.pc', maximize='hc', maximize.args=list(perturb=50, restarts=50))
# iamb.tabu.f <- partial(
#  rsmax2,
#  restrict='inter.iamb',
#  maximize='tabu'
# )
# hybrid.algos <- list(
#  'hybrid.h2pc.hc'=h2pc.f,
#  'hybrid.mmhc.hc'=mmhc.f,
#  'hybrid.hiton.hc'=hiton.f,
#  'hybrid.iamb.tabu'=iamb.tabu.f
# )
# hybrid.const.results = alpha.xval.network(
#  hybrid.algos,
#  alpha.range,
#  n.runs,
#  n.range,
#  ecoli
# )
# save.template <- "data/results/simulated_networks/%s/hybrid/%s_runs_%s.csv"
# write.csv(hybrid.const.results$shd, sprintf(save.template, network.name, 'shd_cpdag', n.runs), row.names = FALSE)
# write.csv(hybrid.const.results$hd, sprintf(save.template, network.name, 'shd_moral', n.runs), row.names = FALSE)
# write.csv(hybrid.const.results$fdr, sprintf(save.template, network.name, 'fdr', n.runs), row.names = FALSE)
# write.csv(hybrid.const.results$sensitivity, sprintf(save.template, network.name, 'sensitivity', n.runs), row.names = FALSE)


# Compare impact of varying network size/data size
pc.f <- partial(pc.stable, test = 'cor', undirected = FALSE, alpha=0.05)
gs.f <- partial(gs, test = 'cor', undirected = FALSE, alpha=0.05)
iamb.f <- partial(iamb, test = 'cor', undirected = FALSE, alpha=0.05)
fast.iamb.f <- partial(fast.iamb, test = 'cor', undirected = FALSE, alpha=0.05)
inter.iamb.f <- partial(inter.iamb, test = 'cor', undirected = FALSE, alpha=0.05)
iamb.fdr.f <- partial(iamb.fdr, test = 'cor', undirected = FALSE, alpha=0.05)
h2pc.f <- partial(h2pc, maximize.args=list(perturb=50, restart=50), restrict.args=list(alpha=0.2))
mmhc.f <- partial(mmhc, maximize.args=list(perturb=50, restart=50), restrict.args=list(alpha=0.2))
hiton.hc.f <- partial(
 rsmax2,
 restrict='si.hiton.pc',
 maximize='hc',
 maximize.args=list(perturb=20, restart=20),
 restrict.args=list(alpha=0.2)
)
iamb.hc.f <- partial(
 rsmax2,
 restrict='inter.iamb',
 maximize='hc',
 restrict.args=list(alpha=0.2),
 maximize.args=list(score='bge', perturb=50, restart=50)
)

# Regression
notears.f <- partial(notears_linear, loss_type='l2', lambda1=0.1)
hybrid.notears.f <- partial(hybrid.notears, params=list(l1=0.1, alpha=0.1, test='cor', B=100))

hc.f <- partial(hc, score='bge', perturb=50, restart=50)

learning.functions <- list(
 "pc"=pc.f,
 "gs"=gs.f,
 "iamb"=iamb.f,
 "fast.iamb"=fast.iamb.f,
 "inter.iamb"= inter.iamb.f,
 "iamb.fdr"=inter.iamb.f,
 'hc'=hc.f,
 'h2pc'=h2pc.f,
 'mmhc'=mmhc.f,
 'hiton.hc'=hiton.hc.f,
 # 'iamb.hc'=iamb.hc.f,
 'notears'=notears.f,
 'hybrid.notears.f'=hybrid.notears.f
)

n.range <- c(50, 100, 250, 500)
args = commandArgs(trailingOnly=TRUE)
n.range <- as.integer(args[1])
print(n.range)
n.runs <- 35
size.data <-data.size.study.network(learning.functions, n.range, n.runs, ecoli)

save.template <- "data/results/simulated_networks/%s/confusion_matrices/%s_%s.csv"
write.csv(size.data$shd, sprintf(save.template, network.name, 'shd', n.range), row.names = FALSE)
write.csv(size.data$hd, sprintf(save.template, network.name, 'hd', n.range), row.names = FALSE)
write.csv(size.data$fdr, sprintf(save.template, network.name, 'fdr', n.range), row.names = FALSE)
write.csv(size.data$sensitivity, sprintf(save.template, network.name, 'sensitivity', n.range), row.names = FALSE)

# # Confusion
# confusion.data = confusion.shd.simulation(learning.functions, 10000, 15, ecoli)
# save.template <- "data/results/simulated_networks/%s/confusion_matrices/%s_500.csv"
# write.csv(
#   apply(confusion.data$shd.relative.cpdag, c(1,2), mean),
#   sprintf(save.template, network.name, 'shd_rel'),
#   row.names = TRUE
# )
# write.csv(
#   apply(confusion.data$shd.relative.moral, c(1,2), mean),
#   sprintf(save.template, network.name, 'moral_shd_rel'),
#   row.names = TRUE
# )
