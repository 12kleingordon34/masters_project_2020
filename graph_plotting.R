library(bnlearn)
library(ggplot2)
library(purrr)
library(plyr)
library(reshape2)

ALPHA.CPDAG.PATH <- './data/results/simulated_gaussian/alpha_validation/shd_cpdag_runs_50.csv'
ALPHA.MORAL.PATH <- './data/results/simulated_gaussian/alpha_validation/shd_moral_runs_50.csv'
ALPHA.FDR.PATH <- './data/results/simulated_gaussian/alpha_validation/fdr_runs_50.csv'
ALPHA.SENSITIVITY.PATH <- './data/results/simulated_gaussian/alpha_validation/sensitivity_runs_50.csv'
ALPHA.PATHS <- list(
  'cpdag'=ALPHA.CPDAG.PATH, 
  'moral'=ALPHA.MORAL.PATH,
  'fdr'=ALPHA.FDR.PATH,
  'sensitivity'=ALPHA.SENSITIVITY.PATH
)

# Plotting histograms
melted.alpha.full <- NULL
mu.full <- NULL
for (i in 1:length(ALPHA.PATHS)){
  path <- ALPHA.PATHS[[i]]
  name <- names(ALPHA.PATHS)[i]
  alpha.sim.gaussian <- read.csv(path)
  melted.alpha <- melt(
    alpha.sim.gaussian,
    id.vars='alpha',
    measure.vars=c('pc', 'gs', 'iamb', 'fast.iamb', 'inter.iamb', 'iamb.fdr')
  )
  melted.alpha['alpha'] = lapply(melted.alpha['alpha'], as.character)
  melted.alpha$data.type <- name
  mu <- ddply(melted.alpha, c('alpha', 'variable'), summarise, value=mean(value))
  mu$data.type <- name
  
  if (is.null(melted.alpha.full)){
    melted.alpha.full = melted.alpha
  } else {
    melted.alpha.full = rbind(melted.alpha.full, melted.alpha)
  }
  if (is.null(mu.full)){
    mu.full <- mu
  } else {
    mu.full <- rbind(mu.full, mu)
  }
}
melted.alpha.full$data.type <- factor(
  melted.alpha.full$data.type,
  levels=c("cpdag", 'moral', 'fdr', 'sensitivity')
)
mu.full$data.type <- factor(
  mu.full$data.type,
  levels=c("cpdag", 'moral', 'fdr', 'sensitivity')
)
alpha.plot <- ggplot(melted.alpha.full, aes(x=value, color=alpha, fill=alpha)) +
  geom_histogram(position='identity', alpha=0.5, bins=20) + 
  geom_vline(data=mu.full, aes(xintercept=value, color=alpha),
             linetype="dashed") +
  facet_grid(
    variable ~ data.type,
    scales='free_x'
  )
alpha.plot

# Plotting TABU variables
TABU.CPDAG.PATH <- './data/results/simulated_gaussian/tabu_validation//shd_cpdag_runs_50.csv'
TABU.MORAL.PATH <- './data/results/simulated_gaussian/tabu_validation/shd_moral_runs_50.csv'
TABU.FDR.PATH <- './data/results/simulated_gaussian/tabu_validation/fdr_runs_50.csv'
TABU.SENSITIVITY.PATH <- './data/results/simulated_gaussian/tabu_validation/sensitivity_runs_50.csv'
TABU.PATHS <- list(
  'cpdag'=TABU.CPDAG.PATH, 
  'moral'=TABU.MORAL.PATH,
  'fdr'=TABU.FDR.PATH,
  'sensitivity'=TABU.SENSITIVITY.PATH
)
melted.tabu.full <- NULL
mu.full <- NULL
for (i in 1:length(TABU.PATHS)){
  path <- TABU.PATHS[[i]]
  name <- names(TABU.PATHS)[i]
  tabu.sim.gaussian <- read.csv(path)
  melted.tabu <- melt(
    tabu.sim.gaussian,
    id.vars=c('tabu.length', 'graph.init'),
    measure.vars=names(tabu.sim.gaussian)[length(names(tabu.sim.gaussian))]
  )
  melted.tabu['tabu.length'] = lapply(melted.tabu['tabu.length'], as.character)
  melted.tabu$data.type <- name
  mu <- ddply(melted.tabu, c('tabu.length', 'graph.init'), summarise, value=mean(value))
  mu$data.type <- name
  
  if (is.null(melted.tabu.full)){
    melted.tabu.full = melted.tabu
  } else {
    melted.tabu.full = rbind(melted.tabu.full, melted.tabu)
  }
  if (is.null(mu.full)){
    mu.full <- mu
  } else {
    mu.full <- rbind(mu.full, mu)
  }
}
melted.tabu.full$data.type <- factor(
  melted.tabu.full$data.type,
  levels=c("cpdag", 'moral', 'fdr', 'sensitivity')
)
mu.full$data.type <- factor(
  mu.full$data.type,
  levels=c("cpdag", 'moral', 'fdr', 'sensitivity')
)
tabu.plot <- ggplot(melted.tabu.full, aes(x=value, color=tabu.length, fill=tabu.length)) +
  geom_histogram(position='identity', alpha=0.5, bins=20) + 
  geom_vline(data=mu.full, aes(xintercept=value, color=tabu.length),
             linetype="dashed") +
  facet_grid(
    . ~ data.type,
    scales='free_x'
  )
tabu.plot

# Plotting HC variables
HC.CPDAG.PATH <- './data/results/simulated_gaussian/hc_validation//shd_cpdag_runs_50.csv'
HC.MORAL.PATH <- './data/results/simulated_gaussian/hc_validation/shd_moral_runs_50.csv'
HC.FDR.PATH <- './data/results/simulated_gaussian/hc_validation/fdr_runs_50.csv'
HC.SENSITIVITY.PATH <- './data/results/simulated_gaussian/hc_validation/sensitivity_runs_50.csv'
HC.PATHS <- list(
  'cpdag'=HC.CPDAG.PATH, 
  'moral'=HC.MORAL.PATH,
  'fdr'=HC.FDR.PATH,
  'sensitivity'=HC.SENSITIVITY.PATH
)
melted.hc.full <- NULL
mu.full <- NULL
for (i in 1:length(HC.PATHS)){
  path <- HC.PATHS[[i]]
  name <- names(HC.PATHS)[i]
  hc.sim.gaussian <- read.csv(path)
  melted.hc <- melt(
    hc.sim.gaussian,
    id.vars=c('restart.no', 'graph.init'),
    measure.vars=names(hc.sim.gaussian)[length(names(hc.sim.gaussian))]
  )
  melted.hc['restart.no'] = lapply(melted.hc['restart.no'], as.character)
  melted.hc$data.type <- name
  mu <- ddply(melted.hc, c('restart.no', 'graph.init'), summarise, value=mean(value))
  mu$data.type <- name
  
  if (is.null(melted.hc.full)){
    melted.hc.full = melted.hc
  } else {
    melted.hc.full = rbind(melted.hc.full, melted.hc)
  }
  if (is.null(mu.full)){
    mu.full <- mu
  } else {
    mu.full <- rbind(mu.full, mu)
  }
}
melted.hc.full$data.type <- factor(
  melted.hc.full$data.type,
  levels=c("cpdag", 'moral', 'fdr', 'sensitivity')
)
mu.full$data.type <- factor(
  mu.full$data.type,
  levels=c("cpdag", 'moral', 'fdr', 'sensitivity')
)
hc.plot <- ggplot(melted.hc.full, aes(x=value, color=restart.no, fill=restart.no)) +
  geom_histogram(position='identity', alpha=0.5, bins=20) + 
  geom_vline(data=mu.full, aes(xintercept=value, color=restart.no),
             linetype="dashed") +
  facet_grid(
    . ~ data.type,
    scales='free_x'
  )
hc.plot


# bn2igraph <- function(g.bn){
#   g <- igraph.from.graphNEL(as.graphNEL(g.bn))
# }
# 
# ecoli <- generate.network.data(N = 500, network.path = 'bnlearn_networks/ecoli70.rds')
# ecoli.true<- ecoli$network
# ecoli.igraph <- bn2igraph(ecoli.true)
# ecoli.df <- as.data.frame(get.edgelist(ecoli.igraph))
# ecoli.df.g <- graph.data.frame(ecoli.df, directed=TRUE)
# iamb.f <- partial(iamb, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
# inter.g <- inter.iamb.f(ecoli$data)
# par(mfrow = c(1, 2))
# graphviz.compare(cpdag(bn.net(ecoli$network)), inter.g, shape='ellipse', main = c("True DAG", "IAMB DAG"))
# graphviz.compare(moral(cpdag(bn.net(ecoli$network))), moral(inter.g), shape='ellipse', main = c("Moral True DAG", " Moral IAMB DAG"))
# print(shd(cpdag(bn.net(ecoli$network)), inter.g))
# print(shd(moral(cpdag(bn.net(ecoli$network))), moral(inter.g)))