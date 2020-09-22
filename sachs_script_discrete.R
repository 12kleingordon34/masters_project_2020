library(bnlearn)
library(fmsb)
library(purrr)

source("./utils.R")

# Load Data
load("./data/sachs/sachs.rda") # Loads true network as `bn`
sachs.df <- read.table("./data/sachs/sachs.data.txt", header=TRUE)
sachs.dis <- discretize(
  sachs.df, method = "hartemink",breaks = 5, ibreaks = 100, idisc = "quantile"
)

sachs.bn <- bn.net(bn)

num.bootstraps <- 50
num.permutations <- 100

algo.params <- list(
  'inter.iamb.mi'=list(
    'algorithm'='inter.iamb',
    'algorithm.args'=list(alpha=0.05, test='mi')
  ),
  'inter.iamb.mi.sh'=list(
    'algorithm'='inter.iamb',
    'algorithm.args'=list(alpha=0.05, test='mi-sh')
  ),
  'inter.iamb.x2'=list(
    'algorithm'='inter.iamb',
    'algorithm.args'=list(alpha=0.05, test='x2')
  ),
  'h2pc.mi' = list(
    'algorithm'='h2pc',
    'algorithm.args'=list(
      restrict.args=list(alpha=0.15, test='mi'),
      maximize.args=list(perturb=4, restart=25, score='bde')
    )
  ),
  'h2pc.mi.sh' = list(
    'algorithm'='h2pc',
    'algorithm.args'=list(
      restrict.args=list(alpha=0.15, test='mi-sh'),
      maximize.args=list(perturb=4, restart=25, score='bde')
    )
  ),
  'h2pc.x2' = list(
    'algorithm'='h2pc',
    'algorithm.args'=list(
      restrict.args=list(alpha=0.15, test='x2'),
      maximize.args=list(perturb=4, restart=25, score='bde')
    )
  ),
  'hc.0'=list(
    'algorithm'='hc',
    'algorithm.args'=list(score='bde', restart=0, perturb=4)
  ),
  'hc.10'=list(
    'algorithm'='hc',
    'algorithm.args'=list(score='bde', restart=10, perturb=4)
  ),
  'tabu.5'=list(
    'algorithm'='tabu',
    'algorithm.args'=list(score='bde', tabu=5)
  ),
  'tabu.10'=list(
    'algorithm'='tabu',
    'algorithm.args'=list(score='bde', tabu=10)
  )
  # 'tabu.15'=list(
  #   'algorithm'='tabu',
  #   'algorithm.args'=list(score='bde', tabu=15)
  # )  
)

averaged.networks <- vector(length(algo.params), mode="list")
for (i in 1:length(algo.params)){
  algo.name <- names(algo.params)[i]
  algo.args <- algo.params[[i]]$algorithm.args
  cat(sprintf(" Time: %s. Algorithm: %s\n", Sys.time(), algo.name))
  boot.graph <- boot.strength(
    sachs.dis,
    R=num.bootstraps,
    m=size.bootstrap,
    algorithm=algo.params[[i]]$algorithm,
    algorithm.args=algo.args
  )
  averaged.networks[[i]] <- averaged.network(boot.graph)
}

# For each network
results <- array(numeric(), c(0, 5))
colnames(results) <- c(
  'algorithm', 
  'full.shd', 
  'full.hd', 
  'full.fdr', 
  'full.sensitivity'
)
for (i in seq(1, length(averaged.networks))){
  learned.graph <- averaged.networks[[i]]
  relevant.nodes = nodes(learned.graph)[sapply(nodes(learned.graph), degree, object = learned.graph) > 0]
  sachs.subgraph <- subgraph(sachs.bn, relevant.nodes); learned.subgraph <- subgraph(learned.graph, relevant.nodes)

  full.shd <- bnlearn::shd(cpdag(learned.subgraph), cpdag(sachs.subgraph))
  full.hd <- bnlearn::hamming(cpdag(learned.subgraph), cpdag(sachs.subgraph))
  full.performance <- calculate.performance.statistics(cpdag(learned.subgraph), cpdag(sachs.subgraph))
  full.fdr <- full.performance$fdr; full.sensitivity <- full.performance$sensitivity
  print(c(
    names(algo.params)[i],
    full.shd,
    full.hd,
    full.fdr,
    full.sensitivity
  ))
  results <- rbind(
    results, 
    c(
      names(algo.params)[i],
      full.shd,
      full.hd,
      full.fdr,
      full.sensitivity
    )
  )
}
write.csv(results, 'data/discrete_sachs_results.csv', row.names=FALSE)


library(ggplot2)
library(ggradar)
library(scales)
library(data.table)


sachs.results <- read.csv('data/discrete_sachs_results.csv', row.names='algorithm')
sachs.results$norm.shd <- sachs.results$full.shd / 20
sachs.results$norm.hd <- sachs.results$full.hd / 20
sachs.results <- sachs.results[, c("norm.shd", "norm.hd", "full.fdr", "full.sensitivity")]

sachs.sub <- sachs.results[rownames(sachs.results) %like% id, ]
sachs.results.t <- t(sachs.results)
headers <- c(colnames(sachs.results.t))

sachs.results.t <- cbind(rownames(sachs.results.t), data.frame(sachs.results.t, row.names=NULL))
sachs.results.t[,1] = c(
  'Normalised SHD', 'Normalised HD', 'FDR', 'Sensitivity'
)
colnames(sachs.results.t) <- c('metric', headers)
(radar.plot <- ggradar(
  sachs.results.t, axis.label.size = 3.5, grid.label.size = 6,
  legend.text.size = 10, grid.max=0.6, grid.mid=0.3,
  values.radar = c('0', '0.3', '0.6')
))
ggsave(
  "./report/figures/radar_dis.pdf",
  radar.plot,
  width=9.5,
  height=5.3,
  device='pdf'
)
