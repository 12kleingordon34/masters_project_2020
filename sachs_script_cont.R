library(bnlearn)
library(fmsb)
library(purrr)
library(reticulate)

use_condaenv("r-reticulate")
source_python('./notears/notears/linear.py')
source_python('./notears/notears/utils.py')
source("./utils.R")

# Load Data
load("./data/sachs/sachs.rda") # Loads true network as `bn`
sachs.df <- read.table("./data/sachs/sachs.data.txt", header=TRUE)


log.sachs.df <- log(sachs.df, base=2)
sachs.bn <- bn.net(bn)

num.bootstraps <- 50
size.bootstrap <- 150
num.permutations <- 100

algo.params <- list(
  'inter.iamb.corr'=list(
    'algorithm'='inter.iamb',
    'algorithm.args'=list(alpha=0.05, test='cor')
  ),
  'inter.iamb.zf'=list(
    'algorithm'='inter.iamb',
    'algorithm.args'=list(alpha=0.05, test='zf')
  ),
  'inter.iamb.mi.g'=list(
    'algorithm'='inter.iamb',
    'algorithm.args'=list(alpha=0.05, test='mi-g')
  ),
  'inter.iamb.mc.cor.100'=list(
    'algorithm'='inter.iamb',
    'algorithm.args'=list(alpha=0.05, test='mc-cor', B=100)
  ),
  'inter.iamb.mc.cor.500'=list(
    'algorithm'='inter.iamb',
    'algorithm.args'=list(alpha=0.05, test='mc-cor', B=500)
  ),
  'inter.iamb.mc.zf.100'=list(
    'algorithm'='inter.iamb',
    'algorithm.args'=list(alpha=0.05, test='mc-zf', B=100)
  ),
  'inter.iamb.mc.zf.500'=list(
    'algorithm'='inter.iamb',
    'algorithm.args'=list(alpha=0.05, test='mc-zf', B=500)
  ),
  'inter.iamb.mc.mi.g.100'=list(
    'algorithm'='inter.iamb',
    'algorithm.args'=list(alpha=0.05, test='mc-mi-g', B=100)
  ),
  'inter.iamb.mc.mi.g.500'=list(
    'algorithm'='inter.iamb',
    'algorithm.args'=list(alpha=0.05, test='mc-mi-g', B=500)
  ),
  'h2pc.cor' = list(
    'algorithm'='h2pc',
    'algorithm.args'=list(
      restrict.args=list(alpha=0.15, test='cor'),
      maximize.args=list(perturb=50, restart=50, score='bge')
    )
  ),
  'h2pc.zf' = list(
    'algorithm'='h2pc',
    'algorithm.args'=list(
      restrict.args=list(alpha=0.15, test='zf'),
      maximize.args=list(perturb=50, restart=50, score='bge')
    )
  ),
  'h2pc.mi.g' = list(
    'algorithm'='h2pc',
    'algorithm.args'=list(
      restrict.args=list(alpha=0.15, test='mi-g'),
      maximize.args=list(perturb=50, restart=50, score='bge')
    )
  ),
  'h2pc.mc.corr.100' = list(
    'algorithm'='h2pc',
    'algorithm.args'=list(
      restrict.args=list(alpha=0.15, test='mc-cor', B=100),
      maximize.args=list(perturb=50, restart=50, score='bge')
    )
  ),
  'h2pc.mc.zf.100' = list(
    'algorithm'='h2pc',
    'algorithm.args'=list(
      restrict.args=list(alpha=0.15, test='mc-zf', B=100),
      maximize.args=list(perturb=50, restart=50, score='bge')
    )
  ),
  'h2pc.mc.corr.500' = list(
    'algorithm'='h2pc',
    'algorithm.args'=list(
      restrict.args=list(alpha=0.15, test='mc-cor', B=500),
      maximize.args=list(perturb=50, restart=50, score='bge')
    )
  ),
  'h2pc.mc.zf.500' = list(
    'algorithm'='h2pc',
    'algorithm.args'=list(
      restrict.args=list(alpha=0.15, test='mc-zf', B=500),
      maximize.args=list(perturb=50, restart=50, score='bge')
    )
  ),
  'h2pc.mc.mi.g.100' = list(
    'algorithm'='h2pc',
    'algorithm.args'=list(
      restrict.args=list(alpha=0.15, test='mc-mi-g', B=100),
      maximize.args=list(perturb=50, restart=50, score='bge')
    )
  ),
  'h2pc.mc.mi.g.500' = list(
    'algorithm'='h2pc',
    'algorithm.args'=list(
      restrict.args=list(alpha=0.15, test='mc-mi-g', B=500),
      maximize.args=list(perturb=50, restart=50, score='bge')
    )
  ),
  'notears.01'=list(
    'algorithm'=partial(notears_linear, loss_type='l2', lambda1=0.1),
    'algorithm.args'=list(hybrid=FALSE)
  ),
  'notears.05'=list(
    'algorithm'=partial(notears_linear, loss_type='l2', lambda1=0.5),
    'algorithm.args'=list(hybrid=FALSE)
  ),
  'hybrid.notears.100'=list(
    'algorithm'=hybrid.notears,
    'algorithm.args'=list(hybrid=TRUE, l1=0.1, alpha=0.15, test='mc-cor', B=100)
  ),
  'hybrid.notears.500'=list(
    'algorithm'=hybrid.notears,
    'algorithm.args'=list(hybrid=TRUE, l1=0.1, alpha=0.15, test='mc-cor', B=500)
  ),
  'hc.0'=list(
    'algorithm'='hc',
    'algorithm.args'=list(score='bge', restart=0, perturb=50)
  ),
  'hc.10'=list(
    'algorithm'='hc',
    'algorithm.args'=list(score='bge', restart=10, perturb=50)
  ),
  'hc.50'=list(
    'algorithm'='hc',
    'algorithm.args'=list(score='bge', restart=50, perturb=50)
  ),
  'tabu.5'=list(
    'algorithm'='tabu',
    'algorithm.args'=list(score='bge', tabu=5)
  ),
  'tabu.10'=list(
    'algorithm'='tabu',
    'algorithm.args'=list(score='bge', tabu=10)
  ),
  'tabu.15'=list(
    'algorithm'='tabu',
    'algorithm.args'=list(score='bge', tabu=15)
  )  
)

averaged.networks <- vector(length(algo.params), mode="list")
for (i in 1:length(algo.params)){
  algo.name <- names(algo.params)[i]
  algo.args <- algo.params[[i]]$algorithm.args
  cat(sprintf(" Time: %s. Algorithm: %s\n", Sys.time(), algo.name))
  if (grepl('notears|lingam', algo.name, fixed=FALSE)){
    boot.graph <- bootstrap.alt(
      algo.params[[i]]$algorithm,
      algo.args,
      log.sachs.df,
      size.bootstrap,
      num.bootstraps
    )
  } else {
    boot.graph <- boot.strength(
      sachs.df,
      R=num.bootstraps,
      m=size.bootstrap,
      algorithm=algo.params[[i]]$algorithm,
      algorithm.args=algo.args
    )
  }
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
  
  forgotten.nodes <- nodes(sachs.bn)[!(nodes(sachs.bn) %in% nodes(learned.graph))]
  for (n in forgotten.nodes) learned.graph = add.node(learned.graph, node=n)
  

  full.shd <- bnlearn::shd(cpdag(learned.graph), cpdag(sachs.bn))
  full.hd <- bnlearn::hamming(cpdag(learned.graph), cpdag(sachs.bn))
  
  full.performance <- calculate.performance.statistics(cpdag(learned.graph), cpdag(sachs.bn))
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
write.csv(results, 'data/sachs_results.csv', row.names=FALSE)

library(ggplot2)
library(ggradar)
library(scales)
library(data.table)

ids <- c('IAMB', 'H2PC', 'NOTEARS', 'HC', 'TABU')
sachs.results <- read.csv('data/sachs_results.csv', row.names = "algorithm")
sachs.results$norm.shd <- sachs.results$full.shd / 20
sachs.results$norm.hd <- sachs.results$full.hd / 20
sachs.results <- sachs.results[, c("norm.shd", "norm.hd", "full.fdr", "full.sensitivity")]

for (id in ids){
  sachs.sub <- sachs.results[rownames(sachs.results) %like% id, ]
  sachs.sub.t <- t(sachs.sub)
  headers <- c(colnames(sachs.sub.t))
  
  sachs.sub.t <- cbind(rownames(sachs.sub.t), data.frame(sachs.sub.t, row.names=NULL))
  sachs.sub.t[,1] = c(
    'Normalised SHD', 'Normalised HD', 'FDR', 'Sensitivity'
  )
  filtered.headers <- unlist(
    lapply(
      headers, 
      function(x) str_match_all(x, '\\(\\s*([\\s\\S]*)\\s*\\)')[[1]][2]
    )
  )
  colnames(sachs.sub.t) <- c('metric', filtered.headers)
  (radar.plot <- ggradar(
    sachs.sub.t, axis.label.size = 3.5,
    grid.label.size = 6, legend.text.size = 10,
    values.radar = c('0', '0.5', '1')
  ))
  ggsave(
    sprintf("./report/figures/radar_cont_%s.pdf", id),
    radar.plot,
    width=9.5,
    height=5.3,
    device='pdf'
  )
}


  # theme_minimal() +
  # theme(text = element_text(size=7), axis.text.y = element_blank())
# log.sachs.results <- read.csv('data/log_sachs_results.csv', row.names='algorithm')
sachs.results$norm.shd <- sachs.results$full.shd / 20
log.sachs.results$norm.shd <- log.sachs.results$full.shd / 20 

sachs.results <- sachs.results[,c('norm.shd', 'full.fdr', 'full.sensitivity')]
log.sachs.results <- log.sachs.results[,c('norm.shd', 'full.fdr', 'full.sensitivity')]


# For discrete models