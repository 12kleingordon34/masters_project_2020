library(bnlearn)
library(ggplot2)
library(purrr)
library(plyr)
library(reshape2)
library(grid)
library(matrixStats)
library(latex2exp)

source('./utils.R')

NETWORK <- 'ecoli'
NUM.RUNS <- 75

ALPHA.CPDAG.PATH <- sprintf('./data/results/simulated_networks/%s/alpha_validation/shd_cpdag_runs_%s.csv', NETWORK, NUM.RUNS)
ALPHA.MORAL.PATH <- sprintf('./data/results/simulated_networks/%s/alpha_validation/shd_moral_runs_%s.csv', NETWORK, NUM.RUNS)
ALPHA.FDR.PATH <- sprintf('./data/results/simulated_networks/%s/alpha_validation/fdr_runs_%s.csv', NETWORK, NUM.RUNS)
ALPHA.SENSITIVITY.PATH <- sprintf('./data/results/simulated_networks/%s/alpha_validation/sensitivity_runs_%s.csv', NETWORK, NUM.RUNS)
ALPHA.PATHS <- list(
  'SHD (CPDAG)'=ALPHA.CPDAG.PATH, 
  'HD (Moralised)'=ALPHA.MORAL.PATH,
  'FDR'=ALPHA.FDR.PATH,
  'Sensitivity'=ALPHA.SENSITIVITY.PATH
)
label.map <- list(
  'iamb.fdr'='IAMB FDR',
  'inter.iamb'='InterIAMB',
  'fast.iamb'='Fast IAMB',
  'iamb'='IAMB',
  'gs'='GS',
  'pc'='PC',
  'N'='N',
  'D'='D',
  'sparsity'='sparsity',
  'alpha'='alpha'
)


# Plotting histograms
melted.alpha.full <- NULL
mu.full <- NULL
for (i in 1:length(ALPHA.PATHS)){
  path <- ALPHA.PATHS[[i]]
  name <- names(ALPHA.PATHS)[i]
  alpha.sim.network <- read.csv(path)
  melted.alpha <- melt(
    alpha.sim.network,
    id.vars=c('alpha', 'D'),
    measure.vars=c('pc', 'gs', 'iamb', 'fast.iamb', 'inter.iamb', 'iamb.fdr')
  )
  for (i in seq(1, length(label.map))){
    melted.alpha$variable <- str_replace_all(
      melted.alpha$variable,
      names(label.map)[i], label.map[i][[1]]
    )
  }
  
  
  if (grepl('HD', name)) melted.alpha$value <- melted.alpha$value / melted.alpha$D 
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
  levels=c("SHD (CPDAG)", 'HD (Moralised)', 'FDR', 'Sensitivity')
)
mu.full$data.type <- factor(
  mu.full$data.type,
  levels=c("SHD (CPDAG)", 'HD (Moralised)', 'FDR', 'Sensitivity')
)
alpha.plot <- ggplot(melted.alpha.full, aes(x=value, fill=alpha)) +
  geom_histogram(position='identity', alpha=0.5, bins=20) + 
  scale_fill_discrete(name = "Alpha") +
  geom_vline(data=mu.full, aes(xintercept=value, color=alpha),
             linetype="dashed", show.legend=FALSE) +
  labs(y='Count', x='')+
  facet_grid(
    variable ~ data.type,
    scales='free_x',
    switch='x'
  ) +
  theme(aspect.ratio=1.2,
        panel.spacing.x=unit(4.5, "mm"), 
        plot.margin=grid::unit(c(0,0,0,0),"mm"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=9))
alpha.plot
ggsave(
  sprintf("./report/figures/alpha_%s_plot.pdf", NETWORK),
  alpha.plot,
  width=8,
  height=12,
  device='pdf'
)

# Plotting TABU variables
TABU.CPDAG.PATH <- sprintf('./data/results/simulated_networks/%s/tabu_validation/shd_cpdag_runs_%s.csv', NETWORK, NUM.RUNS)
TABU.MORAL.PATH <- sprintf('./data/results/simulated_networks/%s/tabu_validation/shd_moral_runs_%s.csv', NETWORK, NUM.RUNS)
TABU.FDR.PATH <- sprintf('./data/results/simulated_networks/%s/tabu_validation/fdr_runs_%s.csv', NETWORK, NUM.RUNS)
TABU.SENSITIVITY.PATH <- sprintf('./data/results/simulated_networks/%s/tabu_validation/sensitivity_runs_%s.csv', NETWORK, NUM.RUNS)
TABU.PATHS <- list(
  'SHD (CPDAG)'=TABU.CPDAG.PATH, 
  'HD (Moralised)'=TABU.MORAL.PATH,
  'FDR'=TABU.FDR.PATH,
  'Sensitivity'=TABU.SENSITIVITY.PATH
)
melted.tabu.full <- NULL
mu.tabu.full <- NULL
mu.graph.full <- NULL
score.mapping <- list('bic-g'='BIC', 'bge'='BGe')

for (i in 1:length(TABU.PATHS)){
  path <- TABU.PATHS[[i]]
  name <- names(TABU.PATHS)[i]
  tabu.sim.network <- read.csv(path)
  melted.tabu <- melt(
    tabu.sim.network,
    id.vars=c('tabu.frac', 'graph.init', 'score', 'D'),
    measure.vars=names(tabu.sim.network)[length(names(tabu.sim.network))]
  )
  for (i in seq(1, length(score.mapping))){
    melted.tabu$score <- str_replace_all(
      melted.tabu$score,
      names(score.mapping)[i], score.mapping[i][[1]]
    )
  }
  
  if (grepl('HD', name)) melted.tabu$value <- melted.tabu$value / melted.tabu$D 
  melted.tabu['tabu.frac'] = lapply(melted.tabu['tabu.frac'], as.character)
  
  melted.tabu$data.type <- name
  mu.tabu <- ddply(melted.tabu, c('tabu.frac','graph.init', 'score'), summarise, value=mean(value))
  mu.tabu$data.type <- name
  mu.graph <- ddply(melted.tabu, c('tabu.frac','graph.init', 'score'), summarise, value=mean(value))
  mu.graph$data.type <- name
  
  
  if (is.null(melted.tabu.full)){
    melted.tabu.full = melted.tabu
  } else {
    melted.tabu.full = rbind(melted.tabu.full, melted.tabu)
  }
  if (is.null(mu.tabu.full)){
    mu.tabu.full <- mu.tabu
    mu.graph.full <- mu.graph
  } else {
    mu.tabu.full <- rbind(mu.tabu.full, mu.tabu)
    mu.graph.full <- rbind(mu.graph.full, mu.graph)
  }
}
melted.tabu.full$data.type <- factor(
  melted.tabu.full$data.type,
  levels=c("SHD (CPDAG)", 'HD (Moralised)', 'FDR', 'Sensitivity')
)
mu.tabu.full$data.type <- factor(
  mu.tabu.full$data.type,
  levels=c("SHD (CPDAG)", 'HD (Moralised)', 'FDR', 'Sensitivity')
)
mu.graph.full$data.type <- factor(
  mu.graph.full$data.type,
  levels=c("SHD (CPDAG)", 'HD (Moralised)', 'FDR', 'Sensitivity')
)

melted.tabu.full <- melted.tabu.full[
  ((melted.tabu.full$value < 2.5) & (melted.tabu.full$data.type =='SHD (CPDAG)'))
  | ((melted.tabu.full$value < 5) & (melted.tabu.full$data.type =='HD (Moralised)'))
  | ((melted.tabu.full$data.type != 'SHD (CPDAG)') & (melted.tabu.full$data.type != 'HD (Moralised)')),
]

tabu.tabu.plot <- ggplot(melted.tabu.full[melted.tabu.full$graph.init=='random',], aes(x=value, fill=tabu.frac)) +
  geom_histogram(position='identity', alpha=0.5, bins=20) + 
  scale_fill_discrete(name = "Tabu Fraction") +
  geom_vline(data=mu.tabu.full[mu.tabu.full$graph.init=='random',], aes(xintercept=value, color=tabu.frac),
             linetype="dashed", show.legend = FALSE) +
  labs(y='Count', x='')+
  facet_grid(
    score ~ data.type,
    scales='free_x',
    switch='x'
  ) +
  theme(aspect.ratio=1.2,
        panel.spacing.x=unit(4.5, "mm"), 
        plot.margin=grid::unit(c(0,0,0,0),"mm"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=9))
tabu.tabu.plot
ggsave(
  sprintf("./report/figures/tabu_tabu_%s_plot.pdf", NETWORK),
  tabu.tabu.plot,
  width=10,
  height=5,
  device='pdf'
)

tabu.graph.plot <- ggplot(melted.tabu.full[melted.tabu.full$tabu.frac == 1, ], aes(x=value, fill=graph.init)) +
  geom_histogram(position='identity', alpha=0.5, bins=20) + 
  scale_fill_discrete(name = "Graph Initialisation") +
  geom_vline(data=mu.graph.full[mu.graph.full$tabu.frac == 1, ], aes(xintercept=value, color=graph.init),
             linetype="dashed", show.legend = FALSE) +
  labs(y='Count', x='')+
  facet_grid(
    score ~ data.type,
    scales='free_x',
    switch='x'
  ) +
  theme(aspect.ratio=1.2,
        panel.spacing.x=unit(4.5, "mm"), 
        plot.margin=grid::unit(c(0,0,0,0),"mm"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=9))
tabu.graph.plot
ggsave(
  sprintf("./report/figures/tabu_graph_%s_plot.pdf", NETWORK),
  tabu.graph.plot,
  width=10,
  height=5,
  device='pdf'
)


# Plotting HC variables
HC.CPDAG.PATH <- sprintf('./data/results/simulated_networks/%s/hc_validation/shd_cpdag_runs_%s.csv', NETWORK, NUM.RUNS)
HC.MORAL.PATH <- sprintf('./data/results/simulated_networks/%s/hc_validation/shd_moral_runs_%s.csv', NETWORK, NUM.RUNS)
HC.FDR.PATH <- sprintf('./data/results/simulated_networks/%s/hc_validation/fdr_runs_%s.csv', NETWORK, NUM.RUNS)
HC.SENSITIVITY.PATH <- sprintf('./data/results/simulated_networks/%s/hc_validation/sensitivity_runs_%s.csv', NETWORK, NUM.RUNS)
HC.PATHS <- list(
  'SHD (CPDAG)'=HC.CPDAG.PATH, 
  'HD (Moralised)'=HC.MORAL.PATH,
  'FDR'=HC.FDR.PATH,
  'Sensitivity'=HC.SENSITIVITY.PATH
)
melted.hc.full <- NULL
mu.restart.full <- NULL
mu.perturb.full <- NULL
mu.graph.init <- NULL
score.mapping <- list('bic-g'='BIC', 'bge'='BGe')


for (i in 1:length(HC.PATHS)){
  path <- HC.PATHS[[i]]
  name <- names(HC.PATHS)[i]
  hc.sim.network <- read.csv(path)
  melted.hc <- melt(
    hc.sim.network,
    id.vars=c('restart.no', 'graph.init', 'perturb', 'score', 'D'),
    measure.vars=names(hc.sim.network)[length(names(hc.sim.network))]
  )
  for (i in seq(1, length(score.mapping))){
    melted.hc$score <- str_replace_all(
      melted.hc$score,
      names(score.mapping)[i], score.mapping[i][[1]]
    )
  }
  
  if (grepl('HD', name)) melted.hc$value <- melted.hc$value / melted.hc$D 
  melted.hc['restart.no'] = lapply(melted.hc['restart.no'], as.character)
  melted.hc['perturb'] = lapply(melted.hc['perturb'], as.character)
  
  melted.hc$data.type <- name
  mu.restart <- ddply(melted.hc, c('perturb', 'restart.no', 'graph.init', 'score'), summarise, value=mean(value))
  mu.restart$data.type <- name
  mu.perturb <- ddply(melted.hc, c('perturb', 'restart.no', 'graph.init', 'score'), summarise, value=mean(value))
  mu.perturb$data.type <- name
  mu.graph <- ddply(melted.hc, c('score', 'graph.init'), summarise, value=mean(value))
  mu.graph$data.type <- name
  
  if (is.null(melted.hc.full)){
    melted.hc.full = melted.hc
  } else {
    melted.hc.full = rbind(melted.hc.full, melted.hc)
  }
  if (is.null(mu.restart.full)){
    mu.restart.full <- mu.restart
    mu.perturb.full <- mu.perturb
    mu.graph.full <- mu.graph
  } else {
    mu.restart.full <- rbind(mu.restart.full, mu.restart)
    mu.perturb.full <- rbind(mu.perturb.full, mu.perturb)
    mu.graph.full <- rbind(mu.graph.full, mu.graph)
  }
}

melted.hc.full <- melted.hc.full[
  ((melted.hc.full$value < 300) & (melted.hc.full$data.type =='SHD (CPDAG)'))
  | ((melted.hc.full$value < 300) & (melted.hc.full$data.type =='HD (Moralised)'))
  | ((melted.hc.full$data.type != 'SHD (CPDAG)') & (melted.hc.full$data.type != 'HD (Moralised)')),
]

melted.hc.full$data.type <- factor(
  melted.hc.full$data.type,
  levels=c("SHD (CPDAG)", 'HD (Moralised)', 'FDR', 'Sensitivity')
)
mu.restart.full$data.type <- factor(
  mu.restart.full$data.type,
  levels=c("SHD (CPDAG)", 'HD (Moralised)', 'FDR', 'Sensitivity')
)
mu.perturb.full$data.type <- factor(
  mu.perturb.full$data.type,
  levels=c("SHD (CPDAG)", 'HD (Moralised)', 'FDR', 'Sensitivity')
)
mu.graph.full$data.type <- factor(
  mu.graph.full$data.type,
  levels=c("SHD (CPDAG)", 'HD (Moralised)', 'FDR', 'Sensitivity')
)

hc.restart.plot <- ggplot(melted.hc.full[(melted.hc.full$perturb == 50)&(melted.hc.full$graph.init == 'random'),], aes(x=value, fill=restart.no)) +
  geom_histogram(position='identity', alpha=0.4, bins=20) + 
  scale_fill_discrete(name = "Restarts", breaks=c("1", "5", "10", "25", "50")) +
  geom_vline(data=mu.restart.full[(mu.restart.full$perturb == 50)&(mu.restart.full$graph.init == 'random'),], aes(xintercept=value, color=restart.no),
             linetype="dashed", show.legend = FALSE) +
  facet_grid(
    score ~ data.type,
    scales='free_x',
    switch='x'
  ) +
  labs(x='', y='Count') +
  theme(aspect.ratio=1.2,
        panel.spacing.x=unit(4.5, "mm"), 
        plot.margin=grid::unit(c(0,0,0,0),"mm"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=9))

hc.perturb.plot <- ggplot(melted.hc.full[(melted.hc.full$restart.no == 50)&(melted.hc.full$graph.init == 'random'),], aes(x=value, fill=perturb)) +
  geom_histogram(position='identity', alpha=0.3, bins=20) + 
  scale_fill_discrete(name = "Perturbations", breaks=c("1", "5", "10", "25", "50")) +
  geom_vline(data=mu.perturb.full[(mu.perturb.full$restart.no == 50)&(mu.perturb.full$graph.init == 'random'),], aes(xintercept=value, color=perturb),
             linetype="dashed", show.legend=FALSE) +
  labs(x='', y='Count') +
  facet_grid(
    score ~ data.type,
    scales='free_x',
    switch='x'
  ) +
  theme(aspect.ratio=1.2,
        panel.spacing.x=unit(4.5, "mm"), 
        plot.margin=grid::unit(c(0,0,0,0),"mm"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=9))

hc.graph.plot <- ggplot(melted.hc.full, aes(x=value, fill=graph.init)) +
  geom_histogram(position='identity', alpha=0.4, bins=20) + 
  scale_fill_discrete(name = "Graph Init.") +
  geom_vline(data=mu.graph.full, aes(xintercept=value, color=graph.init),
             linetype="dashed", show.legend=FALSE) +
  labs(x='Value', y='Count') +
  facet_grid(
    score ~ data.type,
    scales='free_x'
    # space='free'
  ) +
  theme(aspect.ratio=1.2,
        panel.spacing.x=unit(4.5, "mm"), 
        plot.margin=grid::unit(c(0,0,0,0),"mm"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=9))

ggsave(
  sprintf("./report/figures/hc_restart_%s_plot.pdf", NETWORK),
  hc.restart.plot,
  width=10,
  height=5,
  device='pdf'
)
ggsave(
  sprintf("./report/figures/hc_perturb_%s_plot.pdf", NETWORK),
  hc.perturb.plot,
  width=10,
  height=5,
  device='pdf'
)
ggsave(
  sprintf("./report/figures/hc_graph_%s_plot.pdf", NETWORK),
  hc.graph.plot,
  width=10,
  height=5,
  device='pdf'
)


# Plotting NOTEARS variables
NT.CPDAG.PATH <- sprintf('./data/results/simulated_networks/%s/notears_validation/shd_runs_full.csv', NETWORK)
NT.MORAL.PATH <- sprintf('./data/results/simulated_networks/%s/notears_validation/hd_runs_full.csv', NETWORK)
NT.FDR.PATH <- sprintf('./data/results/simulated_networks/%s/notears_validation/fdr_runs_full.csv', NETWORK)
NT.SENSITIVITY.PATH <- sprintf('./data/results/simulated_networks/%s/notears_validation/sensitivity_runs_full.csv', NETWORK)
NT.PATHS <- list(
  'SHD (CPDAG)'=NT.CPDAG.PATH, 
  'HD (Moralised)'=NT.MORAL.PATH,
  'FDR'=NT.FDR.PATH,
  'Sensitivity'=NT.SENSITIVITY.PATH
)
melted.nt.full <- NULL
mu.full <- NULL
for (i in 1:length(NT.PATHS)){
  path <- NT.PATHS[[i]]
  name <- names(NT.PATHS)[i]
  nt.sim.network <- read.csv(path)
  # # Num hyperparams tested
  # num.hyp <- max(table(nt.sim.gaussian$seed))
  # # Count of seeds used
  # t = table(nt.sim.gaussian['seed'])
  # # Analyse seeds which have been assessed by all hyperparam choices
  # nt.sim.network <-nt.sim.network[
  #   is.element(nt.sim.network$seed, names(t[t==num.hyp])),
  # ]
  melted.nt <- melt(
    nt.sim.network,
    id.vars=c('l1', 'D'),
    measure.vars=names(nt.sim.network)[length(names(nt.sim.network))]
  )
  if (grepl('HD', name)) melted.nt$value <- melted.nt$value / melted.nt$D 
  melted.nt['l1'] = lapply(melted.nt['l1'], as.character)
  melted.nt$data.type <- name
  mu <- ddply(melted.nt, c('l1'), summarise, value=mean(value))
  mu$data.type <- name
  
  if (is.null(melted.nt.full)){
    melted.nt.full = melted.nt
  } else {
    melted.nt.full = rbind(melted.nt.full, melted.nt)
  }
  if (is.null(mu.full)){
    mu.full <- mu
  } else {
    mu.full <- rbind(mu.full, mu)
  }
}
# Filter for extreme values (purely for graphical reasons)
melted.nt.full$data.type <- factor(
  melted.nt.full$data.type,
  levels=c("SHD (CPDAG)", 'HD (Moralised)', 'FDR', 'Sensitivity')
)
mu.full$data.type <- factor(
  mu.full$data.type,
  levels=c("SHD (CPDAG)", 'HD (Moralised)', 'FDR', 'Sensitivity')
)
melted.nt.full <- melted.nt.full[
  ((melted.nt.full$value < 6) & (melted.nt.full$data.type =='SHD (CPDAG)'))
  | ((melted.nt.full$value < 12) & (melted.nt.full$data.type =='HD (Moralised)'))
  | ((melted.nt.full$data.type != 'SHD (CPDAG)') & (melted.nt.full$data.type != 'HD (Moralised)')),
]
nt.plot <- ggplot(melted.nt.full, aes(x=value, fill=l1)) +
  geom_histogram(position='identity', alpha=0.4, bins=20) + 
  scale_fill_discrete(name = "L1 Regularisation") +
  geom_vline(data=mu.full, aes(xintercept=value, color=l1),
             linetype="dashed", show.legend=FALSE) +
  facet_grid(
    . ~ data.type,
    scales='free_x',
    switch='x'
  ) +
  labs(x='', y='Count') +
  theme(aspect.ratio=1.2,
        panel.spacing.x=unit(4.5, "mm"), 
        plot.margin=grid::unit(c(0,0,0,0),"mm"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=9))
ggsave(
  sprintf("./report/figures/nt_%s_plot.pdf", NETWORK), 
  nt.plot,
  width=10,
  height=5,
  device='pdf'
)


# Hybrid Plots
HYBRID.CPDAG.PATH <- sprintf('./data/results/simulated_networks/%s/hybrid/shd_cpdag_runs_%s.csv', NETWORK, NUM.RUNS)
HYBRID.MORAL.PATH <- sprintf('./data/results/simulated_networks/%s/hybrid/shd_moral_runs_%s.csv', NETWORK, NUM.RUNS)
HYBRID.FDR.PATH <- sprintf('./data/results/simulated_networks/%s/hybrid/fdr_runs_%s.csv', NETWORK, NUM.RUNS)
HYBRID.SENSITIVITY.PATH <- sprintf('./data/results/simulated_networks/%s/hybrid/sensitivity_runs_%s.csv', NETWORK, NUM.RUNS)
HYBRID.PATHS <- list(
  'SHD (CPDAG)'=HYBRID.CPDAG.PATH, 
  'HD (Moralised)'=HYBRID.MORAL.PATH,
  'FDR'=HYBRID.FDR.PATH,
  'Sensitivity'=HYBRID.SENSITIVITY.PATH
)
hybrid.label.map <- list(
  'hybrid.iamb.hc'='InterIAMB with HC',
  'hybrid.h2pc.hc'='H2PC',
  'hybrid.mmhc.hc'='MMHC',
  'hybrid.hiton.hc'='SI-HITON-PC'
)


# Plotting histograms
melted.hybrid.full <- NULL
mu.full <- NULL
for (i in 1:length(HYBRID.PATHS)){
  path <- HYBRID.PATHS[[i]]
  name <- names(HYBRID.PATHS)[i]
  hybrid.sim.network <- read.csv(path)
  melted.hybrid <- melt(
    hybrid.sim.network,
    id.vars=c('alpha', 'D'),
    measure.vars=names(hybrid.sim.network)[5:length(hybrid.sim.network)]
  )
  for (i in seq(1, length(hybrid.label.map))){
    melted.hybrid$variable <- str_replace_all(
      melted.hybrid$variable,
      names(hybrid.label.map)[i], hybrid.label.map[i][[1]]
    )
  }
  
  if (grepl('HD', name)) melted.hybrid$value <- melted.hybrid$value / melted.hybrid$D 
  melted.hybrid['alpha'] = lapply(melted.hybrid['alpha'], as.character)
  melted.hybrid$data.type <- name
  mu <- ddply(melted.hybrid, c('alpha', 'variable'), summarise, value=mean(value))
  mu$data.type <- name
  
  if (is.null(melted.hybrid.full)){
    melted.hybrid.full = melted.hybrid
  } else {
    melted.hybrid.full = rbind(melted.hybrid.full, melted.hybrid)
  }
  if (is.null(mu.full)){
    mu.full <- mu
  } else {
    mu.full <- rbind(mu.full, mu)
  }
}
melted.hybrid.full$data.type <- factor(
  melted.hybrid.full$data.type,
  levels=c("SHD (CPDAG)", 'HD (Moralised)', 'FDR', 'Sensitivity')
)
mu.full$data.type <- factor(
  mu.full$data.type,
  levels=c("SHD (CPDAG)", 'HD (Moralised)', 'FDR', 'Sensitivity')
)
hybrid.plot <- ggplot(melted.hybrid.full, aes(x=value, fill=alpha)) +
  geom_histogram(position='identity', alpha=0.5, bins=20) + 
  scale_fill_discrete(name = "Alpha") +
  geom_vline(data=mu.full, aes(xintercept=value, color=alpha),
             linetype="dashed", show.legend=FALSE) +
  labs(y='Count', x='')+
  facet_grid(
    variable ~ data.type,
    scales='free_x',
    switch='x'
  ) +
  theme(aspect.ratio=1.2,
        panel.spacing.x=unit(4.5, "mm"), 
        plot.margin=grid::unit(c(0,0,0,0),"mm"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=9))
hybrid.plot
ggsave(
  sprintf("./report/figures/hybrid_%s_plot.pdf", NETWORK),
  hybrid.plot,
  width=8,
  height=8,
  device='pdf'
)

# Data Size Plots
DATA.SIZE.SHD <- sprintf('./data/results/simulated_networks/%s/data_size/shd_full.csv', NETWORK)
DATA.SIZE.HD <- sprintf('./data/results/simulated_networks/%s/data_size/hd_full.csv', NETWORK)
DATA.SIZE.FDR <- sprintf('./data/results/simulated_networks/%s/data_size/fdr_full.csv', NETWORK)
DATA.SIZE.SENSITIVITY <- sprintf('./data/results/simulated_networks/%s/data_size/sensitivity_full.csv', NETWORK)

DATA.SIZE.PATHS <- list(
  'SHD (CPDAG)'=DATA.SIZE.SHD, 
  'HD (Moralised)'=DATA.SIZE.HD,
  'FDR'=DATA.SIZE.FDR,
  'Sensitivity'=DATA.SIZE.SENSITIVITY
)
melted.data.size.full <- NULL
mu.full <- NULL
for (i in 1:length(DATA.SIZE.PATHS)){
  path <- DATA.SIZE.PATHS[[i]]
  name <- names(DATA.SIZE.PATHS)[i]
  data.sim.network <- read.csv(path)
  melted.data <- melt(
    data.sim.network,
    id.vars=c('N', 'D'),
    measure.vars=c('gs', 'inter.iamb', 'h2pc', 'mmhc', 'si.hiton', 'notears', 'hybrid.notears.f', 'hc')
  )
  if (grepl('HD', name)) melted.data$value <- melted.data$value / melted.data$D 
  melted.data['N'] = lapply(melted.data['N'], as.character)
  melted.data$data.type <- name
  mu <- ddply(melted.data, c('N', 'variable'), summarise, value=mean(value))
  mu$data.type <- name
  
  if (is.null(melted.data.size.full)){
    melted.data.size.full = melted.data
  } else {
    melted.data.size.full = rbind(melted.data.size.full, melted.data)
  }
  if (is.null(mu.full)){
    mu.full <- mu
  } else {
    mu.full <- rbind(mu.full, mu)
  }
}
melted.data.size.full$data.type <- factor(
  melted.data.size.full$data.type,
  levels=c("SHD (CPDAG)", 'HD (Moralised)', 'FDR', 'Sensitivity')
)
mu.full$data.type <- factor(
  mu.full$data.type,
  levels=c("SHD (CPDAG)", 'HD (Moralised)', 'FDR', 'Sensitivity')
)
melted.data.size.full$N <- factor(
  melted.data.size.full$N,
  levels=c("50", '100', '250', '500', '1000')
)
mu.full$N <- factor(
  mu.full$N,
  levels=c("50", '100', '250', '500', '1000')
)

melted.data.size.full <- melted.data.size.full[
  ((melted.data.size.full$value < 1.5) & (melted.data.size.full$data.type =='SHD (CPDAG)'))
  | ((melted.data.size.full$value < 1.5) & (melted.data.size.full$data.type =='HD (Moralised)'))
  | ((melted.data.size.full$value < 0.5) & (melted.data.size.full$data.type =='FDR'))
  | ((melted.data.size.full$data.type != 'SHD (CPDAG)') & (melted.data.size.full$data.type != 'HD (Moralised)') & (melted.data.size.full$data.type != 'FDR')),
]
mu.full <- mu.full[
  ((mu.full$value < 1.5) & (mu.full$data.type =='SHD (CPDAG)'))
  | ((mu.full$value < 1.5) & (mu.full$data.type =='HD (Moralised)'))
  | ((mu.full$value < 0.5) & (mu.full$data.type =='FDR'))
  | ((mu.full$data.type != 'SHD (CPDAG)') & (mu.full$data.type != 'HD (Moralised)') & (mu.full$data.type != 'FDR')),
]
(data.size.plot <- ggplot(melted.data.size.full, aes(x=value, fill=variable)) +
  geom_histogram(position='identity', alpha=0.4, bins=30) + 
  scale_fill_discrete(name = "Algorithm", labels=c("GS", "InterIAMB", "H2PC", "MMHC", "SI-HITON-PC", "NOTEARS", "Hybrid NOTEARS", "HC")) +
  geom_vline(data=mu.full, aes(xintercept=value, color=variable),
             linetype="dashed", show.legend=FALSE) +
  facet_grid(
    N ~ data.type,
    scales='free_x',
    switch='x'
  ) +
  labs(y='Count', x='')+
  theme(aspect.ratio=1.2,
        panel.spacing.x=unit(4.5, "mm"), 
        plot.margin=grid::unit(c(0,0,0,0),"mm"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=9)))


ggsave(
  sprintf("./report/figures/data_size_%s_plot.pdf", NETWORK),
  data.size.plot,
  width=8,
  height=9,
  device='pdf'
)

# Compare random DAGs to biological networks
arth <- readRDS('bnlearn_networks/arth150.rds')
ecoli <- readRDS('bnlearn_networks/ecoli70.rds')
ecoli.degrees <- c(); arth.degrees <- c()
for (name in names(ecoli)) ecoli.degrees <- c(ecoli.degrees, bnlearn::degree(ecoli, name))
for (name in names(arth)) arth.degrees <- c(arth.degrees, bnlearn::degree(arth, name))
bio.degrees <- c(arth.degrees, ecoli.degrees)
df.1 <- data.frame(matrix(unlist(bio.degrees)))
df.1$d.frac <- "bio.networks"
df.2 <- NULL
for (d.frac in c(1, 1.5, 3)){
  ran.degrees <- c()
  for (i in seq(1, 5)){
    D <- 75
    ran = generate.dag(1, D, d.frac/D, i)
    for (name in nodes(ran$graph)){
      ran.degrees <- c(ran.degrees, bnlearn::degree(ran$graph, name))
    }
  }
  if (is.null(df.2)) {
    df.2 <- data.frame(matrix(unlist(ran.degrees)))
    df.2$d.frac <- paste('random', as.character(d.frac), sep='.')
  } else {
    temp <- data.frame(matrix(unlist(ran.degrees)))
    temp$d.frac <- paste('random', as.character(d.frac), sep='.')
    df.2 <- rbind(df.2, temp)
  }
}
colnames(df.1) <- c("bio.networks", "d.frac")
colnames(df.2) <- c("random.networks", "d.frac")
full.df <- rbind(melt(df.1), melt(df.2))

full.df <- full.df[full.df$value <= 10, ]

sparsity.plot <- ggplot(full.df, aes(x=value,fill=d.frac))+
  geom_histogram(aes(y=0.5*..density..),
                 alpha=1,position='dodge2',binwidth=0.5) +
  scale_fill_discrete(
    name='Network', 
    labels=unname(TeX(c(
      'Biological','Random: $\\frac{1}{D}$', 
      'Random: $\\frac{3}{2D}$',
      'Random: $\\frac{3}{D}$'
    )))
  ) + labs(x='Degree',y="Normalised Count")
sparsity.plot
ggsave(
  sprintf("./report/figures/sparsity_comparison.pdf", NETWORK),
  sparsity.plot,
  width=8,
  height=5,
  device='pdf'
)

# Confusion Matrix
label.map <- list(
  'hybrid.notears.f'='Hybrid NOTEARS',
  'notears'='NOTEARS',
  'iamb.hc'='SI-HITON-PC',
  'mmhc'='MMHC',
  'h2pc'='H2PC',
  'iamb.fdr'='IAMB FDR',
  'inter.iamb'='InterIAMB',
  'fast.iamb'='Fast IAMB',
  'iamb'='IAMB',
  'gs'='GS',
  'pc'='PC'
)

conf.means <- read.table(
  './data/results/simulated_networks/ecoli/confusion_matrices/shd_rel_500.csv',
  header=TRUE,
  row.names = 1,
  sep = ','
)
melted.conf.means <- melt(as.matrix(conf.means))
colnames(melted.conf.means) <- c('algo.1', 'algo.2', 'mean.shd')
melted.conf.means$int.mean <- round(melted.conf.means$mean.shd, 0)
for (i in seq(1, length(label.map))){
  for (col in c('algo.1', 'algo.2')){
    melted.conf.means[, col] = gsub(names(label.map)[i], label.map[[i]],melted.conf.means[, col])
  }
}

conf.plot <- ggplot(data=melted.conf.means, mapping=aes(x=algo.1, y=algo.2, fill=mean.shd)) +
  geom_tile() +
  geom_text(aes(label=int.mean), vjust=0.5, fontface='bold', alpha=1) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x=element_blank(),
    axis.title.y=element_blank()
  )
ggsave(
  sprintf("./report/figures/confusion_%s.pdf", NETWORK),
  conf.plot,
  width=10,
  height=8,
  device='pdf'
)

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