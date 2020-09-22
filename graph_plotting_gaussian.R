library(bnlearn)
library(ggplot2)
library(purrr)
library(plyr)
library(reshape2)
library(dplyr)
library(grid)

sample.number <- 75
ALPHA.CPDAG.PATH <- sprintf('./data/results/simulated_gaussian/alpha_validation/shd_runs_%s.csv', sample.number)
ALPHA.MORAL.PATH <- sprintf('./data/results/simulated_gaussian/alpha_validation/hd_runs_%s.csv', sample.number)
ALPHA.FDR.PATH <- sprintf('./data/results/simulated_gaussian/alpha_validation/fdr_runs_%s.csv', sample.number)
ALPHA.SENSITIVITY.PATH <- sprintf('./data/results/simulated_gaussian/alpha_validation/sensitivity_runs_%s.csv', sample.number)
ALPHA.PATHS <- list(
  'SHD (CPDAG)'=ALPHA.CPDAG.PATH, 
  'SHD (Moralised)'=ALPHA.MORAL.PATH,
  'FDR'=ALPHA.FDR.PATH,
  'Sensitivity'=ALPHA.SENSITIVITY.PATH
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
    id.vars=c('alpha', 'D'),
    measure.vars=c('pc', 'gs', 'iamb', 'fast.iamb', 'inter.iamb', 'iamb.fdr')
  )
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
  levels=c("SHD (CPDAG)", 'SHD (Moralised)', 'FDR', 'Sensitivity')
)
mu.full$data.type <- factor(
  mu.full$data.type,
  levels=c("SHD (CPDAG)", 'SHD (Moralised)', 'FDR', 'Sensitivity')
)
alpha.plot <- ggplot(melted.alpha.full, aes(x=value, fill=alpha)) +
  geom_histogram(position='identity', alpha=0.3, bins=20) + 
  scale_fill_discrete(name = "Alpha") +
  geom_vline(data=mu.full, aes(xintercept=value, color=alpha),
             linetype="dashed", show.legend=FALSE) +
  labs(y='Count', x='Value')+
  facet_grid(
    variable ~ data.type,
    scales='free_x'
  ) +
  theme(aspect.ratio=1.2,
        panel.spacing.x=unit(4.5, "mm"), 
        plot.margin=grid::unit(c(0,0,0,0),"mm"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=9))

ggsave(
  "./report/figures/alpha_gaussian_plot.pdf",
  alpha.plot,
  width=8,
  height=12,
  device='pdf'
)


# Plotting TABU variables
TABU.CPDAG.PATH <- sprintf('./data/results/simulated_gaussian/tabu_validation/shd_runs_%s.csv', sample.number)
TABU.MORAL.PATH <- sprintf('./data/results/simulated_gaussian/tabu_validation/hd_runs_%s.csv', sample.number)
TABU.FDR.PATH <- sprintf('./data/results/simulated_gaussian/tabu_validation/fdr_runs_%s.csv', sample.number)
TABU.SENSITIVITY.PATH <- sprintf('./data/results/simulated_gaussian/tabu_validation/sensitivity_runs_%s.csv', sample.number)
TABU.PATHS <- list(
  'SHD (CPDAG)'=TABU.CPDAG.PATH, 
  'SHD (Moralised)'=TABU.MORAL.PATH,
  'FDR'=TABU.FDR.PATH,
  'Sensitivity'=TABU.SENSITIVITY.PATH
)
melted.tabu.tabu.full <- NULL
mu.tabu.full <- NULL
melted.tabu.graph.full <- NULL
mu.graph.full <- NULL
for (i in 1:length(TABU.PATHS)){
  path <- TABU.PATHS[[i]]
  name <- names(TABU.PATHS)[i]
  tabu.sim.gaussian <- read.csv(path)
  melted.tabu <- melt(
    tabu.sim.gaussian,
    id.vars=c('tabu.frac', 'graph.init', 'score', 'D'),
    measure.vars=names(tabu.sim.gaussian)[length(names(tabu.sim.gaussian))]
  )
  if (grepl('HD', name)) melted.tabu$value <- melted.tabu$value / melted.tabu$D 
  melted.tabu['tabu.frac'] = lapply(melted.tabu['tabu.frac'], as.character)
  
  melted.tabu$data.type <- name
  mu.tabu <- ddply(melted.tabu, c('graph.init', 'tabu.frac', 'score'), summarise, value=mean(value))
  mu.tabu$data.type <- name
  mu.graph <- ddply(melted.tabu, c('graph.init','tabu.frac', 'score'), summarise, value=mean(value))
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
  levels=c("SHD (CPDAG)", 'SHD (Moralised)', 'FDR', 'Sensitivity')
)
mu.tabu.full$data.type <- factor(
  mu.tabu.full$data.type,
  levels=c("SHD (CPDAG)", 'SHD (Moralised)', 'FDR', 'Sensitivity')
)
mu.graph.full$data.type <- factor(
  mu.graph.full$data.type,
  levels=c("SHD (CPDAG)", 'SHD (Moralised)', 'FDR', 'Sensitivity')
)
melted.tabu.full <- melted.tabu.full[
  ((melted.tabu.full$value < 4) & (melted.tabu.full$data.type =='SHD (CPDAG)'))
  | ((melted.tabu.full$value < 15) & (melted.tabu.full$data.type =='SHD (Moralised)'))
  | ((melted.tabu.full$data.type != 'SHD (CPDAG)') & (melted.tabu.full$data.type != 'SHD (Moralised)')),
]
tabu.tabu.plot <- ggplot(melted.tabu.full[melted.tabu.full$graph.init=='random',], aes(x=value, fill=tabu.frac)) +
  geom_histogram(position='identity', alpha=0.5, bins=20) + 
  scale_fill_discrete(name = "Tabu Fraction") +
  geom_vline(data=mu.tabu.full[mu.tabu.full$graph.init=='random',], aes(xintercept=value, color=tabu.frac),
             linetype="dashed", show.legend = FALSE) +
  labs(x='Value', y='Count') +
  facet_grid(
    score ~ data.type,
    scales='free_x'
  ) + theme(aspect.ratio=1.2,
            panel.spacing.x=unit(4.5, "mm"), 
            plot.margin=grid::unit(c(0,0,0,0),"mm"),
            legend.title=element_text(size=11),
            legend.text=element_text(size=9))
ggsave(
  "./report/figures/tabu_tabu_gaussian_plot.pdf",
  tabu.tabu.plot,
  width=10,
  height=5,
  device='pdf'
)
tabu.graph.plot <- ggplot(melted.tabu.full[melted.tabu.full$tabu.frac == 1, ], aes(x=value, fill=graph.init)) +
  geom_histogram(position='identity', alpha=0.5, bins=20) + 
  scale_fill_discrete(name = "Tabu Fraction") +
  geom_vline(data=mu.graph.full[mu.graph.full$tabu.frac == 1, ], aes(xintercept=value, color=graph.init),
             linetype="dashed", show.legend = FALSE) +
  labs(x='Value', y='Count') +
  facet_grid(
    score ~ data.type,
    scales='free_x'
  ) + theme(aspect.ratio=1.2,
            panel.spacing.x=unit(4.5, "mm"), 
            plot.margin=grid::unit(c(0,0,0,0),"mm"),
            legend.title=element_text(size=11),
            legend.text=element_text(size=9))
ggsave(
  "./report/figures/tabu_graph_gaussian_plot.pdf",
  tabu.graph.plot,
  width=10,
  height=5,
  device='pdf'
)


# Plotting HC variables
HC.CPDAG.PATH <-sprintf('./data/results/simulated_gaussian/hc_validation/shd_runs_%s.csv', sample.number)
HC.MORAL.PATH <-sprintf('./data/results/simulated_gaussian/hc_validation/shd_runs_%s.csv', sample.number)
HC.FDR.PATH <- sprintf('./data/results/simulated_gaussian/hc_validation/fdr_runs_%s.csv', sample.number)
HC.SENSITIVITY.PATH <- sprintf('./data/results/simulated_gaussian/hc_validation/sensitivity_runs_%s.csv', sample.number)
HC.PATHS <- list(
  'SHD (CPDAG)'=HC.CPDAG.PATH, 
  'SHD (Moralised)'=HC.MORAL.PATH,
  'FDR'=HC.FDR.PATH,
  'Sensitivity'=HC.SENSITIVITY.PATH
)
melted.hc.full <- NULL
mu.restart.full <- NULL
mu.perturb.full <- NULL
mu.graph.init <- NULL

for (i in 1:length(HC.PATHS)){
  path <- HC.PATHS[[i]]
  name <- names(HC.PATHS)[i]
  hc.sim.gaussian <- read.csv(path)
  melted.hc <- melt(
    hc.sim.gaussian,
    id.vars=c('restart.no', 'graph.init', 'perturb', 'score', 'D'),
    measure.vars=names(hc.sim.gaussian)[length(names(hc.sim.gaussian))]
  )
  if (grepl('HD', name)) melted.hc$value <- melted.hc$value / melted.hc$D 
  melted.hc['restart.no'] = lapply(melted.hc['restart.no'], as.character)
  melted.hc['perturb'] = lapply(melted.hc['perturb'], as.character)
  
  melted.hc$data.type <- name
  mu.restart <- ddply(melted.hc, c('perturb', 'restart.no', 'score'), summarise, value=mean(value))
  mu.restart$data.type <- name
  mu.perturb <- ddply(melted.hc, c('perturb', 'restart.no', 'score'), summarise, value=mean(value))
  mu.perturb$data.type <- name
  mu.graph <- ddply(melted.hc, c('graph.init', 'score'), summarise, value=mean(value))
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
  ((melted.hc.full$value < 3) & (melted.hc.full$data.type =='SHD (CPDAG)'))
  | ((melted.hc.full$value < 3) & (melted.hc.full$data.type =='SHD (Moralised)'))
  | ((melted.hc.full$data.type != 'SHD (CPDAG)') & (melted.hc.full$data.type != 'SHD (Moralised)')),
]

melted.hc.full$data.type <- factor(
  melted.hc.full$data.type,
  levels=c("SHD (CPDAG)", 'SHD (Moralised)', 'FDR', 'Sensitivity')
)
mu.restart.full$data.type <- factor(
  mu.restart.full$data.type,
  levels=c("SHD (CPDAG)", 'SHD (Moralised)', 'FDR', 'Sensitivity')
)
mu.perturb.full$data.type <- factor(
  mu.perturb.full$data.type,
  levels=c("SHD (CPDAG)", 'SHD (Moralised)', 'FDR', 'Sensitivity')
)
mu.graph.full$data.type <- factor(
  mu.graph.full$data.type,
  levels=c("SHD (CPDAG)", 'SHD (Moralised)', 'FDR', 'Sensitivity')
)
hc.restart.plot <- ggplot(melted.hc.full[mu.restart.full$perturb == 50,], aes(x=value, fill=restart.no)) +
  geom_histogram(position='identity', alpha=0.4, bins=20) + 
  scale_fill_discrete(name = "Restarts") +
  geom_vline(data=mu.restart.full[mu.restart.full$perturb == 50,], aes(xintercept=value, color=restart.no),
             linetype="dashed", show.legend = FALSE) +
  facet_grid(
    score ~ data.type,
    scales='free_x'
  ) +
  labs(x='Value', y='Count') +
  theme(aspect.ratio=1.2,
        panel.spacing.x=unit(4.5, "mm"), 
        plot.margin=grid::unit(c(0,0,0,0),"mm"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=9))
hc.restart.plot

hc.perturb.plot <- ggplot(
    melted.hc.full[melted.hc.full$restart.no == 50,], aes(x=value, fill=perturb)
  ) +
  geom_histogram(position='identity', alpha=0.4, bins=20) + 
  scale_fill_discrete(name = "Perturbations") +
  geom_vline(data=mu.perturb.full[mu.perturb.full$restart.no == 50,], aes(xintercept=value, color=perturb),
             linetype="dashed", show.legend=FALSE) +
  labs(x='Value', y='Count') +
  facet_grid(
    score ~ data.type,
    scales='free_x'
  ) +
  theme(aspect.ratio=1.2,
        panel.spacing.x=unit(4.5, "mm"), 
        plot.margin=grid::unit(c(0,0,0,0),"mm"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=9))
hc.perturb.plot

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
hc.graph.plot
ggsave(
  "./report/figures/hc_restart_gaussian_plot.pdf",
  hc.restart.plot,
  width=10,
  height=5,
  device='pdf'
)
ggsave(
  "./report/figures/hc_perturb_gaussian_plot.pdf",
  hc.perturb.plot,
  width=10,
  height=5,
  device='pdf'
)
ggsave(
  "./report/figures/hc_graph_gaussian_plot.pdf",
  hc.graph.plot,
  width=10,
  height=5,
  device='pdf'
)


# Plotting NOTEARS variables
NT.CPDAG.PATH <- sprintf('./data/results/simulated_gaussian/notears_validation/shd_runs_full.csv', sample.number)
NT.MORAL.PATH <- sprintf('./data/results/simulated_gaussian/notears_validation/hd_runs_full.csv', sample.number)
NT.FDR.PATH <- sprintf('./data/results/simulated_gaussian/notears_validation/fdr_runs_full.csv', sample.number)
NT.SENSITIVITY.PATH <- sprintf('./data/results/simulated_gaussian/notears_validation/sensitivity_runs_full.csv', sample.number)
NT.PATHS <- list(
  'SHD (CPDAG)'=NT.CPDAG.PATH, 
  'SHD (Moralised)'=NT.MORAL.PATH,
  'FDR'=NT.FDR.PATH,
  'Sensitivity'=NT.SENSITIVITY.PATH
)
melted.nt.full <- NULL
mu.full <- NULL
for (i in 1:length(NT.PATHS)){
  path <- NT.PATHS[[i]]
  name <- names(NT.PATHS)[i]
  nt.sim.gaussian <- read.csv(path)
  # # Num hyperparams tested
  # num.hyp <- max(table(nt.sim.gaussian$seed))
  # # Count of seeds used
  # t = table(nt.sim.gaussian['seed'])
  # # Analyse seeds which have been assessed by all hyperparam choices
  # nt.sim.gaussian <-nt.sim.gaussian[
  #   is.element(nt.sim.gaussian$seed, names(t[t==5])),
  # ]
  melted.nt <- melt(
    nt.sim.gaussian,
    id.vars=c('l1', 'D'),
    measure.vars=names(nt.sim.gaussian)[length(names(nt.sim.gaussian))]
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
#Filter for extreme values (purely for graphical reasons)
melted.nt.full <- melted.nt.full[
  ((melted.nt.full$value < 5) & (melted.nt.full$data.type =='SHD (CPDAG)'))
  | ((melted.nt.full$value < 15) & (melted.nt.full$data.type =='SHD (Moralised)'))
  | ((melted.nt.full$data.type != 'SHD (CPDAG)') & (melted.nt.full$data.type != 'SHD (Moralised)')),
]

melted.nt.full$data.type <- factor(
  melted.nt.full$data.type,
  levels=c("SHD (CPDAG)", 'SHD (Moralised)', 'FDR', 'Sensitivity')
)
mu.full$data.type <- factor(
  mu.full$data.type,
  levels=c("SHD (CPDAG)", 'SHD (Moralised)', 'FDR', 'Sensitivity')
)
nt.plot <- ggplot(melted.nt.full, aes(x=value, fill=l1)) +
  geom_histogram(position='identity', alpha=0.4, bins=20) + 
  scale_fill_discrete(name = "L1 Regularisation") +
  geom_vline(data=mu.full, aes(xintercept=value, color=l1),
             linetype="dashed", show.legend=FALSE) +
  facet_grid(
    . ~ data.type,
    scales='free_x'
  ) +
  labs(x='Value', y='Count') +
  theme(aspect.ratio=1.2,
        panel.spacing.x=unit(4.5, "mm"), 
        plot.margin=grid::unit(c(0,0,0,0),"mm"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=9))
ggsave(
  "./report/figures/nt_gaussian_plot.pdf",
  nt.plot,
  width=10,
  height=3,
  device='pdf'
)

# Hybrid Plots
HYBRID.CPDAG.PATH <- sprintf('./data/results/simulated_gaussian/hybrid/shd_runs_%s.csv', NUM.RUNS)
HYBRID.MORAL.PATH <- sprintf('./data/results/simulated_gaussian/hybrid/hd_runs_%s.csv', NUM.RUNS)
HYBRID.FDR.PATH <- sprintf('./data/results/simulated_gaussian/hybrid/fdr_runs_%s.csv', NUM.RUNS)
HYBRID.SENSITIVITY.PATH <- sprintf('./data/results/simulated_gaussian/hybrid/sensitivity_runs_%s.csv', NUM.RUNS)
HYBRID.PATHS <- list(
  'SHD (CPDAG)'=HYBRID.CPDAG.PATH, 
  'SHD (Moralised)'=HYBRID.MORAL.PATH,
  'FDR'=HYBRID.FDR.PATH,
  'Sensitivity'=HYBRID.SENSITIVITY.PATH
)

# Plotting histograms
melted.hybrid.full <- NULL
mu.full <- NULL
for (i in 1:length(HYBRID.PATHS)){
  path <- HYBRID.PATHS[[i]]
  name <- names(HYBRID.PATHS)[i]
  hybrid.sim.gaussian <- read.csv(path)
  melted.hybrid <- melt(
    hybrid.sim.gaussian,
    id.vars=c('alpha', 'D'),
    measure.vars=names(hybrid.sim.gaussian)[5:length(hybrid.sim.gaussian)]
  )
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
  levels=c("SHD (CPDAG)", 'SHD (Moralised)', 'FDR', 'Sensitivity')
)
mu.full$data.type <- factor(
  mu.full$data.type,
  levels=c("SHD (CPDAG)", 'SHD (Moralised)', 'FDR', 'Sensitivity')
)
hybrid.plot <- ggplot(melted.hybrid.full, aes(x=value, fill=alpha)) +
  geom_histogram(position='identity', alpha=0.5, bins=20) + 
  scale_fill_discrete(name = "Alpha") +
  geom_vline(data=mu.full, aes(xintercept=value, color=alpha),
             linetype="dashed", show.legend=FALSE) +
  labs(y='Count', x='Value')+
  facet_grid(
    variable ~ data.type,
    scales='free_x'
  ) +
  theme(aspect.ratio=1.2,
        panel.spacing.x=unit(4.5, "mm"), 
        plot.margin=grid::unit(c(0,0,0,0),"mm"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=9))
hybrid.plot
ggsave(
  sprintf("./report/figures/hybrid_gaussian_plot.pdf", NETWORK),
  hybrid.plot,
  width=8,
  height=8,
  device='pdf'
)

# Data Size Plots
DATA.SIZE.SHD <- sprintf('./data/results/simulated_gaussian/data_size/shd_full.csv', NETWORK)
DATA.SIZE.HD <- sprintf('./data/results/simulated_gaussian/data_size/hd_full.csv', NETWORK)
DATA.SIZE.FDR <- sprintf('./data/results/simulated_gaussian/data_size/fdr_full.csv', NETWORK)
DATA.SIZE.SENSITIVITY <- sprintf('./data/results/simulated_gaussian/data_size/sensitivity_full.csv', NETWORK)

DATA.SIZE.PATHS <- list(
  'SHD (CPDAG)'=DATA.SIZE.SHD, 
  'SHD (Moralised)'=DATA.SIZE.HD,
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
    measure.vars=names(data.sim.network)[4:length(names(data.sim.network))]
  )
  if (grepl('HD', name)) melted.data$value <- melted.data$value / melted.data$D 
  melted.data['N'] = lapply(melted.data['N'], as.character)
  melted.data$data.type <- name
  mu <- ddply(melted.data, c('N', 'variable'), summarise, value=mean(value))
  mu$data.type <- name
  
  if (is.null(melted.nt.full)){
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
  levels=c("SHD (CPDAG)", 'SHD (Moralised)', 'FDR', 'Sensitivity')
)
mu.full$data.type <- factor(
  mu.full$data.type,
  levels=c("SHD (CPDAG)", 'SHD (Moralised)', 'FDR', 'Sensitivity')
)
melted.data.size.full$N <- factor(
  melted.data.size.full$N,
  levels=c("50", '100', '250', '500')
)
mu.full$N <- factor(
  mu.full$N,
  levels=c("50", '100', '250', '500')
)

melted.data.size.full <- melted.data.size.full[
  ((melted.data.size.full$value < 100) & (melted.data.size.full$data.type =='SHD (CPDAG)'))
  | ((melted.data.size.full$value < 100) & (melted.data.size.full$data.type =='SHD (Moralised)'))
  | ((melted.data.size.full$data.type != 'SHD (CPDAG)') & (melted.data.size.full$data.type != 'SHD (Moralised)')),
]
mu.full <- mu.full[
  ((mu.full$value < 100) & (mu.full$data.type =='SHD (CPDAG)'))
  | ((mu.full$value < 100) & (mu.full$data.type =='SHD (Moralised)'))
  | ((mu.full$data.type != 'SHD (CPDAG)') & (mu.full$data.type != 'SHD (Moralised)')),
]

data.size.plot <- ggplot(melted.data.size.full, aes(x=value, fill=variable)) +
  geom_histogram(position='identity', alpha=0.4, bins=50) + 
  scale_fill_discrete(name = "Algorithm") +
  geom_vline(data=mu.full, aes(xintercept=value, color=variable),
             linetype="dashed", show.legend=FALSE) +
  facet_grid(
    N ~ data.type,
    scales='free_x'
  ) +
  theme(aspect.ratio=1.2,
        panel.spacing.x=unit(4.5, "mm"), 
        plot.margin=grid::unit(c(0,0,0,0),"mm"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=9))
data.size.plot
ggsave(
  sprintf("./report/figures/data_size_gaussian_plot.pdf", NETWORK),
  data.size.plot,
  width=8,
  height=12,
  device='pdf'
)


# Graph size Plots
GRAPH.SIZE.SHD <- sprintf('./data/results/simulated_gaussian/data_size/shd_full.csv')
GRAPH.SIZE.HD <- sprintf('./data/results/simulated_gaussian/data_size/hd_full.csv')
GRAPH.SIZE.FDR <- sprintf('./data/results/simulated_gaussian/data_size/fdr_full.csv')
GRAPH.SIZE.SENSITIVITY <- sprintf('./data/results/simulated_gaussian/data_size/sensitivity_full.csv')

GRAPH.SIZE.PATHS <- list(
  'SHD (CPDAG)'=GRAPH.SIZE.SHD, 
  'SHD (Moralised)'=GRAPH.SIZE.HD,
  'FDR'=GRAPH.SIZE.FDR,
  'Sensitivity'=GRAPH.SIZE.SENSITIVITY
)
melted.graph.size.full <- NULL
mu.full <- NULL
for (i in 1:length(GRAPH.SIZE.PATHS)){
  path <- GRAPH.SIZE.PATHS[[i]]
  name <- names(GRAPH.SIZE.PATHS)[i]
  data.sim.network <- read.csv(path)
  melted.data <- melt(
    data.sim.network,
    id.vars=c('N', 'D'),
    measure.vars=names(data.sim.network)[4:length(names(data.sim.network))]
  )
  if (grepl('HD', name)) melted.data$value <- melted.data$value / melted.data$D 
  melted.data['D'] = lapply(melted.data['D'], as.character)
  melted.data$data.type <- name
  mu <- ddply(melted.data, c('D', 'variable'), summarise, value=mean(value))
  mu$data.type <- name
  
  if (is.null(melted.graph.size.full)){
    melted.graph.size.full = melted.data
  } else {
    melted.graph.size.full = rbind(melted.graph.size.full, melted.data)
  }
  if (is.null(mu.full)){
    mu.full <- mu
  } else {
    mu.full <- rbind(mu.full, mu)
  }
}
melted.graph.size.full$data.type <- factor(
  melted.graph.size.full$data.type,
  levels=c("SHD (CPDAG)", 'SHD (Moralised)', 'FDR', 'Sensitivity')
)
mu.full$data.type <- factor(
  mu.full$data.type,
  levels=c("SHD (CPDAG)", 'SHD (Moralised)', 'FDR', 'Sensitivity')
)
melted.graph.size.full$D <- factor(
  melted.graph.size.full$D,
  levels=c("25", '50', '100', '150')
)
mu.full$D <- factor(
  mu.full$D,
  levels=c("25", '50', '100', '150')
)

melted.graph.size.full <- melted.graph.size.full[
  ((melted.graph.size.full$value < 1.5) & (melted.graph.size.full$data.type =='SHD (CPDAG)'))
  | ((melted.graph.size.full$value < 2) & (melted.graph.size.full$data.type =='SHD (Moralised)'))
  | ((melted.graph.size.full$data.type != 'SHD (CPDAG)') & (melted.graph.size.full$data.type != 'SHD (Moralised)')),
]
mu.full <- mu.full[
  ((mu.full$value < 100) & (mu.full$data.type =='SHD (CPDAG)'))
  | ((mu.full$value < 100) & (mu.full$data.type =='SHD (Moralised)'))
  | ((mu.full$data.type != 'SHD (CPDAG)') & (mu.full$data.type != 'SHD (Moralised)')),
]

(graph.size.plot <- ggplot(melted.graph.size.full, aes(x=value, fill=variable)) +
  geom_histogram(position='identity', alpha=0.4, bins=50) + 
  scale_fill_discrete(name = "Algorithm", labels=c("GS", "InterIAMB", "H2PC", "MMHC", "SI-HITON-PC" , "NOTEARS", "Hybrid NOTEARS")) +
  geom_vline(data=mu.full, aes(xintercept=value, color=variable),
             linetype="dashed", show.legend=FALSE) +
  facet_grid(
    D ~ data.type,
    scales='free_x'
  ) +
  theme(aspect.ratio=1.2,
        panel.spacing.x=unit(4.5, "mm"), 
        plot.margin=grid::unit(c(0,0,0,0),"mm"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=9)))

ggsave(
  "./report/figures/graph_size_gaussian_plot.pdf",
  graph.size.plot,
  width=8,
  height=7,
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