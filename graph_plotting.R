library(bnlearn)
require(graph)
require(igraph)
library(purrr)


bn2igraph <- function(g.bn){
  g <- igraph.from.graphNEL(as.graphNEL(g.bn))
}

ecoli <- generate.network.data(N = 500, network.path = 'bnlearn_networks/ecoli70.rds')
ecoli.true<- ecoli$network
ecoli.igraph <- bn2igraph(ecoli.true)
ecoli.df <- as.data.frame(get.edgelist(ecoli.igraph))
ecoli.df.g <- graph.data.frame(ecoli.df, directed=TRUE)
iamb.f <- partial(iamb, test = 'mi-g-sh', undirected = FALSE, alpha=alpha)
inter.g <- inter.iamb.f(ecoli$data)
par(mfrow = c(1, 2))
graphviz.compare(cpdag(bn.net(ecoli$network)), inter.g, shape='ellipse', main = c("True DAG", "IAMB DAG"))
graphviz.compare(moral(cpdag(bn.net(ecoli$network))), moral(inter.g), shape='ellipse', main = c("Moral True DAG", " Moral IAMB DAG"))
print(shd(cpdag(bn.net(ecoli$network)), inter.g))
print(shd(moral(cpdag(bn.net(ecoli$network))), moral(inter.g)))