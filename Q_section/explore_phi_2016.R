library(ALDEx2)
library(CoDaSeq)
library(zCompositions)
library(igraph)
library(car)
source("http://michael.hahsler.net/SMU/ScientificCompR/code/map.R")

data(ak_op)
data(hmpgenera)

f <- codaSeq.filter(ak_op, min.reads=1000, min.prop=0.005, min.occurrence=0.2, samples.by.row=FALSE)

conds <- c(rep("A", 15), rep("O", 15))

f.x <- aldex.clr(f, conds)

x.phi <- codaSeq.phi(f.x)

phi.cutoff <- .3
x.lo.phi <- subset(x.phi, phi <= phi.cutoff)

# generate a graphical object
g <- graph.data.frame(x.lo.phi, directed=FALSE)

# overview of all the proportional relationships
# this can take a long time!!!
V(g)$label.cex <- 1
OTU.names <- V(g)$name

V(g)$name <- as.character(hmpgenera[V(g)$name, "genus"])


plot(g, layout=layout.fruchterman.reingold.grid(g, weight=0.05/E(g)$phi), vertex.size=1, vertex.color="black")

# http://michael.hahsler.net/SMU/LearnROnYourOwn/code/igraph.html
layout <-layout.fruchterman.reingold(g)
plot(g, layout=layout)

# map centrality
plot(g, layout=layout,
  vertex.size=map(betweenness(g),c(1,15)),
  edge.width=map(edge.betweenness(g), c(1,10)))

# plot page rank (connectedness)
pr <- page.rank(g)$vector
plot(g, layout=layout, vertex.size=map(pr, c(1,20)), edge.arrow.size=.2)

# break into component sub-groups
dg <- decompose.graph(g)
plot(dg[[1]],layout=layout.fruchterman.reingold(dg[[1]]),
  vertex.size=map(page.rank(dg[[1]])$vector, c(1,10)),
  edge.arrow.size=.2)


## communities
# https://users.dimi.uniud.it/~massimo.franceschet/R/communities.html
c1 = cluster_fast_greedy(dg[[1]])
plot(dg[[1]], vertex.color=membership(c1), layout=layout_with_fr(dg[[1]]))

c2 = cluster_leading_eigen(g)
plot(g, vertex.color=membership(c2), layout=layout_with_fr(g), vertex.size=map(pr, c(1,20)), edge.arrow.size=.2)

c3 = cluster_edge_betweenness(g)
plot(g, vertex.color=membership(c3), layout=layout_with_fr(g), vertex.size=map(page.rank(g)$vector, c(1,20)), edge.arrow.size=.2)

plot(g, vertex.color=membership(c3), layout=layout_with_fr(g), vertex.size=map(betweenness(g),c(1,15)),
  edge.width=map(edge.betweenness(g), c(1,10)))




pdf("phi_graph.pdf", height=5, width=5)
plot(dg[[21]])
dev.off()

# # get the clusters from the graph object
 g.clust <- clusters(g)
#
# # data frame containing the names and group memberships of each cluster
 g.df.u <- data.frame(Systematic.name=V(g)$name, cluster=g.clust$membership, cluster.size=g.clust$csize[g.clust$membership], stringsAsFactors=FALSE)
g.df.u <- data.frame(Name=V(g)$name, cluster=g.clust$membership, stringsAsFactors=FALSE)

# generate the list of features and clusters in a data frame
g.df <- g.df.u[order(g.df.u[,"cluster"]),]
