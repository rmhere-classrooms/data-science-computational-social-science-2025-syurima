# I Grafy losowe (Erdős-Rényi)

library(igraph)

# Erdos-Renyi graph
# number of nodes 100, rewiring probability 0.05
g <- erdos.renyi.game(p.or.m=0.05, n=100)

# summarize the graph
summary(g) # is the graph weighted? NO

# list all vertices and edges
V(g)
E(g)

# add weights to the edges
E(g)$weight <- runif(length(E(g)), 0.01, 1)

# summarize the graph again
summary(g) # Graph is now weighted (W)

# degree of vertices
degree(g)
# histogram of degrees
hist(degree(g))

# connected components
cl <- clusters(g)
cl
plot(g, vertex.color=cl$membership) # There are 2 clusters (one has only 1 node)

# Visualize graph with vertex size corresponding to PageRank
pr <- page.rank(g)$vector

plot(g, vertex.size=pr*300,vertex.label=NA, edge.arrow.size=.2)
