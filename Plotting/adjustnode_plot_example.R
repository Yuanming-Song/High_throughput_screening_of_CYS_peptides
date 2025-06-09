library(network)
library(sna)

# Sample network
g <- rgraph(10, tprob = 0.2)  # random adjacency matrix
net <- network(g, directed = FALSE)

# Fruchterman-Reingold with stiffer springs
layout_coords <- gplot.layout.fruchtermanreingold(
  net,
  layout.par = list(
    area = 0.31,      # smaller area = nodes pushed closer together (tighter springs)
    repulse = .9,   # lower repulsion for tighter clusters
    niter = 9000,      # more iterations for convergence
    arrowhead.cex = 0
  )
)

# Plot using the stiffer spring layout
gplot(net,
      coord = layout_coords,
      vertex.col = "skyblue",
      vertex.cex = 1,
      displaylabels = FALSE,
      arrowhead.cex = 0)
