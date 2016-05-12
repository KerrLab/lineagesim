#!/usr/bin/env Rscript

library(igraph)

edge_colors <- c("#ED665D", "#67BF5C")
vertex_shapes <- c("circle", "square")
popsize <- 1e6
fixation_threshold <- 1 - 1e-4

g <- read_graph("tree-end.gml", format="gml")
V(g)$color <- V(g)$abundance > 0
V(g)$size <- 4 + ((V(g)$totalabundance / popsize) * 8) 
V(g)$shape <- vertex_shapes[  as.numeric((V(g)$totalabundance / popsize) > fixation_threshold)   + 1]
E(g)$color <- edge_colors[as.numeric(E(g)$fitnesseffect >= 0) + 1]


#png(filename = "~/TEST1.png", res = 150, width = 6000, height = 6000, type="quartz")
pdf(file = "lineage_tree.pdf", width=30, height=30)
plot.igraph(g, layout=layout_as_tree(g), vertex.label.cex = 1.3)
dev.off()
