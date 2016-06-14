#!/usr/bin/env Rscript

library(igraph)
library(ggplot2)
library(ggplot2bdc)

g <- read_graph("tree-end.gml", format="gml")
vd <- data.frame(vertex.attributes(g))

# pZ <- ggplot(data = vd, aes(x = firstseen, y = fitness, color = fitnessdiff)) +
#     geom_segment(aes(xend = lastseen, yend = fitness)) +
#     scale_color_gradient2(name = "Fitness\nDifference") +
#     labs(x = "Time", y = "Fitness") +
#     theme_bdc_grey()

pZ <- ggplot(data = vd, aes(x = firstseen, y = fitness)) +
    geom_segment(aes(xend = lastseen, yend = fitness), alpha = 0.1) +
    labs(x = "Time", y = "Fitness") +
    theme_bdc_grey()
ggsave_golden("fig_lineage_fitness.pdf", plot = gg_rescale_golden(pZ))
