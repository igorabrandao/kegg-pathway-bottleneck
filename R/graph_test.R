library(igraph)
library(lsa)

df <- data.frame(n1 = c("a", "b", "c"), n2 = c("b", "c", "a"))

g <- graph_from_data_frame(df, FALSE)

plot(g)

V(g)

E(g)

g2 <- delete.vertices(g, "a")

plot(g2)

V(g2)
E(g2)
g3 <- add.vertices(g2, 1, attr=list(name="a"))
V(g3)
E(g3)
plot(g3)
g3 <- add.edges(g3, c("b", "a", "a", "c"))
plot(g3)
V(g3)
E(g3)

V(g3)$AP <- FALSE

V(g3)$AP
V(g3)$name

V(g3)$AP[which(V(g3)$name=="a")] <- TRUE

g = make_graph("Zachary")
c1 = cluster_louvain(g)
length(c1)
membership(c1)
coords = layout_with_fr(g)
# plot the graph
plot(g, layout=coords, vertex.label=NA)
plot(g, vertex.color=membership(c1), layout=coords, vertex.label=NA)
plot(c1, g, layout=coords, vertex.label=NA)
crossing(c1, g)

g2 <- delete.vertices(g, "17")
c2 = cluster_louvain(g2)
coords = layout_with_fr(g2)
plot(c2, g2, layout=coords)
