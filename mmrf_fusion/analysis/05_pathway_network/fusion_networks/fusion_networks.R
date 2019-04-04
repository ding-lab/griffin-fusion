# ==============================================================================
# Fusion network graph
# Steven Foltz (smfoltz@wustl.edu), November 2018
# ==============================================================================
library(ggraph)
library(tidygraph)

graph <- as_tbl_graph(highschool) %>% 
  mutate(Popularity = centrality_degree(mode = 'in'))

ggraph(graph, layout = 'kk') + 
  geom_edge_link2() + 
  geom_node_point(aes(size = Popularity)) + 
  facet_edges(~year) + 
  theme_graph(foreground = 'steelblue', fg_text_colour = 'white')
