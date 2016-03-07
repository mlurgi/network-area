

library(rPython)


#### here we load the motif's library
python.exec('from pymfinder import *')
python.load('py_script.py')

load('../networks/network-1')

require(igraph)
                                                                                       
graf <- igraph::graph.adjacency(food_web$M)

a <- walktrap.community(as.undirected(graf))
max(a$modularity)
max(a$membership)

a <- edge.betweenness.community(as.undirected(graf))
max(a$modularity)
max(a$membership)

a <- fastgreedy.community(as.undirected(graf))
max(a$modularity)
max(a$membership)

a <- spinglass.community(graf)
max(a$modularity)
max(a$membership)

write.graph(graf, 'temp.adj', format='edgelist')

file <- 'temp.adj'
python.assign('file', file)
python.exec("m = network_motifs(file)")
motifs <- python.get('m')

python.assign('seed', 345347)
python.assign('iter_factor', 1.0)
python.assign('cooling_factor', 0.993)
python.assign('randoms', 10)
python.exec('mod = modularity_rgraph(file, seed, iter_factor, cooling_factor, randoms)')
modularity <- python.get('mod')




