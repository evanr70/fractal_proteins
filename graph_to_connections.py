import numpy as np
import networkx as nx
import graph_magic
import glob

file_list = glob.glob('data/*.txt')

for file_name in file_list:

    root_name = file_name[5:-3]

    G = graph_magic.get_graph_from_file(file_name)
    with open("connections/" + root_name + "connections", "w+") as f:
        for v in G.nodes:
            print("\t".join([v] + list(G[v]) + ['\t']), file=f)
