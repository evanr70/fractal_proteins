import networkx as nx
import graph_magic

import glob

def create_nodes_and_edges(file_name):
    root_name = file_name[5:-3]

    G = graph_magic.get_graph_from_file(file_name, largest_only=False)

    node_dict = {}

    with open("../node_ints/{}node_dict".format(root_name), "w+") as f:
        for i, node in enumerate(G.nodes):
            node_dict[node] = i
            print("{}\t{}".format(node, i), file=f)

    G = nx.relabel_nodes(G, node_dict)

    with open("../edges/{}tsv".format(root_name), "w+") as f:
        for edge in G.edges:
            print("\t".join(map(str, edge)), file=f)
