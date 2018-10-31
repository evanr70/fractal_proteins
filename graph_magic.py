import numpy as np
import networkx as nx


def get_graph_from_file(file_name):

    data = np.loadtxt(file_name,
                      dtype=str,
                      delimiter="\t",
                      usecols=(7, 8))

    G = nx.Graph()

    if not isinstance(data[0], str):
        G.add_edges_from(data)
    else:
        G.add_edge(*data)

    G.remove_edges_from(G.selfloop_edges())
    return G
