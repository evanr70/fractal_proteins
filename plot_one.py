import numpy as np
import networkx as nx
import glob
import matplotlib.pyplot as plt
import os
import time

# file_list = ['data/BIOGRID-ORGANISM-Human_Immunodeficiency_Virus_1-3.5.165.tab2.txt']
# file_list = ['data/BIOGRID-ORGANISM-Human_papillomavirus_16-3.5.165.tab2.txt']
file_list = ['data/BIOGRID-ORGANISM-Human_Herpesvirus_6B-3.5.165.tab2.txt']

for file_name in file_list:

    root_name = file_name[5:-3]

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

    end_time = time.time()

    plt.subplot(111)
    nx.draw(G, with_labels=True, font_weight='bold')

    plt.show()

    G.clear()
