import collections

import numpy as np
import networkx as nx
import glob
import matplotlib.pyplot as plt
import os
import operator
import time

# file_name = 'data/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.5.165.tab2.txt'
# file_name = 'data/BIOGRID-ORGANISM-Homo_sapiens-3.5.165.tab2.txt'
# file_list = ['data/BIOGRID-ORGANISM-Rattus_norvegicus-3.5.165.tab2.txt']
# file_list = ['data/BIOGRID-ORGANISM-Caenorhabditis_elegans-3.5.165.tab2.txt']

file_list = glob.glob('data/*.txt')

fig = plt.figure(figsize=(1920, 949))

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

    # plt.subplot(111)
    # nx.draw(G, with_labels=True, font_weight='bold')
    #
    # plt.show(block=False)
    # plt.pause(1)
    # plt.cla()

    # STATS

    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
    degree_count = collections.Counter(degree_sequence)

    degree_count = {x: degree_count[x] for x in degree_count if degree_count[x] >= 1}

    deg, cnt = zip(*degree_count.items())

    print(degree_count)

    plt.bar(deg, cnt, width=0.80, color='b')

    plt.title("Degree histogram")
    plt.ylabel("Count")
    plt.xlabel("Degree")
    plt.xlim(0, 550)

    plt.show(block=False)
    plt.pause(0.1)

    plt.savefig('analysis/deg_dist/' + root_name + 'png')

    G.clear()
    plt.cla()
