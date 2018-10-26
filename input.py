import numpy as np
import networkx as nx
import glob
import matplotlib.pyplot as plt
import os

file_list = glob.glob('data/*.txt')

for file_name in file_list:

    file_size = os.path.getsize(file_name)

    if file_size > 5000:
        print("skipping: {}\nsize = {}".format(file_name, file_size))
        continue

    root_name = file_name[6:-3]

    data = np.loadtxt(file_name,
                      dtype=str,
                      delimiter="\t",
                      usecols=(7, 8))

    G = nx.Graph()

    print('\n\n')

    print(data)

    print('\n\n')

    if not isinstance(data[0], str):
        G.add_edges_from(data)
    else:
        G.add_edge(*data)

    G.remove_edges_from(G.selfloop_edges())

    plt.subplot(111)
    nx.draw(G, with_labels=True, font_weight='bold')

    plt.show(block=False)
    plt.pause(1)
    plt.cla()

    print(G.edges, file=open('outputs' + '\\' + root_name + 'edgeList', 'w+'))
    print(G.nodes, file=open('outputs' + '\\' + root_name + 'proteinSymbols', 'w+'))

    G.clear()
