import numpy as np
import networkx as nx
import glob
import matplotlib.pyplot as plt
import os
import time

plot = False

file_list = glob.glob('data/*.txt')

with open('analysis/time_scale.dat', 'w+') as f:
    print('#filename\tN\tM\tT', file=f)

for file_name in file_list:

    start_time = time.time()

    file_size = os.path.getsize(file_name)

    # if file_size > 5000:
    #     print("skipping: {}\nsize = {}".format(file_name, file_size))
    #     continue

    root_name = file_name[5:-3]

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

    end_time = time.time()

    # if plot:
    #     plt.subplot(111)
    #     nx.draw(G, with_labels=True, font_weight='bold')
    #
    #     plt.show(block=False)
    #     plt.pause(1)
    #     plt.cla()
    #
    # print(G.edges, file=open('outputs' + '\\' + root_name + 'edgeList', 'w+'))
    # print(G.nodes, file=open('outputs' + '\\' + root_name + 'proteinSymbols', 'w+'))

    num_nodes = str(G.number_of_nodes())
    num_edges = str(G.number_of_edges())
    total_time = str(end_time - start_time)

    with open('analysis/time_scale.dat', 'a') as f:
        print("\t".join([root_name, num_nodes, num_edges, total_time]), file=f)

    G.clear()
