import numpy as np
import networkx as nx
import glob
import matplotlib.pyplot as plt
import os
import operator
import time

# file_list = ['data/BIOGRID-ORGANISM-Homo_sapiens-3.5.165.tab2.txt']
# file_list = ['data/BIOGRID-ORGANISM-Rattus_norvegicus-3.5.165.tab2.txt']
# file_list = ['data/BIOGRID-ORGANISM-Human_papillomavirus_16-3.5.165.tab2.txt']
# file_list = ['data/BIOGRID-ORGANISM-Caenorhabditis_elegans-3.5.165.tab2.txt']
file_list = ['data/BIOGRID-ORGANISM-Human_Immunodeficiency_Virus_1-3.5.165.tab2.txt']

stats = []

for file_name in file_list:

    file_size = os.path.getsize(file_name)

    root_name = file_name[6:-3]

    data = np.loadtxt(file_name,
                      dtype=str,
                      delimiter="\t",
                      usecols=(7, 8))

    G = nx.Graph()

    print('\n')

    print(data)

    print('\n')

    print("There are {} edges in this file.".format(data.shape[0]))

    # if not isinstance(data[0], str):
    #     G.add_edges_from(data)
    # else:
    #     G.add_edge(*data)

    for i, line in enumerate(data):
        if i % 10 == 0:
            print("completed {} / {}".format(i, data.shape[0]))
        G.add_edge(*line)

    G.remove_edges_from(G.selfloop_edges())

    # plt.subplot(111)
    # nx.draw(G, with_labels=True, font_weight='bold')
    #
    # plt.show(block=False)
    # plt.pause(1)
    # plt.cla()

    print(G.edges, file=open('outputs' + '\\' + root_name + 'edgeList', 'w+'))
    print(G.nodes, file=open('outputs' + '\\' + root_name + 'proteinSymbols', 'w+'))

    # STATS

    number_of_nodes = G.number_of_nodes()
    number_of_edges = G.number_of_edges()

    best = 3

    start_time = time.time()

    # keep
    print('calculating degree centrality')
    degree_centrality = list(nx.degree_centrality(G).items())
    degree_centrality.sort(key=operator.itemgetter(1), reverse=True)

    new_time = time.time()
    print("time taken = {}".format(new_time - start_time))
    start_time = time.time()

    # keep
    print('calculating eigenvector centrality')
    eigenvector_centrality = list(nx.eigenvector_centrality(G).items())
    eigenvector_centrality.sort(key=operator.itemgetter(1), reverse=True)

    new_time = time.time()
    print("time taken = {}".format(new_time - start_time))
    start_time = time.time()

    # very slow
    print('calculating closeness centrality')
    closeness_centrality = list(nx.closeness_centrality(G).items())
    closeness_centrality.sort(key=operator.itemgetter(1), reverse=True)

    new_time = time.time()
    print("time taken = {}".format(new_time - start_time))
    start_time = time.time()

    # very slow
    print('calculating betweenness centrality')
    betweenness_centrality = list(nx.betweenness_centrality(G).items())
    betweenness_centrality.sort(key=operator.itemgetter(1), reverse=True)

    new_time = time.time()
    print("time taken = {}".format(new_time - start_time))
    start_time = time.time()

    # print('calculating load centrality')
    # load_centrality = list(nx.load_centrality(G).items())
    # load_centrality.sort(key=operator.itemgetter(1), reverse=True)

    # print('calculating harmonic centrality')
    # harmonic_centrality = list(nx.harmonic_centrality(G).items())
    # harmonic_centrality.sort(key=operator.itemgetter(1), reverse=True)
    #
    # print('calculating percolation centrality')
    # percolation_centrality = list(nx.percolation_centrality(G).items())
    # percolation_centrality.sort(key=operator.itemgetter(1), reverse=True)
    #
    # print('calculating second order centrality')
    # second_order_centrality = list(nx.second_order_centrality(G).items())
    # second_order_centrality.sort(key=operator.itemgetter(1), reverse=True)

    centrality_measures = [
        'degree',
        'eigenvector',
        'closeness',
        'betweenness',
        # 'load',
        # 'harmonic',
        # 'percolation',
        # 'second_order'
    ]

    centrality_measures_values = [
        degree_centrality,
        eigenvector_centrality,
        closeness_centrality,
        betweenness_centrality,
        # load_centrality,
        # harmonic_centrality,
        # percolation_centrality,
        # second_order_centrality
    ]

    with open('analysis/{}.dat'.format(root_name), 'w+') as f:
        print('#{}'.format('\t'.join(centrality_measures)))
        for i in range(best):
            row = [measure[i][0] for measure in centrality_measures_values]
            print('{}'.format('\t'.join(row)))

    G.clear()

    print(betweenness_centrality)
