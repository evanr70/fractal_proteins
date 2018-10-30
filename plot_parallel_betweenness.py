"""
====================
Parallel Betweenness
====================

Example of parallel implementation of betweenness centrality using the
multiprocessing module from Python Standard Library.

The function betweenness centrality accepts a bunch of nodes and computes
the contribution of those nodes to the betweenness centrality of the whole
network. Here we divide the network in chunks of nodes and we compute their
contribution to the betweenness centrality of the whole network.

This doesn't work in python2.7.13. It does work in 3.6, 3.5, 3.4, and 3.3.

It may be related to this:
https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map
"""
import glob
from multiprocessing import Pool
import time
import itertools
import os

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx


def chunks(l, n):
    """Divide a list of nodes `l` in `n` chunks"""
    l_c = iter(l)
    while 1:
        x = tuple(itertools.islice(l_c, n))
        if not x:
            return
        yield x


def _betmap(G_normalized_weight_sources_tuple):
    """Pool for multiprocess only accepts functions with one argument.
    This function uses a tuple as its only argument. We use a named tuple for
    python 3 compatibility, and then unpack it when we send it to
    `betweenness_centrality_source`
    """
    return nx.betweenness_centrality_source(*G_normalized_weight_sources_tuple)


def betweenness_centrality_parallel(G, processes=None):
    """Parallel betweenness centrality  function"""
    p = Pool(processes=processes)
    node_divisor = len(p._pool) * 4
    node_chunks = list(chunks(G.nodes(), int(G.order() / node_divisor)))
    num_chunks = len(node_chunks)
    bt_sc = p.map(_betmap,
                  zip([G] * num_chunks,
                      [True] * num_chunks,
                      [None] * num_chunks,
                      node_chunks))

    # Reduce the partial solutions
    bt_c = bt_sc[0]
    for bt in bt_sc[1:]:
        for n in bt:
            bt_c[n] += bt[n]
    return bt_c


if __name__ == "__main__":

    with open('analysis/betweenness_time.dat', 'w+') as f:
        print('#filename\tfilesize\tnodes\tedges\ttime', file=f)

    # file_list = glob.glob('data/*.txt')
    # file_list = ['data/BIOGRID-ORGANISM-Human_papillomavirus_16-3.5.165.tab2.txt']
    file_list = ['data/BIOGRID-ORGANISM-Human_Immunodeficiency_Virus_1-3.5.165.tab2.txt']
    for file_name in file_list:

        root_name = file_name[5:-3]

        file_size = os.path.getsize(file_name)

        # if file_size > 800000:
        #     print("skipping {}".format(file_name))
        #     continue

        data = np.loadtxt(file_name,
                          dtype=str,
                          delimiter="\t",
                          usecols=(7, 8))

        G = nx.Graph()

        if not isinstance(data[0], str):
            G.add_edges_from(data)
        else:
            G.add_edge(*data)

        # if G.number_of_nodes() < 15 or G.number_of_edges() < 15:
        #     continue

        print("")
        print("Computing betweenness centrality for:")
        print(nx.info(G))
        print("\tParallel version")
        start = time.time()
        bt = list(betweenness_centrality_parallel(G).items())
        total_time = time.time() - start
        print("\t\tTime: %.4F" % (time.time() - start))
        # print("\t\tBetweenness centrality for nodes: {}".format(bt))
        print('top three nodes by betweenness:')
        bts = sorted(bt, key=lambda x: x[1], reverse=True)[:3]
        print(bts)

        print('\n\n')
        print(bt)
        with open('analysis/betweenness_time.dat', 'a') as f:
            print("{}\t{}\t{}\t{}\t{}".format(root_name, file_size, G.number_of_nodes(), G.number_of_edges(), total_time), file=f)

