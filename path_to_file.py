from collections import Counter

import networkx as nx
import time
import matplotlib.pyplot as plt

import graph_magic

# G = graph_magic.get_graph_from_file(r"data/BIOGRID-ORGANISM-Caenorhabditis_elegans-3.5.165.tab2.txt")

G = nx.Graph()

# G.add_path(["2", "3", "4", "5", "6", "7", "5"])
# G.add_edge("1", "3")

G.add_path(["9", "7", "8", "5", "4", "3", "2", "1", "3"])
G.add_edge("7", "6")
G.add_edge("7", "5")
G.add_edge("8", "6")
G.add_edge("5", "6")
G.add_edge("4", "1")
G.add_edge("6", "4")

G.remove_edges_from(G.selfloop_edges())

print(print(nx.info(G)))

shortest_paths = nx.all_pairs_shortest_path(G)

sp = nx.all_pairs_shortest_path(G)
steps = []
i = 0

start_time = time.time()
with open(r'C:\Users\Evan\PycharmProjects\fractal_proteins\paths\human.path', 'w+') as f:
    for source in sp:
        if i % 10 == 0:
            print(i)
        for path in list(source[1].items())[:5]:
            if len(path[1]) > 1:
                steps += path[1][1:-1]
                print("\t".join(path[1][1:-1]), file=f)
        i += 1

c = Counter(steps)
n = len(G)

scale = 2 / ((n - 1) * (n - 2))

for count in c.most_common()[:10]:
    print("{}: {}".format(count[0], count[1] * scale))

plt.subplot(111)
nx.draw(G, with_labels=True, font_weight='bold')

plt.show()

print(c.most_common())
