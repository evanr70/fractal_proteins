import matplotlib.pyplot as plt
import networkx as nx

import graph_magic

G = graph_magic.get_graph_from_file("data/BIOGRID-ORGANISM-Bos_taurus-3.5.165.tab2.txt")

plt.subplot(111)
nx.draw(G, with_labels=True, font_weight='bold')

plt.show()

G.clear()
