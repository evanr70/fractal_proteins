import networkx as nx
import numpy as np
from scipy import stats


def colorMatrix(inputGraph):
    lMax = nx.diameter(inputGraph)
    numberOfNodes = nx.number_of_nodes(inputGraph)
    color = np.arange(1, numberOfNodes*lMax + 1).reshape(numberOfNodes, lMax)
    coloredMatrix = np.zeros((numberOfNodes, lMax))
    for i in range(1, numberOfNodes):
        for j in range(i):
            pathLength = len(nx.bidirectional_shortest_path(G, j + 1, i + 1)) - 1
            for lCurrent in range(lMax):
                if(pathLength > lCurrent + 1):
                    coloredMatrix[i][lCurrent] = color[j][pathLength - 1]


    lnBoxNumbers = np.zeros(lMax)
    for i in range(lMax):
        lnBoxNumbers[i] = np.log((len(set(coloredMatrix[:, i]))))

    lnBoxLength = np.log(np.arange(1, 4))
    gradient = stats.spearmanr(np.log(np.arange(1, 4)), lnBoxNumbers)[0]*(np.std(lnBoxNumbers)/np.std(lnBoxLength))
    return -gradient


if __name__ == "__main__":
    G = nx.Graph()
    G.add_nodes_from([1, 2, 3, 4, 5])
    G.add_edge(1, 2)
    G.add_edge(2, 3)
    G.add_edge(3, 4)
    G.add_edge(4, 5)
    G.add_edge(3, 5)
    #plt.subplot(121)
    #nx.draw(G, with_labels=True, font_weight='bold')
    #plt.show()
    print(colorMatrix(G))