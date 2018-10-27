import networkx as nx
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


def colorGraph(inputGraph):
    subGraphs = list(nx.connected_component_subgraphs(inputGraph))
    firstSubGraph = True
    for x in subGraphs:
        if(firstSubGraph):
            firstSubGraphTotals = colorSubGraph(x)
            totalBoxes = addVectors(np.zeros(len(firstSubGraphTotals)), firstSubGraphTotals)
            firstSubGraph = False
            print(totalBoxes)
        else:
            totalBoxes = addVectors(totalBoxes, colorSubGraph(x))
    #return calculateFractalDimension(np.log(totalBoxes), len(totalBoxes))
    return totalBoxes

def colorSubGraph(inputSubGraph):
    lMax = nx.diameter(inputSubGraph)
    numberOfNodes = nx.number_of_nodes(inputSubGraph)
    color = np.arange(1, numberOfNodes*lMax + 1).reshape(numberOfNodes, lMax)
    coloredMatrix = np.zeros((numberOfNodes, lMax))
    for i in range(1, numberOfNodes):
        for j in range(i):
            pathLength = len(nx.bidirectional_shortest_path(G, j + 1, i + 1)) - 1
            for lCurrent in range(lMax):
                if(pathLength > lCurrent + 1):
                    coloredMatrix[i][lCurrent] = color[j][pathLength - 1]

    boxNumbers = np.zeros(lMax)
    for i in range(lMax):
        boxNumbers[i] = (len(set(coloredMatrix[:, i])))

    return boxNumbers


def calculateFractalDimension(lnBoxNumber, maxLength):
    lnBoxLength = np.log(np.arange(1, maxLength + 1))
    gradient = stats.spearmanr(np.log(np.arange(1, maxLength + 1)), lnBoxNumber)[0] * (np.std(lnBoxNumber) / np.std(lnBoxLength))
    return -gradient

def addVectors(vec1, vec2):
    if(len(vec1) > len(vec2)):
        vec2 = np.append(vec2, np.array([1]*(len(vec1) - len(vec2))))
        return np.add(vec1, vec2)
    if(len(vec2) > len(vec1)):
        return addVectors(vec2, vec1)
    else:
        return np.add(vec1, vec2)


if __name__ == "__main__":
    G = nx.Graph()
    G.add_nodes_from([1, 2, 3, 4, 5, 6, 7])
    G.add_edge(1, 2)
    G.add_edge(2, 3)
    G.add_edge(3, 4)
    G.add_edge(4, 5)
    G.add_edge(3, 5)
    G.add_edge(6, 7)
    plt.subplot(121)
    graphs = list(nx.connected_component_subgraphs(G))
    print(len(graphs))
    nx.draw(G, with_labels=True, font_weight='bold')
    plt.show()
    print(colorGraph(G))