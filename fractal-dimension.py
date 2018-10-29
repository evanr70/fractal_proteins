import networkx as nx
import numpy as np
import math
from scipy import stats
import matplotlib.pyplot as plt
from random import shuffle
import itertools


def colorGraph(inputGraph):
    subGraphs = list(nx.connected_component_subgraphs(inputGraph))
    totalBoxes = colorSubGraph(subGraphs[0])
    for x in range(1, len(subGraphs)):
        totalBoxes = addVectors(totalBoxes, colorSubGraph(subGraphs[x]))
    return calculateFractalDimension(np.log(totalBoxes), len(totalBoxes))


def colorSubGraph(inputSubGraph):
    nodeNames = list(nx.nodes(inputSubGraph))
    lMax = nx.diameter(inputSubGraph) + 1
    numberOfNodes = nx.number_of_nodes(inputSubGraph)
    color = np.arange(1, numberOfNodes*lMax + 1).reshape(numberOfNodes, lMax)
    coloredMatrix = np.zeros((numberOfNodes, lMax))
    if(numberOfNodes < 7):
        iterations = math.factorial(numberOfNodes)
    else:
        iterations = 1000
    for k in range(iterations):
        for i in range(1, numberOfNodes):
            for j in range(i):
                pathLength = len(nx.bidirectional_shortest_path(inputSubGraph, source=nodeNames[j], target=nodeNames[i])) - 1
                for lCurrent in range(lMax):
                    if(pathLength > lCurrent):
                        coloredMatrix[i][lCurrent] = color[j][pathLength - 1]

        currentBoxNumbers = np.zeros(lMax)
        for i in range(lMax):
            currentBoxNumbers[i] = (len(set(coloredMatrix[:, i])))

        if(k < 1):
            boxNumbers = currentBoxNumbers
        else:
            boxNumbers = np.fromiter(map(min, boxNumbers.tolist(), currentBoxNumbers.tolist()), int)

        if(numberOfNodes < 7):
            itertools.permutations(nodeNames)
        else:
            shuffle(nodeNames)
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
    G.add_nodes_from(["a", "b", "c", "d", "e"])
    G.add_edge("a", "b")
    G.add_edge("b", "c")
    G.add_edge("c", "d")
    G.add_edge("d", "e")
    G.add_edge("c", "e")
    plt.subplot(121)
    nx.draw(G, with_labels=True, font_weight='bold')
    plt.show()
    print(colorGraph(G))

    #H = nx.Graph()
    #H.add_nodes_from([1, 2, 3, 4, 5, 6, 7, 8, 9])
    #H.add_edge(1, 2)
    #H.add_edge(2, 3)
    #H.add_edge(3, 4)
    #H.add_edge(4, 5)
    #H.add_edge(3, 5)
    #H.add_edge(6, 7)
    #H.add_edge(8, 9)
    #plt.subplot(121)
    #nx.draw(H, with_labels=True, font_weight='bold')
    #plt.show()
    #print(colorGraph(H))

    #J = nx.Graph()
    #J.add_nodes_from([1, 2, 3, 4, 5, 6, 7, 8, 9])
    #J.add_edge(1, 2)
    #J.add_edge(2, 3)
    #J.add_edge(3, 4)
    #J.add_edge(4, 5)
    #J.add_edge(5, 6)
    #J.add_edge(6, 7)
    #J.add_edge(7, 8)
    #J.add_edge(8, 9)
    #J.add_edge(1, 6)
    #J.add_edge(2, 5)
    #J.add_edge(4, 9)
    #J.add_edge(5, 8)
    #plt.subplot(121)
    #nx.draw(J, with_labels=True, font_weight='bold')
    #plt.show()
    #print(colorSubGraph(J))