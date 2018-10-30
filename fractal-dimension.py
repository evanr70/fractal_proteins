import networkx as nx
import numpy as np
import math
from scipy import stats
import matplotlib.pyplot as plt
import random
import itertools
import time


# ---------------------------------------------------------------------------------------------------------
# These sets of functions are used to implement the Greedy Colouring algorithm to count the number of boxes
# required by each length
def colourGraph(inputGraph):
    subGraphs = list(nx.connected_component_subgraphs(inputGraph))
    totalBoxes = colourSubGraph(subGraphs[0])
    for x in range(1, len(subGraphs)):
        totalBoxes = addVectors(totalBoxes, colourSubGraph(subGraphs[x]))
    return calculateFractalDimension(np.log(totalBoxes), len(totalBoxes))


def colourSubGraph(inputSubGraph):
    nodeNames = list(nx.nodes(inputSubGraph))
    lMax = nx.diameter(inputSubGraph) + 1
    numberOfNodes = nx.number_of_nodes(inputSubGraph)
    colour = np.arange(1, numberOfNodes*lMax + 1).reshape(numberOfNodes, lMax)
    colouredMatrix = np.zeros((numberOfNodes, lMax))
    if(numberOfNodes < 8):
        iterations = math.factorial(numberOfNodes)
    else:
        iterations = 10000
    for k in range(iterations):
        for i in range(1, numberOfNodes):
            for j in range(i):
                pathLength = len(nx.bidirectional_shortest_path(inputSubGraph, source=nodeNames[j], target=nodeNames[i])) - 1
                for lCurrent in range(lMax):
                    if(pathLength > lCurrent):
                        colouredMatrix[i][lCurrent] = colour[j][pathLength - 1]

        currentBoxNumbers = np.zeros(lMax)
        for i in range(lMax):
            currentBoxNumbers[i] = (len(set(colouredMatrix[:, i])))

        if(k < 1):
            boxNumbers = currentBoxNumbers
        else:
            boxNumbers = np.fromiter(map(min, boxNumbers.tolist(), currentBoxNumbers.tolist()), int)

        if(numberOfNodes < 8):
            itertools.permutations(nodeNames)
        else:
            random.shuffle(nodeNames)
    return boxNumbers

#-----------------------------------------------------------------------------------------------------
# These are the set of functions which implement the Compact Box Burning algorithm.
def compactBoxBurnig(inputGraph):
    subGraphs = list(nx.connected_component_subgraphs(inputGraph))
    maxLength = nx.diameter(subGraphs[0])
    numberOfNodes = nx.number_of_nodes(subGraphs[0])
    totalBoxes = np.array([numberOfNodes])
    for i in range(2, maxLength + 2):
        totalBoxes = np.append(totalBoxes, countBoxesOfLength(subGraphs[0], i, maxLength))
    for x in range(1, len(subGraphs)):
        maxLength = nx.diameter(subGraphs[x])
        numberOfNodes = nx.number_of_nodes(subGraphs[x])
        currentBoxes = np.array([numberOfNodes])
        for i in range(2, maxLength + 2):
            currentBoxes = np.append(currentBoxes, countBoxesOfLength(subGraphs[0], i, maxLength))
        totalBoxes = addVectors(totalBoxes, currentBoxes)
    return calculateFractalDimension(np.log(totalBoxes), len(totalBoxes))


def countBoxesOfLength(inputSubGraph, length, maxLength):
    if(length == maxLength + 1):
        return 1
    else:
        numberOfBoxes = 0
        remainingNodes = np.array(list(nx.nodes(inputSubGraph)))
        while(len(remainingNodes) > 0):
            startNode = random.choice(remainingNodes)
            evaluationSet = np.array([startNode])
            candidateSet = np.delete(remainingNodes, np.where(remainingNodes == startNode))
            while(len(candidateSet) > 0):
                candidateSet = evaluatedCandidateSet(evaluationSet, candidateSet, inputSubGraph, length)
                if(len(candidateSet) == 1):
                    evaluationSet = np.append(evaluationSet, candidateSet[0])
                    candidateSet = []
                if(len(candidateSet) > 1):
                    newTarget = random.choice(candidateSet)
                    evaluationSet = np.append(evaluationSet, newTarget)
                    candidateSet = np.delete(candidateSet, np.where(candidateSet == newTarget))
            numberOfBoxes += 1
            remainingNodes = removeElementsFromArray(remainingNodes, evaluationSet)
        return numberOfBoxes


def evaluatedCandidateSet(evaluationNodes, candiateSet, subGraph, length):
    return np.array([x for x in candiateSet if inUnion(evaluationNodes, x, subGraph, length)])


def inUnion(evaluationNodes, targetNode, subGraph, length):
    for currentNode in evaluationNodes:
        if(len(nx.bidirectional_shortest_path(subGraph, currentNode, targetNode)) > length):
            return False
    return True


def removeElementsFromArray(firstArray, secondArray):
    secondSet = set(secondArray)
    firstSet = set(firstArray)
    return np.array(list(firstSet.difference(secondSet)))


#------------------------------------------------------------------------------------------------------------
# This is the set of functions which are required by both algorithms.
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


#------------------------------------------------------------------------------------------------------------
# These functions measure the time taken for a given algorithm to run.
def testAlgorithm(func, numberOfIter, numberOfNodes):
    returnTimes = np.array(testAlgorithmNNodes(func, numberOfIter, 2))
    for i in range(3, numberOfNodes + 1):
        returnTimes = np.append(returnTimes, testAlgorithmNNodes(func, numberOfIter, i))
    return returnTimes


def testAlgorithmNNodes(func, numberOfIter, numberOfNodes):
    totalTime = 0
    for i in range(numberOfIter):
        G = nx.fast_gnp_random_graph(numberOfNodes, 0.5)
        currentTime = time.time()
        func(G)
        totalTime += (time.time() -  currentTime)
    return totalTime/numberOfIter


if __name__ == "__main__":
    #G = nx.Graph()
    #G.add_nodes_from(["a", "b", "c", "d", "e"])
    #G.add_edge("a", "b")
    #G.add_edge("b", "c")
    #G.add_edge("c", "d")
    #G.add_edge("d", "e")
    #G.add_edge("c", "e")
    #plt.subplot(121)
    #nx.draw(G, with_labels=True, font_weight='bold')
    #plt.show()
    #print(colourGraph(G))
    #print(countBoxesOfLength(G, 2, nx.diameter(G)))
    #print(compactBoxBurnig(G))

    H = nx.Graph()
    H.add_nodes_from([1, 2, 3, 4, 5, 6, 7, 8, 9])
    H.add_edge(1, 2)
    H.add_edge(2, 3)
    H.add_edge(3, 4)
    H.add_edge(4, 5)
    H.add_edge(3, 5)
    H.add_edge(6, 7)
    H.add_edge(8, 9)
    plt.subplot(121)
    nx.draw(H, with_labels=True, font_weight='bold')
    plt.show()
    print(compactBoxBurnig(H))
    nodeNumber = 16
    iterNumber = 20
    # greedyAlTimes = testAlgorithm(colourGraph, iterNumber, nodeNumber)
    # burningAlTimes = testAlgorithm(compactBoxBurnig, iterNumber, nodeNumber)
    # plt.plot(np.arange(2, nodeNumber + 1), greedyAlTimes)
    # plt.xlabel("Number Of Nodes")
    # plt.ylabel("Time Taken (s)")
    # plt.title("Measurements for Greedy Colouring Algorithm")
    # plt.show()
    #
    # plt.plot(np.arange(2, nodeNumber + 1), burningAlTimes)
    # plt.xlabel("Number Of Nodes")
    # plt.ylabel("Time Taken (s)")
    # plt.title("Measurements for Compact Box Burning Algorithm")
    # plt.show()

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
    #print(colourSubGraph(J))