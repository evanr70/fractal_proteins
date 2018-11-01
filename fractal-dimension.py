import networkx as nx
import numpy as np
import math
from scipy import stats
import matplotlib.pyplot as plt
import random
import itertools
import time

import graph_magic


# -----------------------------------------------------------------------------------------------------
# These are the set of functions which implement the Compact Box Burning algorithm.

def compact_box_burning(input_graph):
    sub_graphs = list(nx.connected_component_subgraphs(input_graph))
    total_boxes = []

    for sub_graph in sub_graphs:
        max_length = nx.diameter(sub_graph)

        current_boxes = [count_boxes_of_length(sub_graph, i) for i in range(2, max_length + 2)]
        total_boxes = [m + n for m, n in itertools.zip_longest(total_boxes, current_boxes, fillvalue=1)]

    return calculate_fractal_dimension(np.log(total_boxes), len(total_boxes))


def count_boxes_of_length(input_sub_graph, length):
    if length == nx.diameter(input_sub_graph) + 1:
        return 1
    else:
        number_of_boxes = 0
        remaining_nodes = list(nx.nodes(input_sub_graph))

        while remaining_nodes:

            start_node = random.choice(remaining_nodes)
            evaluation_set = [start_node]
            candidate_set = remaining_nodes.copy()
            candidate_set.remove(start_node)

            while candidate_set:

                candidate_set = evaluated_candidate_set(evaluation_set, candidate_set, input_sub_graph, length)
                new_target = random.choice(candidate_set)
                evaluation_set += new_target
                candidate_set.remove(new_target)

            number_of_boxes += 1
            remaining_nodes = remove_elements_from_array(remaining_nodes, evaluation_set)

        return number_of_boxes


def evaluated_candidate_set(evaluation_nodes, candidate_set, sub_graph, length):
    return [x for x in candidate_set if in_union(evaluation_nodes, x, sub_graph, length)]


def in_union(evaluation_nodes, target_node, sub_graph, length):
    for current_node in evaluation_nodes:
        if len(nx.bidirectional_shortest_path(sub_graph, current_node, target_node)) > length:
            return False
    return True


def remove_elements_from_array(first_array, second_array):
    second_set = set(second_array)
    first_set = set(first_array)
    return list(first_set.difference(second_set))


# ------------------------------------------------------------------------------------------------------------
# Implementation of the Maximum-Excluded-Mass-Burning Algorithm
def maximum_excluded_mass_burning_sub_graph(input_sub_graph):
    total_boxes = [nx.number_of_nodes(input_sub_graph)]
    max_path = nx.diameter(input_sub_graph) + 1
    for i in range(2, max_path + 1):
        if(i == max_path):
            total_boxes.append(1)
        if(i % 2 == 0 and i != max_path):
            total_boxes.append(count_boxes_of_length(input_sub_graph, i))
        if(i % 2 == 1 and i != max_path):
            total_boxes.append(count_centers(input_sub_graph, (i - 1)/2))
    return total_boxes


def count_centers(input_sub_graph, radius):
    nx.set_node_attributes(input_sub_graph, False, "covered")
    nx.set_node_attributes(input_sub_graph, False, "center")
    covered = nx.get_node_attributes(input_sub_graph, "covered")
    centers = 0
    while(all(value == True for value in covered.values()) is False):
        (input_sub_graph, max_node) = excluded_masses(input_sub_graph, radius)
        input_sub_graph = cover_nodes(input_sub_graph, radius, max_node)
        covered = nx.get_node_attributes(input_sub_graph, "covered")
        centers += 1
    return centers


def cover_nodes(input_sub_graph, radius, max_node):
    input_sub_graph.node[max_node]["center"] = True
    input_sub_graph.node[max_node]["excluded_mass"] = 0
    paths = nx.single_source_shortest_path(input_sub_graph, max_node, cutoff=radius)
    for path in paths:
        input_sub_graph.node[path]["covered"] = True
    return input_sub_graph


def excluded_masses(input_sub_graph, radius):
    nodes = list(nx.nodes(input_sub_graph))
    max_node_mass = 0
    max_node = 0
    for x in nodes:
        if(input_sub_graph.nodes[x]["center"] is not True):
            input_sub_graph.nodes[x]["excluded_mass"] = excluded_mass(input_sub_graph, radius, x)
            if(input_sub_graph.nodes[x]["excluded_mass"] > max_node_mass):
                max_node = x
                max_node_mass = input_sub_graph.nodes[x]["excluded_mass"]
    return (input_sub_graph, max_node)


def excluded_mass(input_sub_graph, radius, node):
    paths = nx.single_source_shortest_path(input_sub_graph, node, cutoff=radius)
    excluded_mass = 0
    for path in paths:
        if(input_sub_graph.node[path]["covered"] is False and
                input_sub_graph.node[path]["center"] is False):
            excluded_mass += 1
    return excluded_mass


# ------------------------------------------------------------------------------------------------------------
# This is the set of functions which are required by both algorithms.
def calculate_fractal_dimension(ln_box_number, max_length):
    ln_box_length = np.log(np.arange(1, max_length + 1))
    gradient = stats.spearmanr(np.log(np.arange(1, max_length + 1)), ln_box_number)[0] * (
            np.std(ln_box_number) / np.std(ln_box_length))
    return -gradient


# ------------------------------------------------------------------------------------------------------------
# These functions measure the time taken for a given algorithm to run.
def test_algorithm(func, number_of_iterations, number_of_nodes):
    return_times = np.array(test_algorithm_n_nodes(func, number_of_iterations, 2))
    for i in range(3, number_of_nodes + 1):
        return_times = np.append(return_times, test_algorithm_n_nodes(func, number_of_iterations, i))
    return return_times


def test_algorithm_n_nodes(func, number_of_iterations, number_of_nodes):
    total_time = 0
    for i in range(number_of_iterations):
        G = nx.fast_gnp_random_graph(number_of_nodes, 0.5)
        current_time = time.time()
        func(G)
        total_time += (time.time() - current_time)
    return total_time / number_of_iterations


if __name__ == "__main__":
    # H = nx.barabasi_albert_graph(50, 20)

    # H = graph_magic.get_graph_from_file("data/BIOGRID-ORGANISM-Rattus_norvegicus-3.5.165.tab2.txt")
    #
    # # plt.subplot(111)
    # # nx.draw(H, with_labels=False, font_weight='bold')
    # # plt.show()
    # print(compact_box_burning(H))
    node_number = 16
    iteration_number = 20
    G = nx.Graph()
    G.add_nodes_from(["a", "b", "c", "d", "e"])
    G.add_edge("a", "b")
    G.add_edge("b", "c")
    G.add_edge("c", "d")
    G.add_edge("d", "e")
    G.add_edge("c", "e")
    # print(evaluated_candidate_set(["a"], ["b", "c", "d", "e"], G, 2))
    # print(count_boxes_of_length(G, 2))
    # print(count_centers(G, 1))
    print(maximum_excluded_mass_burning_sub_graph(G))
    # plt.subplot(121)
    # nx.draw(G, with_labels=True, font_weight='bold')
    # plt.show()
    # print(colourGraph(G))
    # print(countBoxesOfLength(G, 2, nx.diameter(G)))
    # print(compact_box_burning(G))
    # nx.set_node_attributes(G, False, "covered")
    # nx.set_node_attributes(G, False, "center")
    # covered = nx.get_node_attributes(G, "covered")
    # center = nx.get_node_attributes(G, "center")
    # print(covered)
    # print(all(value == False for value in covered.values()))
    # (G, max_node) = excluded_masses(G, 1)
    # # print(G.node["a"]["excluded_mass"])
    # # print(G.node["b"]["excluded_mass"])
    # # print(G.node["c"]["excluded_mass"])
    # # print(G.node["d"]["excluded_mass"])
    # # print(G.node["e"]["excluded_mass"])
    # # print(max_node)
    # G = cover_nodes(G, 1, "c")
    # covered = nx.get_node_attributes(G, "covered")
    # print(covered)
    # print(all(value == False for value in covered.values()))
    # G.node[max_node]["center"] = True
    # G.node[max_node]["excluded_mass"] = 0
    # G.node["d"]["covered"] = True
    # G.node["e"]["covered"] = True
    # G.node["b"]["covered"] = True



    # I = nx.grid_2d_graph(3, 3)
    # I = excluded_mass(I, 2)
    # for x in nx.nodes(I):
    #     print(x)
    #     print(I.node[x]["excluded_mass"])
    # plt.subplot(121)
    # nx.draw(I, with_labels=True, font_weight='bold')
    # plt.show()

    # H = nx.Graph()
    # H.add_nodes_from([1, 2, 3, 4, 5, 6, 7, 8, 9])
    # H.add_edge(1, 2)
    # H.add_edge(2, 3)
    # H.add_edge(3, 4)
    # H.add_edge(4, 5)
    # H.add_edge(3, 5)
    # H.add_edge(6, 7)
    # H.add_edge(8, 9)
    # plt.subplot(121)
    # nx.draw(H, with_labels=True, font_weight='bold')
    # plt.show()
    # print(compact_box_burning(H))

    # node_number = 16
    # iteration_number = 20
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

    # J = nx.Graph()
    # J.add_nodes_from([1, 2, 3, 4, 5, 6, 7, 8, 9])
    # J.add_edge(1, 2)
    # J.add_edge(2, 3)
    # J.add_edge(3, 4)
    # J.add_edge(4, 5)
    # J.add_edge(5, 6)
    # J.add_edge(6, 7)
    # J.add_edge(7, 8)
    # J.add_edge(8, 9)
    # J.add_edge(1, 6)
    # J.add_edge(2, 5)
    # J.add_edge(4, 9)
    # J.add_edge(5, 8)
    # plt.subplot(121)
    # nx.draw(J, with_labels=True, font_weight='bold')
    # plt.show()
    # print(colourSubGraph(J))
