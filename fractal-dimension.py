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
            candidate_set = evaluation_set.copy()
            candidate_set.remove(start_node)

            while candidate_set:

                candidate_set = evaluated_candidate_set(evaluation_set, candidate_set, input_sub_graph, length)

                if len(candidate_set) == 1:

                    evaluation_set += candidate_set.pop()

                else:

                    new_target = random.choice(candidate_set)
                    evaluation_set += new_target
                    candidate_set = evaluation_set.copy()
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

    H = graph_magic.get_graph_from_file("data/BIOGRID-ORGANISM-Rattus_norvegicus-3.5.165.tab2.txt")

    # plt.subplot(111)
    # nx.draw(H, with_labels=False, font_weight='bold')
    # plt.show()
    print(compact_box_burning(H))
    node_number = 16
    iteration_number = 20
