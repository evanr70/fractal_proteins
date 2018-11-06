import fractal_dimension as fd
import matplotlib.pyplot as plt
import numpy as np
import time
import networkx as nx


def test_algorithm(func, number_of_iterations, number_of_nodes):
    return [test_algorithm_n_nodes(func, number_of_iterations, i)
            for i in number_of_nodes]


def test_algorithm_n_nodes(func, number_of_iterations, number_of_nodes):
    total_time = 0
    for i in range(number_of_iterations):
        random_graph = nx.fast_gnp_random_graph(number_of_nodes, 0.5)
        current_time = time.time()
        func(random_graph)
        total_time += (time.time() - current_time)
    return total_time/number_of_iterations



if __name__ == "__main__":
    iter_number = 10
    node_number = np.append(2, np.append(np.arange(10, 90, 10), np.arange(100, 1001, 50)))
    burning_al_times = test_algorithm(fd.compact_box_burning, iter_number, node_number)
    mass_al_times = test_algorithm(fd.maximum_excluded_mass_burning, iter_number, node_number)
    plt.plot(node_number, burning_al_times, label="Box Burning Algorithm")
    plt.plot(node_number, mass_al_times, label="Maximum Mass Algorithm")
    plt.xlabel("Number Of Nodes")
    plt.ylabel("Time Taken (s)")
    plt.title("Measurements for the two different box counting Algorithms")
    plt.legend()
    plt.savefig("../graphs/box_counting.png")