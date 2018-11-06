import fractal_dimension as fd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import time
import itertools


def add_vectors(vec1, vec2):
    return [a_i + b_i for a_i, b_i in itertools.zip_longest(vec1, vec2, fillvalue=0)]

def divide(a, b):
    return a/b


def density(mean, number_of_nodes):
    return 2*mean/(number_of_nodes*(number_of_nodes - 1))


def test_density(number_of_nodes, iterations):
    xx = np.arange(10, int((number_of_nodes - 1) * number_of_nodes * 0.5) + 1, 10)
    densities = list(map(density, xx, [number_of_nodes]*len(xx)))
    fractal_dimensions = []
    for i in range(iterations):
        print(i)
        graphs = list(map(nx.fast_gnp_random_graph, [number_of_nodes] * len(densities), densities))
        boxes = list(map(fd.maximum_excluded_mass_burning, graphs))
        log_boxes = list(map(np.log, boxes))
        current_fractal_dimensions = list(map(fd.calculate_fractal_dimension, log_boxes))
        fractal_dimensions = add_vectors(fractal_dimensions, current_fractal_dimensions)
    return_variable = list(map(divide, fractal_dimensions, [iterations]*len(fractal_dimensions)))
    return [0] + return_variable


if __name__ == "__main__":
    node_number = 200
    iterations = 10
    fractal_dimensions = test_density(node_number, iterations)
    means = np.arange(0, int((node_number - 1) * node_number * 0.5) + 10, 10)
    densities = list(map(density, means, [node_number]*(len(means) + 1)))
    plt.plot(densities, fractal_dimensions)
    plt.xlabel("Densities")
    plt.ylabel("Fractal Dimensions")
    plt.savefig("../graphs/fractal_dimensions.png")