import fractal_dimension as fd
import networkx as nx
import graph_magic
import itertools
import numpy as np
import matplotlib.pyplot as plt


def return_variance(total_squared, total, iterations):
    mean = total/iterations
    return total_squared/iterations - pow(mean, 2)


def calculated_variance(func, graph, iterations):
    for i in range(iterations):
        current_boxes = func(graph)
        current_boxes.pop()
        current_boxes.pop(0)
        max_length = len(current_boxes)
        if i == 0:
            total_squared = [0] * max_length
            total = [0] * max_length
        total = [m + n for m, n in itertools.zip_longest(total, current_boxes)]
        total_squared = [m + pow(n, 2) for m, n in itertools.zip_longest(total_squared, current_boxes)]
    return list(map(return_variance, total_squared, total, [iterations]*max_length))

if __name__ == "__main__":
    E_Coli_graph = graph_magic.get_graph_from_file(
        "../../BIOGRID-ORGANISM-3.5.165.tab2/BIOGRID-ORGANISM-Escherichia_coli_K12_MG1655-3.5.165.tab2.txt")
    iterations = 10000
    box_burning_variance = calculated_variance(fd.compact_box_burning, E_Coli_graph, iterations)
    excluded_mass_variance = calculated_variance(fd.compact_box_burning, E_Coli_graph, iterations)
    lengths = np.linspace(2, len(box_burning_variance) + 1, len(box_burning_variance))
    plt.plot(lengths, box_burning_variance, label="Variance of the Box Burning Algorithm")
    plt.plot(lengths, excluded_mass_variance, label="Variance of the Excluded Mass Algorithm")
    plt.xlabel("Lengths of Boxes")
    plt.ylabel("Variance of number of Boxes")
    plt.savefig("../graphs/no_of_boxes_variance.png")