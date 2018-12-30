import network_tools.fractal_dimension as fd
import networkx as nx
import network_tools.graph_magic
import itertools
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import pandas as pd
import time


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


def create_data_sets(func, graph, iterations):
    box_numbers = []
    box_lengths = []
    for i in range(iterations):
        current_boxes = func(graph)
        current_boxes.pop()
        current_boxes.pop(0)
        box_numbers.extend(current_boxes)
        box_lengths.extend(np.linspace(2, len(current_boxes) + 1, len(current_boxes)))
    return box_numbers, box_lengths


if __name__ == "__main__":
    E_Coli_graph = network_tools.graph_magic.get_graph_from_file(
        "../../BIOGRID-ORGANISM-3.5.165.tab2/BIOGRID-ORGANISM-Escherichia_coli_K12_MG1655-3.5.165.tab2.txt")
    iterations = 10000
    box_burning = create_data_sets(fd.compact_box_burning, E_Coli_graph, iterations)
    excluded_mass = create_data_sets(fd.maximum_excluded_mass_burning, E_Coli_graph, iterations)
    size_of_boxes = len(fd.compact_box_burning(E_Coli_graph)) - 2
    total_burning = [0]*size_of_boxes
    total_mass = [0]*size_of_boxes
    total_burning_squared = [0]*size_of_boxes
    total_mass_squared = [0]*size_of_boxes
    for i in range(0, len(box_burning[0])):
        total_burning[i % size_of_boxes] += box_burning[0][i % size_of_boxes]
        total_mass[i % size_of_boxes] += excluded_mass[0][i % size_of_boxes]
        total_burning_squared[i % size_of_boxes] += pow(box_burning[0][i % size_of_boxes], 2)
        total_mass_squared[i % size_of_boxes] += pow(excluded_mass[0][i % size_of_boxes], 2)

    for i in range(0, size_of_boxes):
        total_mass[i] = total_mass[i]/iterations
        total_burning[i] = total_burning[i]/iterations
        total_mass_squared[i] = total_mass_squared[i]/iterations
        total_burning_squared[i] = total_burning_squared[i]/iterations

    for i in range(0, size_of_boxes):
        total_burning_squared[i] = total_burning_squared[i] - pow(total_burning[i], 2)
        total_mass_squared[i] = total_mass_squared[i] - pow(total_mass[i], 2)

    save_data = pd.DataFrame(data={"Mean Box Numbers Burning Box": total_burning, "Mean Box Excluded Mass": total_mass,
                                   "Variance Box Numbers Burning Box": total_burning_squared, "Variance Box Numbers Excluded Mass": total_mass_squared})
    print(save_data)
    save_data.to_csv("Al_data.csv", sep='\t')
    data = pd.DataFrame(data={"Number of Boxes": box_burning[0], "Length of Boxes": box_burning[1]})
    sb.violinplot(y="Number of Boxes", x="Length of Boxes", data=data)
    plt.savefig("../graphs/Box_Burning_Violin_plot.png")
    #
    # data = pd.DataFrame(data={"Number of Boxes": excluded_mass[0], "Length of Boxes": excluded_mass[1]})
    # sb.violinplot(y="Number of Boxes", x="Length of Boxes", data=data)
    # plt.savefig("../graphs/Excluded_Mass_Violin_plot.png")

    test_graph = network_tools.graph_magic.get_graph_from_file(
        "../../BIOGRID-ORGANISM-3.5.165.tab2/BIOGRID-ORGANISM-Bos_taurus-3.5.165.tab2.txt")
    print(nx.number_of_nodes(test_graph))
    current_time = time.time()
    print(fd.maximum_excluded_mass_burning(test_graph))
    print(time.time() - current_time)
