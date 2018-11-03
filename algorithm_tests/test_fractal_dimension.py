import fractal_dimension as fd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx


def box_theory(length, fractal_dimension, number_of_nodes):
    return number_of_nodes*pow(length, -fractal_dimension)


if __name__ == "__main__":
    node_number = 10000
    mean = 3000
    graph = nx.fast_gnp_random_graph(node_number, (2*mean)/(node_number*(node_number - 1)))
    mass_boxes = fd.maximum_excluded_mass_burning(graph)
    print(mass_boxes)
    # plt.subplot(121)
    # nx.draw(graph, with_labels=True, font_weight='bold')
    # plt.show()
    max_length = len(mass_boxes)
    lengths = np.arange(1, max_length + 1)
    graph_fd_spearman = fd.calculate_fractal_dimension_spearman(np.log(mass_boxes), len(mass_boxes))
    graph_fd = fd.calculate_fractal_dimension(np.log(mass_boxes), len(mass_boxes))
    boxes_spearman_theory = []
    boxes_pearson_theory = []
    print(graph_fd_spearman)
    print(graph_fd)
    for length in lengths:
        boxes_spearman_theory.append(box_theory(length, graph_fd_spearman, node_number))
        boxes_pearson_theory.append(box_theory(length, graph_fd, node_number))
    plt.plot(lengths, mass_boxes, label="Actual")
    plt.plot(lengths, boxes_spearman_theory, label="Theoretical based on Spearman")
    plt.plot(lengths, boxes_pearson_theory, label="Theoretical based on Pearson")
    plt.xlabel("Lengths of Boxes")
    plt.ylabel("Number of Boxes")
    plt.legend()
    plt.savefig("../graphs/no_of_boxes.png")

    plt.plot(np.log(lengths), np.log(mass_boxes), label="Actual")
    plt.plot(np.log(lengths), np.log(boxes_spearman_theory), label="Theoretical based on Spearman")
    plt.plot(np.log(lengths), np.log(boxes_pearson_theory), label="Theoretical based on Pearson")
    plt.xlabel("Log Lengths of Boxes")
    plt.ylabel("Log Number of Boxes")
    plt.legend()
    plt.savefig("../graphs/log_no_of_boxes.png")