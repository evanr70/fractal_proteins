import fractal_dimension as fd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from scipy import stats


def box_theory(length, gradient, number_of_nodes):
    return number_of_nodes*pow(length, gradient)

def linear_fit(length, gradient, constant):
    return constant + gradient*length

def fit_spearmans(mass_boxes, range):
    ln_box_number = np.log(mass_boxes)
    ln_range = np.log(range)
    gradient = stats.spearmanr(ln_range, ln_box_number)[0] * (
            np.std(ln_box_number) / np.std(ln_range))
    constant = ln_box_number[0] - gradient*ln_range[0]
    return (gradient, constant)


if __name__ == "__main__":
    node_number = 10000
    graph = nx.fast_gnp_random_graph(node_number, 0.01)
    mass_boxes = fd.maximum_excluded_mass_burning(graph)
    print(mass_boxes)
    # plt.subplot(121)
    # nx.draw(graph, with_labels=True, font_weight='bold')
    # plt.show()
    max_length = len(mass_boxes)
    lengths = np.arange(1, max_length + 1)
    (spearman_gradient, spearman_constant) = fit_spearmans(mass_boxes, lengths)
    (pearson_gradient, pearson_constant) = np.polyfit(np.log(lengths), np.log(mass_boxes), deg=1)
    boxes_spearman_theory = list(map(box_theory, lengths, [spearman_gradient]*max_length,
                                 [node_number]*max_length))
    boxes_pearson_theory = list(map(box_theory,
                                     lengths, [pearson_gradient]*max_length,
                                [node_number]*max_length))
    ln_boxes_spearman_theory = list(map(linear_fit, np.log(lengths), [spearman_gradient]*max_length,
                                    [spearman_constant]*max_length))
    ln_boxes_pearson_theory = list(map(linear_fit, np.log(lengths), [pearson_gradient] * max_length,
                                    [pearson_constant] * max_length))
    plt.figure()
    plt.plot(lengths, mass_boxes, label="Actual")
    plt.plot(lengths, boxes_spearman_theory, label="Theoretical based on Spearman")
    plt.plot(lengths, boxes_pearson_theory, label="Theoretical based on Pearson")
    plt.xlabel("Lengths of Boxes")
    plt.ylabel("Number of Boxes")
    plt.legend()
    plt.savefig("../graphs/no_of_boxes_from_all_boxes.png")

    plt.figure()
    plt.scatter(np.log(lengths), np.log(mass_boxes), label="Actual")
    plt.plot(np.log(lengths), ln_boxes_spearman_theory, label="Theoretical based on Spearman")
    plt.plot(np.log(lengths), ln_boxes_pearson_theory, label="Theoretical based on Pearson")
    plt.xlabel("Log Lengths of Boxes")
    plt.ylabel("Log Number of Boxes")
    plt.legend()
    plt.savefig("../graphs/log_no_of_boxes_from_all_boxes.png")
    print("Difference between Spearman's approximation and actual number of boxes.")
    print(np.subtract(mass_boxes, boxes_spearman_theory))
    print("Difference between Pearson's approximation and actual number of boxes.")
    print(np.subtract(mass_boxes, boxes_pearson_theory))
    print("Difference between Spearman's approximation and actual ln number of boxes.")
    print(np.subtract(np.log(mass_boxes), ln_boxes_spearman_theory))
    print("Difference between Pearson's approximation and actual lnnumber of boxes.")
    print(np.subtract(np.log(mass_boxes), ln_boxes_pearson_theory))

    mass_boxes.pop(0)
    max_length = len(mass_boxes)
    lengths = np.arange(2, max_length + 2)
    (spearman_gradient, spearman_constant) = fit_spearmans(mass_boxes, lengths)
    (pearson_gradient, pearson_constant) = np.polyfit(np.log(lengths), np.log(mass_boxes), deg=1)
    boxes_spearman_theory = list(map(box_theory, lengths, [spearman_gradient] * max_length,
                                 [node_number] * max_length))
    boxes_pearson_theory = list(map(box_theory,
                                    lengths, [pearson_gradient] * max_length, [node_number] * max_length))
    ln_boxes_spearman_theory = list(map(linear_fit, np.log(lengths),
                                        [spearman_gradient] * max_length, [spearman_constant] * max_length))
    ln_boxes_pearson_theory = list(map(linear_fit, np.log(lengths), [pearson_gradient] * max_length,
                                   [pearson_constant] * max_length))
    plt.figure()
    plt.plot(lengths, mass_boxes, label="Actual")
    plt.plot(lengths, boxes_spearman_theory, label="Theoretical based on Spearman")
    plt.plot(lengths, boxes_pearson_theory, label="Theoretical based on Pearson")
    plt.xlabel("Lengths of Boxes")
    plt.ylabel("Number of Boxes")
    plt.legend()
    plt.savefig("../graphs/no_of_boxes.png")

    plt.figure()
    plt.scatter(np.log(lengths), np.log(mass_boxes), label="Actual")
    plt.plot(np.log(lengths), ln_boxes_spearman_theory, label="Theoretical based on Spearman")
    plt.plot(np.log(lengths), ln_boxes_pearson_theory, label="Theoretical based on Pearson")
    plt.xlabel("Log Lengths of Boxes")
    plt.ylabel("Log Number of Boxes")
    plt.legend()
    plt.savefig("../graphs/log_no_of_boxes.png")
    print("Difference between Spearman's approximation and actual number of boxes.")
    print(np.subtract(mass_boxes, boxes_spearman_theory))
    print("Difference between Pearson's approximation and actual number of boxes.")
    print(np.subtract(mass_boxes, boxes_pearson_theory))
    print("Difference between Spearman's approximation and actual ln number of boxes.")
    print(np.subtract(np.log(mass_boxes), ln_boxes_spearman_theory))
    print("Difference between Pearson's approximation and actual lnnumber of boxes.")
    print(np.subtract(np.log(mass_boxes), ln_boxes_pearson_theory))