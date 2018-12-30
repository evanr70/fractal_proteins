import network_tools.fractal_dimension as fd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from scipy import stats
import network_tools.graph_magic
import csv


def box_theory(length, gradient, number_of_nodes):
    return number_of_nodes * pow(length, gradient)


def linear_fit(length, gradient, constant):
    return constant + gradient * length


def fit_spearmans(mass_boxes, range):
    ln_box_number = np.log(mass_boxes)
    ln_range = np.log(range)
    gradient = stats.spearmanr(ln_range, ln_box_number)[0] * (
            np.std(ln_box_number) / np.std(ln_range))
    constant = ln_box_number[0] - gradient * ln_range[0]
    return (gradient, constant)


if __name__ == "__main__":
    graph = network_tools.graph_magic.get_graph_from_file(
        "../../BIOGRID-ORGANISM-3.5.165.tab2/BIOGRID-ORGANISM-Human_Immunodeficiency_Virus_1-3.5.165.tab2.txt")
    # graph = graph_magic.get_graph_from_file(
    #     "../../BIOGRID-ORGANISM-3.5.165.tab2/BIOGRID-ORGANISM-Escherichia_coli_K12_MC4100_BW2952-3.5.165.tab2.txt")
    node_number = nx.number_of_nodes(graph)
    print(node_number)
    mass_boxes = fd.compact_box_burning(graph)
    print(mass_boxes)
    # plt.subplot(121)
    # nx.draw(graph, with_labels=True, font_weight='bold')
    # plt.show()
    max_length = len(mass_boxes)
    lengths = np.arange(1, max_length + 1)
    (spearman_gradient, spearman_constant) = fit_spearmans(mass_boxes, lengths)
    (pearson_gradient, pearson_constant) = np.polyfit(np.log(lengths), np.log(mass_boxes), deg=1)
    boxes_spearman_theory = list(map(box_theory, lengths, [spearman_gradient] * max_length,
                                     [node_number] * max_length))
    boxes_pearson_theory = list(map(box_theory,
                                    lengths, [pearson_gradient] * max_length,
                                    [node_number] * max_length))
    ln_boxes_spearman_theory = list(map(linear_fit, np.log(lengths), [spearman_gradient] * max_length,
                                        [spearman_constant] * max_length))
    ln_boxes_pearson_theory = list(map(linear_fit, np.log(lengths), [pearson_gradient] * max_length,
                                       [pearson_constant] * max_length))
    plt.figure()
    plt.plot(lengths, mass_boxes, label="Actual")
    plt.plot(lengths, boxes_spearman_theory, label="Theoretical based on Spearman")
    plt.plot(lengths, boxes_pearson_theory, label="Theoretical based on Pearson")
    plt.xlabel("Lengths of Boxes")
    plt.ylabel("Number of Boxes")
    plt.legend()
    plt.savefig("../graphs/no_of_boxes_from_all_boxes1.png")

    plt.figure()
    plt.scatter(np.log(lengths), np.log(mass_boxes), label="Actual")
    plt.plot(np.log(lengths), ln_boxes_spearman_theory, label="Theoretical based on Spearman")
    plt.plot(np.log(lengths), ln_boxes_pearson_theory, label="Theoretical based on Pearson")
    plt.xlabel("Log Lengths of Boxes")
    plt.ylabel("Log Number of Boxes")
    plt.legend()
    plt.savefig("../graphs/log_no_of_boxes_from_all_boxes1.png")
    print("Difference between Spearman's approximation and actual number of boxes.")
    print(np.subtract(mass_boxes, boxes_spearman_theory))
    print("Difference between Pearson's approximation and actual number of boxes.")
    print(np.subtract(mass_boxes, boxes_pearson_theory))
    print("Difference between Spearman's approximation and actual ln number of boxes.")
    print(np.subtract(np.log(mass_boxes), ln_boxes_spearman_theory))
    print("Difference between Pearson's approximation and actual lnnumber of boxes.")
    print(np.subtract(np.log(mass_boxes), ln_boxes_pearson_theory))

    with open('boxes21.csv', mode='w') as employee_file:
        employee_writer = csv.writer(employee_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        employee_writer.writerow(["Difference between Spearman's approximation and actual number of boxes."])
        employee_writer.writerow(np.subtract(mass_boxes, boxes_spearman_theory))
        employee_writer.writerow(["Mean and Variance of difference"])
        employee_writer.writerow([np.mean(np.subtract(mass_boxes, boxes_spearman_theory)),
                                  np.var(np.subtract(mass_boxes, boxes_spearman_theory))])
        employee_writer.writerow([])

        employee_writer.writerow(["Difference between Pearson's approximation and actual number of boxes."])
        employee_writer.writerow(np.subtract(mass_boxes, boxes_pearson_theory))
        employee_writer.writerow(["Mean and Variance of difference"])
        employee_writer.writerow([np.mean(np.subtract(mass_boxes, boxes_pearson_theory)),
                                  np.var(np.subtract(mass_boxes, boxes_pearson_theory))])
        employee_writer.writerow([])

        employee_writer.writerow(["Difference between Spearman's approximation and actual ln number of boxes."])
        employee_writer.writerow(np.subtract(np.log(mass_boxes), ln_boxes_spearman_theory))
        employee_writer.writerow(["Mean and Variance of difference"])
        employee_writer.writerow([np.mean(np.subtract(np.log(mass_boxes), np.log(boxes_spearman_theory))),
                                  np.var(np.subtract(np.log(mass_boxes), np.log(boxes_spearman_theory)))])
        employee_writer.writerow([])

        employee_writer.writerow(["Difference between Pearson's approximation and actual ln number of boxes."])
        employee_writer.writerow(np.subtract(np.log(mass_boxes), ln_boxes_pearson_theory))
        employee_writer.writerow(["Mean and Variance of difference"])
        employee_writer.writerow([np.mean(np.subtract(np.log(mass_boxes), np.log(boxes_pearson_theory))),
                                  np.var(np.subtract(np.log(mass_boxes), np.log(boxes_pearson_theory)))])
        employee_writer.writerow([])

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
    plt.savefig("../graphs/no_of_boxes1.png")

    plt.figure()
    plt.scatter(np.log(lengths), np.log(mass_boxes), label="Actual")
    plt.plot(np.log(lengths), ln_boxes_spearman_theory, label="Theoretical based on Spearman")
    plt.plot(np.log(lengths), ln_boxes_pearson_theory, label="Theoretical based on Pearson")
    plt.xlabel("Log Lengths of Boxes")
    plt.ylabel("Log Number of Boxes")
    plt.legend()
    plt.savefig("../graphs/log_no_of_boxes1.png")
    print("Difference between Spearman's approximation and actual number of boxes.")
    print(np.subtract(mass_boxes, boxes_spearman_theory))
    print("Difference between Pearson's approximation and actual number of boxes.")
    print(np.subtract(mass_boxes, boxes_pearson_theory))
    print("Difference between Spearman's approximation and actual ln number of boxes.")
    print(np.subtract(np.log(mass_boxes), ln_boxes_spearman_theory))
    print("Difference between Pearson's approximation and actual ln number of boxes.")
    print(np.subtract(np.log(mass_boxes), ln_boxes_pearson_theory))
    with open('boxes22.csv', mode='w') as employee_file:
        employee_writer = csv.writer(employee_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        employee_writer.writerow(["Difference between Spearman's approximation and actual number of boxes."])
        employee_writer.writerow(np.subtract(mass_boxes, boxes_spearman_theory))
        employee_writer.writerow(["Mean and Variance of difference"])
        employee_writer.writerow([np.mean(np.subtract(mass_boxes, boxes_spearman_theory)),
                                  np.var(np.subtract(mass_boxes, boxes_spearman_theory))])
        employee_writer.writerow([])

        employee_writer.writerow(["Difference between Pearson's approximation and actual number of boxes."])
        employee_writer.writerow(np.subtract(mass_boxes, boxes_pearson_theory))
        employee_writer.writerow(["Mean and Variance of difference"])
        employee_writer.writerow([np.mean(np.subtract(mass_boxes, boxes_pearson_theory)),
                                  np.var(np.subtract(mass_boxes, boxes_pearson_theory))])
        employee_writer.writerow([])

        employee_writer.writerow(["Difference between Spearman's approximation and actual ln number of boxes."])
        employee_writer.writerow(np.subtract(np.log(mass_boxes), ln_boxes_spearman_theory))
        employee_writer.writerow(["Mean and Variance of difference"])
        employee_writer.writerow([np.mean(np.subtract(np.log(mass_boxes), np.log(boxes_spearman_theory))),
                                  np.var(np.subtract(np.log(mass_boxes), np.log(boxes_spearman_theory)))])
        employee_writer.writerow([])

        employee_writer.writerow(["Difference between Pearson's approximation and actual ln number of boxes."])
        employee_writer.writerow(np.subtract(np.log(mass_boxes), ln_boxes_pearson_theory))
        employee_writer.writerow(["Mean and Variance of difference"])
        employee_writer.writerow([np.mean(np.subtract(np.log(mass_boxes), np.log(boxes_pearson_theory))),
                                  np.var(np.subtract(np.log(mass_boxes), np.log(boxes_pearson_theory)))])
        employee_writer.writerow([])
