import fractal_dimension as fd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import time
import itertools
from scipy import stats
import math


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


def gradients(fractal_dimensions, densities):
    density_gradients = []
    for i in range(len(densities) - 1):
        density_gradients.append((fractal_dimensions[i + 1] - fractal_dimensions[i])
           /(densities[i + 1] - densities[i]))
    return density_gradients


def lower_gradients(gradients, densities):
    mean_gradients = np.mean(gradients)
    lower_gradients = [densities]
    for i in range(len(gradients)):
        if math.abs(gradients[i]) < 0.1*mean_gradients:
            lower_gradients.append(densities[i])
    return lower_gradients


if __name__ == "__main__":
    correlations = []
    plt.figure()
    fractal_dimensions = []
    node_numbers = []
    for i in range(1, 11):
        node_number = 20*i
        iterations = 10
        fractal_dimension = 0
        for j in range(iterations):
            graph = nx.fast_gnp_random_graph(node_number, 0.85)
            boxes = fd.compact_box_burning(graph)
            fractal_dimension += fd.calculate_fractal_dimension(np.log(boxes))
        fractal_dimension = fractal_dimension/iterations
        fractal_dimensions.append(fractal_dimension)
        node_numbers.append(node_number)
    coefficients = np.polyfit(np.log(node_numbers), fractal_dimensions, deg=1)
    print(coefficients)
    plt.plot(np.log(node_numbers), fractal_dimensions)
    plt.show()
        # fractal_dimensions = test_density(node_number, iterations)
        # means = np.arange(0, int((node_number - 1) * node_number * 0.5) + 10, 10)
        # densities = list(map(density, means, [node_number]*(len(means) + 1)))
        # correlations.append(stats.spearmanr(densities, fractal_dimensions)[0])
        # plt.plot(densities, fractal_dimensions)
        # if i == 1:
        #     lower_number_hist_data = lower_gradients(gradients(fractal_dimensions, densities),densities)
        # if i == 5:
        #     mid_number_hist_data = lower_gradients(gradients(fractal_dimensions, densities),densities)
        # if i == 10:
        #     high_number_hist_data = lower_gradients(gradients(fractal_dimensions, densities),densities)
    # plt.xlabel("Densities")
    # plt.ylabel("Fractal Dimensions")
    # plt.savefig("../graphs/fractal_dimensions.png")
    #
    # plt.figure()
    # plt.hist(lower_number_hist_data, 20)
    # plt.xlabel("Densities with little change in the TFD")
    # plt.savefig("../graphs/20_node_histogram.png")
    #
    # plt.figure()
    # plt.hist(mid_number_hist_data, 20)
    # plt.xlabel("Densities with little change in the TFD")
    # plt.savefig("../graphs/10_node_histogram.png")
    #
    # plt.figure()
    # plt.hist(high_number_hist_data, 20)
    # plt.xlabel("Densities with little change in the TFD")
    # plt.savefig("../graphs/200_node_histogram.png")
    #
    # with open('correlations.csv', mode='w') as employee_file:
    #     employee_writer = csv.writer(employee_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    #     employee_writer.writerow(correlations)
    # G = nx.fast_gnp_random_graph(node_number, 0.05)
    # current_time = time.time()
    # print(fd.maximum_excluded_mass_burning(G))
    # print(time.time() - current_time)

    #     for i in range(len(densities) - 1):
    #         density_gradients.append((fractal_dimensions[i + 1] - fractal_dimensions[i])
    #            /(densities[i + 1] - densities[i]))
    # mean_gradient = np.mean(density_gradients)
    # std_gradient = np.std(density_gradients)
    # print(mean_gradient, std_gradient)
    # density_range = []
    # for i in range(len(densities) - 1):
    #     if abs((fractal_dimensions[i + 1] - fractal_dimensions[i])
    #            /(densities[i + 1] -  densities[i])) < 0.1*mean_gradient:
    #         density_range.append(densities[i])
    # print(min(density_range), max(density_range))
