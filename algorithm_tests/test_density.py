import fractal_dimension as fd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


def density(mean, number_of_nodes):
    return 2*mean/(number_of_nodes*(number_of_nodes - 1))


def test_density(number_of_nodes, iterations):
    fractal_dimensions = []
    for j in range(1, int((number_of_nodes - 1)*number_of_nodes*0.5 + 1)):
        total_fd = 0
        for i in range(0, iterations):
            random_ER_graph = nx.fast_gnp_random_graph(number_of_nodes, density(j, number_of_nodes))
            if len(list(nx.connected_component_subgraphs(random_ER_graph))) != number_of_nodes:
                current_boxes = fd.maximum_excluded_mass_burning(random_ER_graph)
                total_fd += fd.calculate_fractal_dimension(np.log(np.array(current_boxes)), len(current_boxes))
            else:
                i -= 1
        fractal_dimensions.append(total_fd/iterations)
    return fractal_dimensions

if __name__ == "__main__":
    node_number = 200
    iterations = 100
    fractal_dimensions = test_density(node_number, iterations)
    means = np.linspace(1, node_number*(node_number - 1), node_number*(node_number - 1))
    densities = list(map(density, means, [node_number]*(len(means) + 1)))
    plt.plot(densities, fractal_dimensions)
    plt.xlabel("Densities")
    plt.ylabel("Fractal Dimensions")
    plt.savefig("../graphs/fractal_dimensions.png")