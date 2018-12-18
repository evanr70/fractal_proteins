import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import datetime
import time
import itertools
from scipy import stats
import math
import random
import shutil
import os
import multiprocessing

import sys
sys.path.insert(0, '../network_tools')

import fractal_dimension as fd
import get_edges

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
###############################################################################
# Desnity experiments
def experiment_1(node_list, number_of_iterations):
    start_of_experiment = datetime.datetime.now().strftime("%a_%d_%b_%Y_%H:%M:%S")
    start_time = time.time()

    print("Started experiment 1: \t\t\t\t\t\t\t\t\t" + datetime.datetime.now().strftime("%a, %d %b %Y %H:%M:%S.%f"))
    os.makedirs("../output_data/density_experiments/" + start_of_experiment + "/boxes")
    if os.path.exists("../edges"):
        shutil.rmtree("../edges")
    if os.path.exists("../node_ints"):
        shutil.rmtree("../node_ints")

    print("Creating edges for test: \t\t\t\t\t\t\t\t" + datetime.datetime.now().strftime("%a, %d %b %Y %H:%M:%S.%f"))
    for node_number in node_list:
        edges = choose_steps(node_number)
        for edge in edges:
            networks = generate_random_networks(node_number, edge, number_of_iterations)
            convert_networks_to_files(networks, node_number, edge)

    print("Finished creating edges files: \t\t\t\t\t\t\t" + datetime.datetime.now().strftime("%a, %d %b %Y %H:%M:%S.%f"))
    print("Starting to count the boxes: \t\t\t\t\t\t\t" + datetime.datetime.now().strftime("%a, %d %b %Y %H:%M:%S.%f"))
    running_time_stamp = start_time
    for node_number in node_list:
        os.makedirs("../output_data/density_experiments/" + start_of_experiment + "/boxes/" + "networks_of_size_"
                    + str(node_number))
        edges = os.listdir("../edges/" + str(node_number))
        for edge in edges:
            if time.time() - running_time_stamp > 3600:
                running_time_stamp = time.time()
                print("Experiment has been running for more than " +
                      str(hours_running(start_time, running_time_stamp)) + " hours: \t\t"
                      + datetime.datetime.now().strftime("%a, %d %b %Y %H:%M:%S.%f"))
            networks = os.listdir("../edges/" + str(node_number) + "/" + str(edge))
            networks = list(["../edges/" + str(node_number) + "/" + str(edge) + "/" +
                             network for network in networks])
            with multiprocessing.Pool() as p:
                p.map(fd.maximum_excluded_mass_burning, networks)
            jlogs = os.listdir("jlog")





###############################################################################
# Functions for running the experiment
def choose_steps(node_number):
    if (node_number * (node_number - 1))/400 < 10:
        steps = int(math.floor((node_number * (node_number - 1))/200))
        return list(np.arange(steps, int((node_number - 1) * node_number * 0.5), steps))
    else:
        steps = int(math.floor((node_number * (node_number - 1)) / 200))
        return list([10, 20, 30, 40, 50, 60, 70, 80, 90]) + \
               list(np.arange(100, int((node_number - 1) * node_number * 0.5), steps))


def hours_running(start_time, current_time):
    return int(((current_time - start_time) - ((current_time - start_time) % 3600))/3600)


###############################################################################
# Functions to set up the tests
def correct_edge_number(G, edges, node_number):
    if G.number_of_edges() < edges:
        first_random_node = random.randint(0, node_number - 1)
        non_connected_nodes = list(nx.non_neighbors(G, first_random_node))
        if len(non_connected_nodes) > 0:
            second_random_node = random.choice(non_connected_nodes)
            nx.add_path(G, [first_random_node, second_random_node])
            return correct_edge_number(G, edges, node_number)
        else: return correct_edge_number(G, edges, node_number)
    if G.number_of_edges() > edges:
        first_random_node = random.randint(0, node_number - 1)
        connected_nodes = list(nx.neighbors(G, first_random_node))
        if len(connected_nodes) > 0:
            second_random_node = random.choice(connected_nodes)
            G.remove_edge(first_random_node, second_random_node)
            return correct_edge_number(G, edges, node_number)
        else: return correct_edge_number(G, edges, node_number)
    else: return G


def create_network_of_density(node_number, edges):
    density = (2*edges)/(node_number*(node_number - 1))
    G = nx.fast_gnp_random_graph(node_number, density)
    return correct_edge_number(G, edges, node_number)


def generate_random_networks(node_number, edges, number_of_iterations):
    node_numbers = [node_number]*number_of_iterations
    edge_numbers = [edges]*number_of_iterations
    return list(map(create_network_of_density, node_numbers, edge_numbers))


def convert_networks_to_files(network_list, node_number, edge_number):
    for network in network_list:
        get_edges.create_nodes_and_edges_from_network(str(node_number) + "/"
                                                      + str(edge_number) + "/" + str(time.time()), network)


if __name__ == "__main__":
    experiment_1([50], 10)
    # print(datetime.datetime.now().strftime("%a, %d %b %Y %H:%M:%S.%f"))
    # G = nx.fast_gnp_random_graph(50, 0.1)
    # plt.subplot(111)
    # nx.draw(G, with_labels=True, font_weight='bold')
    #
    # plt.show()
    #
    # convert_networks_to_files(list([G]), 50, nx.number_of_edges(G))
    # G = nx.Graph()
    # G.add_nodes_from([0,1,2,3])
    # G.add_path([0,1,2,3])
    # G = correct_edge_number(G, 2, 4)
    # print(type(G))
    # plt.subplot(111)
    # nx.draw(G, with_labels=True, font_weight='bold')
    #
    # plt.show()
