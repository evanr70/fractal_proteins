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
from functools import partial
import pandas as pd

import sys

sys.path.insert(0, '../network_tools')
sys.setrecursionlimit(10000)

import fractal_dimension as fd
import get_edges
import read_jlogs

###############################################################################
# Desnity experiments
def experiment_1(node_list, number_of_iterations):
    start_of_experiment = datetime.datetime.now().strftime("%a_%d_%b_%Y_%H:%M:%S")
    start_time = time.time()

    print("Starting to count the boxes: \t\t\t\t\t\t\t" + datetime.datetime.now().strftime("%a, %d %b %Y %H:%M:%S.%f"))

    running_time_stamp = start_time
    for node_number in node_list:
        ave_edges = list(map(round, list(np.arange(1, 100) * 0.01*(node_number - 1)), [3] * 99))
        if os.path.exists("jlog"):
            shutil.rmtree("jlog")

        for ave in ave_edges:
            os.makedirs("../output_data/density_experiments/"
                        + start_of_experiment
                        + "/number_of_vertices_" + str(node_number)
                        + "/number_of_edges_" + str(round(ave*node_number/2, 3)))
            box_number_data = pd.DataFrame()
            tfd_data = pd.DataFrame()
            if time.time() - running_time_stamp > 3600:
                running_time_stamp = time.time()
                print("Experiment has been running for more than " +
                      str(hours_running(start_time, running_time_stamp)) + " hours: \t\t"
                      + datetime.datetime.now().strftime("%a, %d %b %Y %H:%M:%S.%f"))
                print("Currently calculating the number of boxes for a network with " + str(node_number)
                      + "vertices and " + str(ave) + " number of average edges.")

            func = partial(fd.maximum_excluded_mass_burning_erdos, node_number)
            with multiprocessing.Pool() as p:
                p.map(func, [ave]*number_of_iterations)

            jlogs = os.listdir("jlog")
            if len(jlogs) < number_of_iterations:
                print("There are fewer jlogs than expected: \t\t\t\t\t"
                      + datetime.datetime.now().strftime("%a, %d %b %Y %H:%M:%S.%f"))
                raise Exception('too few jlogs')
            if len(jlogs) > number_of_iterations:
                print("There are more jlogs than expected: \t\t\t\t\t"
                      + datetime.datetime.now().strftime("%a, %d %b %Y %H:%M:%S.%f"))
                raise Exception('too many jlogs')
            current_data = read_jlogs.read_jlogs()
            box_number_data = add_box_data_to_DataFrame(current_data, box_number_data, node_number, round(ave*node_number/2, 3))
            tfd_data = add_tfd_data_to_DataFrame(current_data, tfd_data, node_number, round(ave*node_number/2, 3))
            box_number_data.to_csv("../output_data/density_experiments/"
                                   + start_of_experiment + "/number_of_vertices_"
                                   + str(node_number)
                                   + "/number_of_edges_" + str(round(ave*node_number/2, 3)) + "/raw_number_of_boxes.csv",
                                   sep=",")
            tfd_data.to_csv("../output_data/density_experiments/"
                            + start_of_experiment
                            + "/number_of_vertices_"
                            + str(node_number)
                            + "/number_of_edges_" + str(round(ave*node_number/2, 3)) + "/tfd_data.csv",
                            sep=",")
            shutil.rmtree("jlog")


###############################################################################
# Functions for running the experiment
def choose_steps(node_number):
    return list(map(generate_edges_of_specific_densities,
                    list(map(round, list(np.arange(1, 100) * 0.01), [3] * 99)),
                    [node_number]*99))


def hours_running(start_time, current_time):
    return int(((current_time - start_time) - ((current_time - start_time) % 3600)) / 3600)


def add_box_data_to_DataFrame(raw_data, data_frame, node_number, edge):
    for x in raw_data:
        new_data_frame = pd.DataFrame(data={"boxes": x[0], "lengths": [2*r+1 for r in x[1]], "vertices": node_number, "edges": edge})
        data_frame = data_frame.append(new_data_frame)
    return data_frame


def add_tfd_data_to_DataFrame(raw_data, data_frame, node_number, edge):
    tfds = list(fd.calculate_fractal_dimension([2 * r + 1 for r in data[1]], data[0]) for data in raw_data)
    indx = 0
    for tfd in tfds:
        data_frame = data_frame.append(pd.DataFrame(data={"tfd": tfd,
                                  "vertices": node_number,
                                  "density": round(2*edge/(node_number*(node_number - 1)), 3)},
                          index=[indx]))
        indx += 1
    return data_frame

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
        else:
            return correct_edge_number(G, edges, node_number)
    if G.number_of_edges() > edges:
        first_random_node = random.randint(0, node_number - 1)
        connected_nodes = list(nx.neighbors(G, first_random_node))
        if len(connected_nodes) > 0:
            second_random_node = random.choice(connected_nodes)
            G.remove_edge(first_random_node, second_random_node)
            return correct_edge_number(G, edges, node_number)
        else:
            return correct_edge_number(G, edges, node_number)
    else:
        return G


def create_network_of_density(node_number, edges):
    density = (2 * edges) / (node_number * (node_number - 1))
    G = nx.fast_gnp_random_graph(node_number, density)
    return correct_edge_number(G, edges, node_number)


def generate_random_networks(node_number, edges, number_of_iterations):
    node_numbers = [node_number] * number_of_iterations
    edge_numbers = [edges] * number_of_iterations
    return list(map(create_network_of_density, node_numbers, edge_numbers))


def convert_networks_to_files(network_list):
    for network in network_list:
        get_edges.create_nodes_and_edges_from_network(network)


def generate_edges_of_specific_densities(density, nodes):
    return int(round(density * nodes * (nodes - 1) / 2, 0))


def coverting_nodes_to_edges(number_of_iterations, node_number):
    edges = choose_steps(node_number)
    for edge in edges:
        networks = generate_random_networks(node_number, edge, number_of_iterations)
        convert_networks_to_files(networks)


if __name__ == "__main__":
    # try:
    #     experiment_1([100, 200], 50)
    #     os.system('shutdown')
    # except:
    #     os.system('shutdown')
    experiment_1([100, 200], 1)