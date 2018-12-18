import networkx as nx
import graph_magic
import os

def create_nodes_and_edges_from_file(file_name, folder_path=None):
    file_path= file_name.split("/")
    root_name = file_path[-1]
    root_name = root_name[:(len(root_name) - 3)]
    if folder_path:
        root_name = folder_path + "/" + root_name
    print(root_name)

    if folder_path:
        if not os.path.exists("../node_ints/" + folder_path):
            os.makedirs("../node_ints/" + folder_path)

        if not os.path.exists("../edges/" + folder_path):
            os.makedirs("../edges/" + folder_path)
    else:
        if not os.path.exists("../node_ints"):
            os.makedirs("../node_ints")

        if not os.path.exists("../edges"):
            os.makedirs("../edges")

    G = graph_magic.get_graph_from_file(file_name, largest_only=False)

    node_dict = {}

    with open("../node_ints/{}node_dict".format(root_name), "w+") as f:
        for i, node in enumerate(G.nodes):
            node_dict[node] = i
            print("{}\t{}".format(node, i), file=f)

    G = nx.relabel_nodes(G, node_dict)

    with open("../edges/{}tsv".format(root_name), "w+") as f:
        for edge in G.edges:
            print("\t".join(map(str, edge)), file=f)


def create_nodes_and_edges_from_network(file_name, G):
    dir_name_nodes = "../node_ints/" + os.path.dirname(file_name)
    dir_name_edges = "../edges/" + os.path.dirname(file_name)

    if not os.path.exists(dir_name_nodes):
        os.makedirs(dir_name_nodes)

    if not os.path.exists(dir_name_edges):
        os.makedirs(dir_name_edges)

    node_dict = {}

    with open("../node_ints/{}node_dict".format(file_name + "."), "w+") as f:
        for i, node in enumerate(G.nodes):
            node_dict[node] = i
            print("{}\t{}".format(node, i), file=f)

    G = nx.relabel_nodes(G, node_dict)

    with open("../edges/{}tsv".format(file_name + "."), "w+") as f:
        for edge in G.edges:
            print("\t".join(map(str, edge)), file=f)


if __name__ == "__main__":
    # create_nodes_and_edges_from_file("../data/BIOGRID-ORGANISM-Oryza_sativa_Japonica-3.5.165.tab2.txt", "bilogical_networks")
    G = nx.fast_gnp_random_graph(20, 0.5)
    create_nodes_and_edges_from_network("test/example", G)