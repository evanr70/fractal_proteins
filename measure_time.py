import fractal_dimension as fd
import matplotlib.pyplot as plt
import numpy as np


if __name__ == "__main__":
    node_number = 100
    iter_number = 1000
    burning_al_times = fd.test_algorithm(fd.compact_box_burning, iter_number, node_number)
    mass_al_times = fd.test_algorithm(fd.maximum_excluded_mass_burning, iter_number, node_number)
    plt.plot(np.arange(2, node_number + 1), burning_al_times, label="Box Burning Algorithm")
    plt.plot(np.arange(2, node_number + 1), mass_al_times, label="Maximum Mass Algorithm")
    plt.xlabel("Number Of Nodes")
    plt.ylabel("Time Taken (s)")
    plt.title("Measurements for the two different box counting Algorithms")
    plt.legend()
    plt.savefig("graphs/box_counting.png")