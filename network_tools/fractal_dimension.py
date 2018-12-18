import numpy as np
import subprocess
import os

# ------------------------------------------------------------------------------------------------------------
# Implementation of the Maximum-Excluded-Mass-Burning Algorithm
def maximum_excluded_mass_burning(file_name):
    command = ["../graph-sketch-fractality/bin/box_cover", "-type=tsv", "-method=memb", "-graph=" + file_name,
               " -seed=47398"]
    FNULL = open(os.devnull, 'w')
    subprocess.call(command, stdout=FNULL, stderr=FNULL)


# ------------------------------------------------------------------------------------------------------------
# This is the set of functions which are required by both algorithms.
def calculate_fractal_dimension(box_sizes, number_of_boxes):
    return -np.polyfit(np.log(box_sizes), np.log(number_of_boxes), deg=1)[0]


if __name__ == "__main__":
    maximum_excluded_mass_burning("../edges/example.tsv")
    # subprocess.call(["ls", "-l"])