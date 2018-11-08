import numpy as np
import networkx as nx
import json
import glob
import matplotlib.pyplot as plt
from pprint import pprint
from collections import Counter

if __name__ == "__main__":
    jlogs = glob.glob("jlogs2/*")
    memb = []
    cbb = []
    for jlog in jlogs:
        with open(jlog, "r") as f:
            data = json.load(f)
            for line in data['run']['args']:
                if 'method' in line:
                    method = line.split('=')[1]
                    if method == "memb":
                        memb.append(data)
                    if method == "cbb":
                        cbb.append(data)
    print(memb[0])
    print(cbb[0])
