import json
import os


def read_jlogs(file_path="."):
    jlogs = os.listdir(file_path + "/jlog")
    memb = []
    for jlog in jlogs:
        with open(file_path + "/jlog/" + jlog, "r") as f:
            data = json.load(f)
            memb.append((data["size"], data["radius"]))
    return memb


if __name__ == "__main__":
    print(read_jlogs("../algorithm_tests"))
