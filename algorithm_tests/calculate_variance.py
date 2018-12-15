import csv
import sys
import seaborn as sb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def divide(y, x):
    return y/x


def variance(z, y, x):
    return z/x - pow(y, 2)


def process_data(input_data):
    lengths = []
    boxes = []
    means = [0]*len(input_data[0])
    mean_squared = [0]*len(input_data[0])
    iterater = 0
    for x in input_data:
        if iterater < 800:
            for i in range(1, len(x)):
                boxes.append(int(x[i]))
                lengths.append(2*i + 1)
                means[i] += int(x[i])
                mean_squared[i] += pow(int(x[i]), 2)
        iterater += 1
    means = list(map(divide, means, [800]*len(means)))
    mean_squared = list(map(variance, mean_squared, means, [800]*len(means)))
    return lengths, boxes, means, mean_squared


if __name__ == "__main__":
    with open("cbb_boxes.csv", 'r') as f:  #opens PW file
        reader = csv.reader(f)
        data = list(list(rec) for rec in csv.reader(f, delimiter=','))
    processed_data = process_data(data)
    save_data = pd.DataFrame(data={"length": processed_data[0], "boxes": processed_data[1]})

    with open('means_and_variance.csv', mode='w') as employee_file:
        employee_writer = csv.writer(employee_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        employee_writer.writerow(["Means of the number of boxes for each length."])
        employee_writer.writerow(processed_data[2])
        employee_writer.writerow([])

        employee_writer.writerow(["variance of the number of boxes for each length."])
        employee_writer.writerow(processed_data[3])

    # plt.figure()
    # sb.violinplot(y="boxes", x="length", data=save_data)
    # plt.savefig("../graphs/memb_plot.png")