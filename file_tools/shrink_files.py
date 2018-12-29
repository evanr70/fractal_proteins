import pandas as pd
import os


def shrink(input_file, output_file):
    if os.path.isfile(output_file):
        return

    df = pd.read_table(input_file)

    df = df.loc[:, ['Official Symbol Interactor A', 'Official Symbol Interactor B']]

    df = df[df["Official Symbol Interactor A"] != df["Official Symbol Interactor B"]]

    df.to_csv(path_or_buf=output_file, sep="\t", index=False, header=False)


archive_directory = "../biogrid_archive"
archives = map(lambda x: archive_directory + '/' + x, os.listdir(archive_directory))

files = []
for archive in archives:
    files.extend(map(lambda x: archive + '/' + x, os.listdir(archive)))

for file in files:
    save_path = file.replace('tab2.txt', 'csv').replace('biogrid_archive', 'shrunk')
    print(file)
    shrink(file, save_path)
