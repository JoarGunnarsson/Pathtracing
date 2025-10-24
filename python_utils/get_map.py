import argparse
import os
import pathlib

import matplotlib.pyplot as plt
import numpy as np


MAP_FILE_PATH = pathlib.Path(__file__).parent.parent / "maps"


class Modes:
    ALBEDO = "albedo"


def open_file(file_name):
    file_path = MAP_FILE_PATH / file_name
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"{file_path} does not exist")
    return plt.imread(file_path)


def create_albedo_map(input_file, output_file):
    img = open_file(input_file)
    height, width, dim = img.shape[0], img.shape[1], 3
    print(f"Creating new albedo map based on file with width={width} and height={height}.")

    img = np.reshape(img[:, :, 0:3], (width * height * dim,))

    if np.max(img) > 1:
        norm_val = 255
    else:
        norm_val = 1
    img = img / norm_val
    img = np.hstack((np.array([height, width, dim]), img))

    img.tofile(MAP_FILE_PATH / output_file, sep="")


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("in_file", help="input file path, relative to maps/")
    parser.add_argument("out_file", help="output file path, relative to maps/")
    parser.add_argument("-m", "--mode", default=Modes.ALBEDO, help="which mode to use, available modes are: 'albedo'")
    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.mode == Modes.ALBEDO:
        create_albedo_map(args.in_file, args.out_file)
    else:
        raise ValueError(f"{args.mode} is not a valid mode!")


if __name__ == "__main__":
    main()
