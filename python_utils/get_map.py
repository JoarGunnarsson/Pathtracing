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
    height, width, dim = img.shape[0], img.shape[1], img.shape[2]
    print(f"Creating new albedo map, based on file with width={width} and height={height}.")

    img = np.reshape(img, (width * height, dim))

    if np.max(img) > 1:
        norm_val = 255
    else:
        norm_val = 1
    img = img / norm_val

    img_data = []
    for rgb in img:
        for x in rgb[:3]:
            img_data.append(str(x))
    print("write")
    write_data_to_file(output_file, img_data, width, height, 3)


def write_data_to_file(file_name, img_data, width, height, dimension):
    with open(file_name, "w") as file:
        file.write(str(width) + "\n")
        file.write(str(height) + "\n")
        file.write(str(dimension) + "\n")
        file.write("\n".join(img_data))


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("in_file", help="input file path, relative to maps/")
    parser.add_argument("out_file", help="output file path, relative to maps/")
    parser.add_argument("-m", "--mode", default=Modes.TEXTURE, help="which mode to use, available modes are: 'albedo'")
    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.mode == Modes.ALBEDO:
        create_albedo_map(args.in_file, args.out_file)
    else:
        raise ValueError(f"{args.mode} is not a valid mode!")


if __name__ == "__main__":
    main()
