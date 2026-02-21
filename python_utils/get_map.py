import argparse
import os
import pathlib

import matplotlib.pyplot as plt
import numpy as np


PROJECT_PATH = pathlib.Path(__file__).parent.parent


class Modes:
    ALBEDO = "albedo"
    TRANSPARENCY = "transparency"
    OPACITY = "opacity"


def open_file(file_name):
    file_path = PROJECT_PATH / file_name
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"{file_path} does not exist")
    return plt.imread(file_path)


def normalize_image(image):
    if np.max(image) > 1:
        norm_val = 255
    else:
        norm_val = 1

    return image / norm_val


def create_albedo_map(input_file, output_file):
    img = open_file(input_file)
    height, width, dim = img.shape[0], img.shape[1], 3
    print(f"Creating new albedo map based on file with width: {width} and height: {height}.")

    img = np.reshape(img[:, :, 0:dim], (width * height * dim,))

    img = normalize_image(img)
    img = np.hstack((np.array([width, height, dim]), img))

    img.tofile(PROJECT_PATH / output_file, sep="")


def create_1D_map(input_file, output_file, invert):
    img = open_file(input_file)
    height, width, dim = img.shape[0], img.shape[1], 1
    print(f"Creating new 1D map based on file with width: {width} and height: {height}.")

    img = np.max(img[:, :, 0:3], axis=-1)
    img = img.flatten()

    img = normalize_image(img)

    if invert:
        img = 1.0 - img

    img = np.hstack((np.array([width, height, dim]), img))

    img.tofile(PROJECT_PATH / output_file, sep="")


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("in_file", help="input file path, relative to the project directory")
    parser.add_argument("out_file", help="output file path, relative to the project directory")
    parser.add_argument(
        "-m",
        "--mode",
        default=Modes.ALBEDO,
        help="which mode to use, available modes are: 'albedo', 'transparency', and 'opacity'",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    if args.mode == Modes.ALBEDO:
        create_albedo_map(args.in_file, args.out_file)
    elif args.mode == Modes.TRANSPARENCY:
        create_1D_map(args.in_file, args.out_file, invert=False)
    elif args.mode == Modes.OPACITY:
        create_1D_map(args.in_file, args.out_file, invert=True)
    else:
        raise ValueError(f"{args.mode} is not a valid mode!")


if __name__ == "__main__":
    main()
