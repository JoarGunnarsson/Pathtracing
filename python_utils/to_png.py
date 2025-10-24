import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np


IMAGES_FILE_PATH = pathlib.Path(__file__).parent.parent / "Images"
TMP_FILES_PATH = pathlib.Path(__file__).parent.parent / "temp"


def load_image_data(file_name, width):
    image = np.fromfile(file_name, dtype=np.float64, sep="")
    height = image.shape[0] // (3 * width)
    image = np.reshape(image, (height, width, 3))

    image_min = np.min(image)
    image_max = np.max(image)
    if image_min < 0:
        print(f"Image minimum is negative ({image_min})!")

    if image_max > 1:
        print(f"Image maximum is greater than 1 ({image_max})!")
    return np.clip(image, 0, 1)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", type=str, help="name of the resulting image.")
    parser.add_argument("--width", type=int, help="the width of the image.")
    return parser.parse_args()


def main():
    args = parse_arguments()
    default_image = load_image_data(TMP_FILES_PATH / "raw.dat", args.width)
    plt.imsave(IMAGES_FILE_PATH / args.name, default_image)

    load_denoising_file = True
    if load_denoising_file:
        denoised_image = load_image_data(TMP_FILES_PATH / "raw_denoised.dat", args.width)
        plt.imsave(IMAGES_FILE_PATH / "denoised" / args.name, denoised_image)


if __name__ == "__main__":
    main()
