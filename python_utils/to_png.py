import argparse

import matplotlib.pyplot as plt
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("--name", type=str, help="Name of the resulting image.")
parser.add_argument("--width", type=int, help="The width of the image.")
args = parser.parse_args()


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


default_image = load_image_data("temp/raw.dat", args.width)
plt.imsave(f"Images/{args.name}", default_image)

load_denoising_file = True
if load_denoising_file:
    denoised_image = load_image_data("temp/raw_denoised.dat", args.width)
    plt.imsave(f"Images/denoised/{args.name}", denoised_image)
