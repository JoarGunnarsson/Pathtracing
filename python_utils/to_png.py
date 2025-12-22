import argparse
import json
import pathlib

import matplotlib.pyplot as plt
import numpy as np

PROJECT_PATH = pathlib.Path(__file__).parent.parent
IMAGES_FILE_PATH = PROJECT_PATH / "images"
TMP_FILES_PATH = PROJECT_PATH / "temp"


def load_settings_file(file_name):
    with open(file_name) as f:
        return json.load(f)


def load_image_data(file_name, width, height):
    image = np.fromfile(file_name, dtype=np.float64, sep="")
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
    parser.add_argument("--name", type=str, default="result.png", help="name of the resulting image.")
    parser.add_argument(
        "--settings_file",
        type=str,
        default="scenes/settings.json",
        help="the settings file path, relative to main project directory.",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()
    settings = load_settings_file(PROJECT_PATH / args.settings_file)
    WIDTH, HEIGHT = settings.get("WIDTH", 1000), settings.get("HEIGHT", 1000)

    default_image = load_image_data(TMP_FILES_PATH / "raw.dat", WIDTH, HEIGHT)
    plt.imsave(IMAGES_FILE_PATH / args.name, default_image)

    if settings.get("enable_denoising", False):
        denoised_image = load_image_data(TMP_FILES_PATH / "raw_denoised.dat", WIDTH, HEIGHT)
        plt.imsave(IMAGES_FILE_PATH / "denoised" / args.name, denoised_image)


if __name__ == "__main__":
    main()
