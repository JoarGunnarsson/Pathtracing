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


def tone_map(image):
    return image / (np.max(image, axis=-1) + 1)[:, :, None]


def load_image_data(file_name, width, height):
    image = np.fromfile(file_name, dtype=np.float64, sep="")
    image = np.reshape(image, (height, width, 3))
    image = tone_map(image)

    image_min = np.min(image)
    image_max = np.max(image)
    if image_min < 0:
        print(f"Image minimum is negative ({image_min})!")

    if image_max > 1:
        print(f"Image maximum is greater than 1 ({image_max})!")
    return np.clip(image, 0, 1)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", type=str, default="result.png", help="name of the resulting image.", required=True)
    parser.add_argument(
        "--scene_directory",
        type=str,
        default="scenes/example",
        help="path to the scene directory, relative to the main project directory.",
        required=True,
    )
    return parser.parse_args()


def gamma_correction(image):
    threshold = 0.0031308
    image = np.where(image <= threshold, image * 12.92, 1.055 * image ** (1.0 / 2.4) - 0.055)
    return image


def main():
    args = parse_arguments()
    settings = load_settings_file(PROJECT_PATH / args.scene_directory / "settings.json")
    width, height = settings.get("WIDTH", 1000), settings.get("HEIGHT", 1000)

    use_gamma_correction = settings.get("use_gamma_correction", False)
    result_image = load_image_data(TMP_FILES_PATH / "raw_pixel.dat", width, height)
    if use_gamma_correction:
        result_image = gamma_correction(result_image)

    image_path = IMAGES_FILE_PATH / args.name
    pathlib.Path(image_path.parent).mkdir(exist_ok=True)
    plt.imsave(image_path, result_image)

    if settings.get("enable_atrous_filtering", False) or settings.get("enable_median_filtering", False):
        denoised_image = load_image_data(TMP_FILES_PATH / "raw_denoised.dat", width, height)
        if use_gamma_correction:
            denoised_image = gamma_correction(denoised_image)

        denoised_image_path = IMAGES_FILE_PATH / "denoised" / args.name
        pathlib.Path(denoised_image_path.parent).mkdir(exist_ok=True)
        plt.imsave(denoised_image_path, denoised_image)


if __name__ == "__main__":
    main()
