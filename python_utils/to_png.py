import matplotlib.pyplot as plt
import numpy as np
import sys
from filtering import gaussfft, medfilt


def read_data():
    arr = []
    with open("temp/result_data.txt", "r") as file:
        file_contents = file.read()
        for iter, line in enumerate(file_contents.split("\n")):
            if line != "":
                arr.append(line)
    return arr


def filter_one_image(img, scale=1):
    filtered_channels = []
    for channel in range(3):
        filtered_image = medfilt(img[:, :, channel], scale)
        filtered_channels.append(filtered_image[:, :, None])

    result = np.concatenate(filtered_channels, axis=-1)
    return result


def create_image():
    image_data = []

    data_array = read_data()
    width_and_height = data_array[0][5:].split(" ")

    width, height = int(width_and_height[0]), int(width_and_height[1])

    for line in data_array[1:]:
        cleaned_line = line
        values = cleaned_line.split(" ")
        rgb_values = [float(x) for x in values if x != " " and x != ""]
        image_data.append(rgb_values)
    
    if len(image_data) < height * width:
        image_data[-1] = [0,0,0]
    

    image = np.zeros((width * height, 3))
    image[:len(image_data), :] = np.array(image_data)

    image = np.reshape(image, (height, width, 3))

    print(np.min(image), np.max(image))
    return np.clip(image, 0, 1)


image_name = sys.argv[1]
image = create_image()
plt.imsave(f"Images/{image_name}", image)
