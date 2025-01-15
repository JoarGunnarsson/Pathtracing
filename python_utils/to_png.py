import matplotlib.pyplot as plt
import numpy as np
import sys


def read_data(file_name):
    arr = []
    with open(file_name, "r") as file:
        file_contents = file.read()
        for iter, line in enumerate(file_contents.split("\n")):
            if line != "":
                arr.append(line)
    return arr


def load_image_data(file_name):
    image_data = []

    data_array = read_data(file_name)
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

    image_min = np.min(image)
    image_max = np.max(image)
    if image_min < 0:
        print(f"Image minimum is negative ({image_min})!")

    if image_max > 1:
        print(f"Image maximum is greater than 1 ({image_max})!")

    return np.clip(image, 0, 1)


image_name = sys.argv[1]

denoised_image = load_image_data("temp/raw_data.txt")
plt.imsave(f"Images/raw_{image_name}", denoised_image)

denoised_image = load_image_data("temp/denoised_data.txt")
plt.imsave(f"Images/{image_name}", denoised_image)
