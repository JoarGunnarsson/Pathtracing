import matplotlib.pyplot as plt
import numpy as np
import sys

def read_data():
    arr = []
    with open("temp/result_data.txt", "r") as file:
        file_contents = file.read()
        for iter, line in enumerate(file_contents.split("\n")):
            if line != "":
                arr.append(line)
    return arr


def create_image():
    image = []

    data_array = read_data()
    size_and_width = data_array[0][5:].split(" ")

    size, width = int(size_and_width[0]), int(size_and_width[1])
    for line in data_array[1:]:
        cleaned_line = line
        rgb_values = cleaned_line.split(" ")
        rgb_values = [float(x) for x in rgb_values if x != " " and x != ""]
        image.append(rgb_values)
    
    if len(image) < size * width:
        image[-1] = [0,0,0]
        image.extend([[0, 0, 0]] * (size * width - len(image)))
    image = np.array(image)
    return np.reshape(image, (size, width, 3))


    
def tone_map(image):
    w, h = image.shape[0], image.shape[1]
    image = np.reshape(image, (w*h, 3))

    max_white = np.max(image, axis=0) + 0.0001
    
    image = (image * (1 + image/max_white ** 2)) / (1 + image)
    return np.reshape(image, (w, h, 3))


image_name = sys.argv[1]
image = create_image()
image = tone_map(image)
plt.imsave(f"Images/{image_name}", image)
