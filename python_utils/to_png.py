import matplotlib.pyplot as plt
import numpy as np
import sys

image_name = sys.argv[1]
arr = []
with open("temp/result_data.txt", "r") as file:
    file_contents = file.read()
    for iter, line in enumerate(file_contents.split("\n")):
        if line != "":
            arr.append(line)


image = []

size_and_width = arr[0][5:].split(" ")

size, width = int(size_and_width[0]), int(size_and_width[1])
for line in arr[1:]:
    cleaned_line = line
    rgb_values = cleaned_line.split(" ")
    rgb_values = [int(x) for x in rgb_values if x != " " and x != ""]
    image.append(rgb_values)
 
if len(image) < size * width:
    image[-1] = [0,0,0]
    image.extend([[0, 0, 0]] * (size * width - len(image)))


image = np.array(image, dtype=np.uint8)
image = np.reshape(image, (size, width, 3))
plt.imsave(f"Images/{image_name}", image)
