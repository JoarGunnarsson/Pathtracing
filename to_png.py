import sys
import matplotlib.pyplot as plt
import numpy as np
arr = []
for line in sys.stdin:
    arr.append(line)


image = []
size_and_width = arr[1][:-1].split(" ")
size, width = int(size_and_width[0]), int(size_and_width[1])
for line in arr[3:]:
    cleaned_line = line[:-1]
    rgb_values = cleaned_line.split(" ")
    rgb_values = [int(x) for x in rgb_values if x != " " and x != ""]
    image.append(rgb_values)

image = np.array(image, dtype=np.uint8)
image = np.reshape(image, (size, width, 3))
plt.imsave("Images/result.png", image)

