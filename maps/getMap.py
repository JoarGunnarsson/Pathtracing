import numpy as np
import matplotlib.pyplot as plt


img = plt.imread("WorldMap.jpg")


imgData = []

for x in img.flatten():
    imgData.append(str(float(x)/255))


y = img[:, :, 2]

landIdx = np.max(img, axis=-1) > y
img[landIdx] = 100

img[np.logical_not(landIdx)] = 0

plt.imshow(img)
plt.show()
y = img[:, :, 2]
print(y.shape)

reflectionData = []
for x in y.flatten():
    reflectionData.append(str(float(x)/255))
    

print(len(reflectionData))
with open("data.txt", "w") as file:
    file.write(", ".join(imgData))


with open("reflectionData.txt", "w") as file:
    file.write(", ".join(reflectionData))

