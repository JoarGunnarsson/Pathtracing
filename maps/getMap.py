import numpy as np
import matplotlib.pyplot as plt


def convertImageToMap(inputFile, outputFile):
    img = plt.imread(inputFile)
    height, width, dim = img.shape[0], img.shape[1], img.shape[2]
    print(width, height)

    img = np.reshape(img, (width*height, dim))

    img_data = []
    maxValue = np.max(img)
    for rgb in img:
        for x in rgb[:3]:
            img_data.append(str(float(x)/maxValue))

    write_data_to_file(outputFile, img_data, width, height, 3)


def world_map_mask():
    img = np.array(plt.imread("WorldMap.jpg"))
    height, width = img.shape[0], img.shape[1]
    print(width, height)

    y = img[:, :, 2]
    landIdx = np.max(img, axis=-1) > y
    img[landIdx] = 1
    img[np.logical_not(landIdx)] = 0

    plt.imshow(img*255)
    plt.show()

    reflectionData = []
    for x in y.flatten():
        reflectionData.append(str(float(x) / 10))
        
    write_data_to_file("world_roughness.map", reflectionData, width, height, 1)


def write_data_to_file(file_name, img_data, width, height, dimension):
    with open(file_name, "w") as file:
        file.write(str(width) + "\n")
        file.write(str(height) + "\n")
        file.write(str(dimension) + "\n")
        file.write("\n".join(img_data))


convertImageToMap("sandstone_floor.webp", "sandstone_floor.map")
#world_map_mask()
