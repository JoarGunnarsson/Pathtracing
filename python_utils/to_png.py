import matplotlib.pyplot as plt
import numpy as np
import sys
import cv2

load_denoising_file = True

def load_image_data(file_name):
    image = cv2.imread(file_name)
    return cv2.cvtColor(image, cv2.COLOR_BGR2RGB)


image_name = sys.argv[1]

default_image = load_image_data("temp/raw_data.PPM")
plt.imsave(f"Images/{image_name}", default_image)

if load_denoising_file:
    denoised_image = load_image_data("temp/denoised_data.PPM")
    plt.imsave(f"Images/denoised/{image_name}", denoised_image)
