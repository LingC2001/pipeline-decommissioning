"This file is used to resize the images to a certain size without affecting aspect ratio via padding"

from PIL import Image, ImageOps
import os

def padding(img, expected_size):
    desired_size = expected_size
    delta_width = desired_size - img.size[0]
    delta_height = desired_size - img.size[1]
    pad_width = delta_width // 2
    pad_height = delta_height // 2
    padding = (pad_width, pad_height, delta_width - pad_width, delta_height - pad_height)
    return ImageOps.expand(img, padding)


def resize_with_padding(img, expected_size):
    img.thumbnail((expected_size[0], expected_size[1]))
    # print(img.size)
    delta_width = expected_size[0] - img.size[0]
    delta_height = expected_size[1] - img.size[1]
    pad_width = delta_width // 2
    pad_height = delta_height // 2
    padding = (pad_width, pad_height, delta_width - pad_width, delta_height - pad_height)
    return ImageOps.expand(img, padding)

file_path = "./"
output_path = "resized/"
for img in os.listdir(file_path):
    if img.endswith('.jpg') or img.endswith('.png') or img.endswith('.JPG'):
        image = Image.open(file_path + img)
        image = resize_with_padding(image, (512, 512))
        image.save(output_path+ img)