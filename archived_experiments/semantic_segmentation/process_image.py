import cv2
import os
import numpy as np

input_filepath = "images_to_process"
output_filepath = "images"

for img in os.listdir(input_filepath):
    if img.endswith('.jpg') or img.endswith('.png') or img.endswith('.JPG'):
        # Read image
        image = cv2.imread(input_filepath + "/" + img)
        
        # Crop image to square
        crop_size = min(image.shape[0], image.shape[1])
        start_row = (image.shape[0]-crop_size)//2
        end_row = (image.shape[0]-crop_size)//2 + crop_size
        start_col = (image.shape[1]-crop_size)//2
        end_col = (image.shape[1]-crop_size)//2 + crop_size
        image = image[start_row:end_row, start_col:end_col]
        
        # Rescale image
        image = cv2.resize(image, (512,512), interpolation=cv2.INTER_CUBIC)

        cv2.imwrite(output_filepath+"/"+img, image)