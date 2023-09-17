"""
This script is used to resize and format images to use for the model
"""

import cv2
import torchvision.transforms as T
import numpy as np

def resize_for_model(img):
    """
    resize image into appropriate type for model
    :param img: image to resize
    :returns img: resized and formatted image
    """
    no_transforms = T.Compose([T.ToTensor()])

    # Crop image to square
    crop_size = min(img.shape[0], img.shape[1])
    start_row = (img.shape[0]-crop_size)//2
    end_row = (img.shape[0]-crop_size)//2 + crop_size
    start_col = (img.shape[1]-crop_size)//2
    end_col = (img.shape[1]-crop_size)//2 + crop_size
    img = img[start_row:end_row, start_col:end_col]
    
    # Rescale image
    img = cv2.resize(img, (256,256), interpolation=cv2.INTER_CUBIC)
    img = no_transforms(np.array(img))
    img = np.moveaxis(img.numpy(), 0, 2)

    return img
