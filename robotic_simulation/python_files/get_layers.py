from robodk import robolink    # RoboDK API
from robodk import robomath    # Robot toolbox
from robodk.robomath import *

from python_files.layer_segmentation import layer_segmentation
import numpy as np
import cv2

def get_layers(image, cam_height, fov):
    # Getting pipe dimensions from object measurement
    real_frame_width = 2*cam_height*np.tan((fov/2)*np.pi/180)
    image, radii = layer_segmentation(image, real_frame_width)
    radii = [x*1000 for x in radii]
    # diameter = (0.1836+0.168393829)*1000/2
    # length = 960
    #print("length: ",length,"\ndiameter: ",diameter)
    return image, radii