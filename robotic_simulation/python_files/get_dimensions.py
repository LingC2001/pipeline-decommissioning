from robodk import robolink    # RoboDK API
from robodk import robomath    # Robot toolbox
from robodk.robomath import *

from python_files.object_measurement import object_measurement
import numpy as np
import cv2

def get_dimensions(image, cam_height, fov):
    # Getting pipe dimensions from object measurement
    real_frame_width = 2*cam_height*np.tan((fov/2)*np.pi/180)
    image, length, diameter = object_measurement(image, real_frame_width, cam_height)
    length = length*1000
    diameter = diameter*1000
    # diameter = (0.1836+0.168393829)*1000/2
    # length = 960
    #print("length: ",length,"\ndiameter: ",diameter)
    return image, length, diameter