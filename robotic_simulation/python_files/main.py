from robodk import robolink    # RoboDK API
from robodk import robomath    # Robot toolbox
from robodk.robomath import *

import sys

# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, 'D:\pipe_repo\pipe-recognisation\robotic_simulation')
sys.path.insert(2, "C:\\Users\\Hp\\OneDrive - Monash University\\Documents\\GitHub\\pipe-recognisation\\robotic_simulation")

from python_files.object_measurement import object_measurement
from python_files.layer_segmentation import layer_segmentation
from python_files.operate_side_cam import get_side_cam
from python_files.init_cams import init_cams
from python_files.init_conv_pos import init_conv_pos
from python_files.move_conv import move_conv
from python_files.operate_section_cam import get_section_cam
from python_files.robot_move_dimensions import robot_move_dimensions
from python_files.robot_move_layers import robot_move_layers
from python_files.get_dimensions import get_dimensions
from python_files.get_layers import get_layers
from python_files.init_robot import init_robot

import numpy as np
import cv2
RDK = robolink.Robolink()
SIDE_CAM_HEIGHT = 1.1
SIDE_CAM_FOV = 60
SECTION_CAM_DIS = 0.2
SECTION_CAM_FOV = 60

init_robot()
init_cams()
init_conv_pos()
move_conv()
image_side = get_side_cam()
image_section = get_section_cam()
dimension_image, length, side_diameter = get_dimensions(image_side, SIDE_CAM_HEIGHT, SIDE_CAM_FOV)
layers_image, radii = get_layers(image_section, SECTION_CAM_DIS, SECTION_CAM_FOV)
section_diameter = radii[-1]*2
avg_diameter = (section_diameter + side_diameter)/2


cv2.imshow("dim", dimension_image)
cv2.imshow("layers", layers_image)
cv2.waitKey(0) # PRESS A KEY ON THE IMAGES WINDOW TO CONTINUE
cv2.destroyAllWindows()

robot_move_dimensions(length, avg_diameter)
robot_move_layers(radii)

