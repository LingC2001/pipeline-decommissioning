from robodk import robolink    # RoboDK API
from robodk import robomath    # Robot toolbox
from robodk.robomath import *

import numpy as np
import cv2
RDK = robolink.Robolink()

def robot_move_layers(layers):
    # Defining robot
    ROBOT_NAME = "KUKA LBR iiwa 14 R820"
    robot = RDK.Item(ROBOT_NAME, robolink.ITEM_TYPE_ROBOT)
    moving_pose = robot.Pose()

    # Defining base targets
    home_target = RDK.Item('home_target')

    for i in range(len(layers)-1, 0, -1):
        thickness = layers[i] - layers[i-1]
        moving_pose = moving_pose*transl(0, 0, thickness)
        robot.MoveL(moving_pose)
        cv2.waitKey(1000)
    
    # Move back to home position
    robot.MoveJ(home_target)
