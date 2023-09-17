from robodk import robolink    # RoboDK API
from robodk import robomath    # Robot toolbox
from robodk.robomath import *

import numpy as np
import cv2
RDK = robolink.Robolink()

def init_robot():
    # Defining robot
    ROBOT_NAME = "KUKA LBR iiwa 14 R820"
    robot = RDK.Item(ROBOT_NAME, robolink.ITEM_TYPE_ROBOT)
    # Resetting Robot
    HOME = [0.000, 0.000, 0.000, -90.000, 0.000, 0.000, 0.000]
    robot.setJoints(HOME)