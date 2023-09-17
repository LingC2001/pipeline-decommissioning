from robodk import robolink    # RoboDK API
from robodk import robomath    # Robot toolbox
from robodk.robomath import *

import numpy as np
import cv2
RDK = robolink.Robolink()

def robot_move_dimensions(length, diameter):
    # Defining robot
    ROBOT_NAME = "KUKA LBR iiwa 14 R820"
    robot = RDK.Item(ROBOT_NAME, robolink.ITEM_TYPE_ROBOT)
    # Resetting Robot
    HOME = [0.000, 0.000, 0.000, -90.000, 0.000, 0.000, 0.000]
    robot.setJoints(HOME)

    # Defining base targets
    home_target = RDK.Item('home_target')
    target = RDK.Item('above_sensor_line')


    # Move tool to target above sensor line
    robot.MoveJ(target)

    # Compute robot path
    target_pose = target.Pose()
    xyz_ref = target_pose.Pos()
    # print(xyz_ref)

    # moving to touch pipe
    moving_pose = target.Pose()*transl(0,0,xyz_ref[2]-100-diameter)
    robot.MoveL(moving_pose)

    # moving to left side of pipe
    moving_pose = moving_pose*transl(0, 0, diameter/2)
    prev_pose = moving_pose*transl(-diameter*np.cos(np.pi/4)/2, 0, -diameter*np.sin(np.pi/4)/2)*roty(np.pi/8)
    moving_pose = prev_pose*roty(-np.pi/8)*transl(diameter*np.cos(np.pi/4)/2, 0, diameter*np.sin(np.pi/4)/2)
    moving_pose = moving_pose*transl(-diameter/2, 0, 0)*roty(np.pi/8)
    robot.MoveC(prev_pose, moving_pose)

    # move along the length of the pipe
    moving_pose = moving_pose*transl(0, -length, 0)
    robot.MoveL(moving_pose)

    # move to the top of pipe
    moving_pose = moving_pose*roty(-np.pi/8)*transl(diameter/2, 0, 0)
    prev_pose =  moving_pose*transl(-diameter*np.cos(np.pi/4)/2, 0, -diameter*np.sin(np.pi/4)/2)*roty(np.pi/8)
    moving_pose = prev_pose*roty(-np.pi/8)*transl(diameter*np.cos(np.pi/4)/2, 0, diameter*np.sin(np.pi/4)/2)
    moving_pose = moving_pose*transl(0, 0, -diameter/2)
    robot.MoveC(prev_pose, moving_pose)

    # move back along the length of the pipe
    moving_pose = moving_pose*transl(0, length, 0)
    robot.MoveL(moving_pose)




# ref_frame = RDK.Item('Conveyor Belt (2m) Base', robolink.ITEM_TYPE_FRAME)
# print(ref_frame.Pose())
# print(robot.Tool())