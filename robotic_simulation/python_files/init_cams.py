# You can also use the new version of the API:
from robodk import robolink    # RoboDK API
from robodk import robomath    # Robot toolbox
import cv2
import numpy as np
RDK = robolink.Robolink()

def init_cams():

    # Open side cam
    CAM_NAME = "side_camera"
    cam_item = RDK.Item(CAM_NAME,robolink.ITEM_TYPE_CAMERA)
    cam_item.setParam('Open', 1)

    # Open section cam
    CAM_NAME = "section_camera"
    cam_item = RDK.Item(CAM_NAME,robolink.ITEM_TYPE_CAMERA)
    cam_item.setParam('Open', 1)
