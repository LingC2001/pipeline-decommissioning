# You can also use the new version of the API:
from robodk import robolink    # RoboDK API
from robodk import robomath    # Robot toolbox
import cv2
import numpy as np

def get_section_cam():
    RDK = robolink.Robolink()
    # Close any open 2D camera views
    CAM_NAME = "section_camera"
    cam_item = RDK.Item(CAM_NAME,robolink.ITEM_TYPE_CAMERA)
    cam_item.setParam('Open', 1)

    bytes_img = RDK.Cam2D_Snapshot("", cam_item)
    if isinstance(bytes_img, bytes) and bytes_img != b'':
        nparr = np.frombuffer(bytes_img, np.uint8)
        img_socket = cv2.imdecode(nparr, cv2.IMREAD_COLOR)

    return img_socket
# cv2.imshow("side_image", img_socket)
# cv2.imwrite("side.jpg", img_socket)
# cv2.waitKey(0)
# cv2.destroyAllWindows()