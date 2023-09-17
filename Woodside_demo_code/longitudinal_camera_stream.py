import cv2
import numpy as np
from object_measurement import object_measurement
import time

# Define system constants
real_frame_width = 0.75
cam_dis = 0.9
camera = cv2.VideoCapture(0)

def ResizeWithAspectRatio(image, width=None, height=None, inter=cv2.INTER_AREA):
    dim = None
    (h, w) = image.shape[:2]

    if width is None and height is None:
        return image
    if width is None:
        r = height / float(h)
        dim = (int(w * r), height)
    else:
        r = width / float(w)
        dim = (width, int(h * r))

    return cv2.resize(image, dim, interpolation=inter)



while cv2.waitKey(1) & 0xFF != ord('q'):
    ret, img = camera.read()
    # img = cv2.imread("test_images/side_90.jpg")

    output_image, obj_length, obj_width = object_measurement(img, real_frame_width, cam_dis, background=True, very_long=False, is_horizontal = True)

    print("=======================================================")
    print('The length of the pipe is: %.2f mm'%(obj_length))
    print('The diameter of the pipe is: %.2f mm'%(obj_width))
    cv2.imshow("output", ResizeWithAspectRatio(output_image, height=768))
    
camera.release()