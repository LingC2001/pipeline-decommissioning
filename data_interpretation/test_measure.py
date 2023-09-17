import cv2
import numpy as np
from object_measurement import object_measurement
from layer_segmentation import layer_segmentation

img_long = cv2.imread("pipe_short.png")
img_cross = cv2.imread("section.jpg")

real_frame_width = 0.5
cam_dis = 0.2
output_image, obj_length, obj_width = object_measurement(img_long, real_frame_width, cam_dis)

print('The length of the pipe is: %.2f mm'%(obj_length))
print('The diameter of the pipe is: %.2f mm'%(obj_width))
cv2.imshow("output", output_image)
cv2.waitKey(0)
cv2.destroyAllWindows()



cam_dis = 0.2
real_frame_dis = 2*cam_dis*np.tan(30*np.pi/180)
image_layers , radii = layer_segmentation(img_cross, real_frame_dis)

# Showing results
cv2.imshow("Pipe Layers", image_layers)
print("Layers Radii: " + str(radii))