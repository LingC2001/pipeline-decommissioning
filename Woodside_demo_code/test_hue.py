import cv2
import numpy as np

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

image = cv2.imread("test_images/side.jpg")

# Mask dark pixels?
hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
mask = cv2.inRange(hsv, (0, 30, 150), (255, 255, 255))
mask = cv2.bitwise_not(mask)
image[np.where(mask)] = (0,0,0)


new_img = image.copy()
new_img[np.where(mask)] = (0,0,0)

cv2.imshow("mask", ResizeWithAspectRatio(new_img, width=768))
cv2.waitKey(0)