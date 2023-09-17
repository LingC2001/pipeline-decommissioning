import matplotlib.pyplot as plt

from skimage import data, color, img_as_ubyte
from skimage.feature import canny
from skimage.transform import hough_ellipse
from skimage.draw import ellipse_perimeter
import cv2
import numpy as np

# Load picture, convert to grayscale and detect edges
image = cv2.imread('p39.jpg')
image = cv2.resize(image, (1024,1024), interpolation=cv2.INTER_CUBIC)

kernel = np.array([[0, -1, 0],
                    [-1, 5,-1],
                    [0, -1, 0]])
image_sharp = cv2.filter2D(src=image, ddepth=-1, kernel=kernel)

blurred = cv2.GaussianBlur(image_sharp,(11,11),0)

canny_edges = cv2.Canny(blurred, 45, 135)

save_path = "canny_images/"
cv2.imwrite(save_path+"s3.jpg", image_sharp)

"""
image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
image_gray = cv2.cvtColor(image_rgb, cv2.COLOR_BGR2GRAY)
edges = canny(image_gray, sigma=2.0,
              low_threshold=0.55, high_threshold=0.8)
# Perform a Hough Transform
# The accuracy corresponds to the bin size of a major axis.
# The value is chosen in order to get a single high accumulator.
# The threshold eliminates low accumulators
result = hough_ellipse(edges, accuracy=20, threshold=250,
                       min_size=100, max_size=120)
result.sort(order='accumulator')

# Estimated parameters for the ellipse
best = list(result[-1])
yc, xc, a, b = [int(round(x)) for x in best[1:5]]
orientation = best[5]

# Draw the ellipse on the original image
cy, cx = ellipse_perimeter(yc, xc, a, b, orientation)
image_rgb[cy, cx] = (0, 0, 255)
# Draw the edge (white) and the resulting ellipse (red)
edges = color.gray2rgb(img_as_ubyte(edges))
edges[cy, cx] = (250, 0, 0)

cv2.imshow("img_rbg", image_rgb)
cv2.imshow("edges", edges)
"""


