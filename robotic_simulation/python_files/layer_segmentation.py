"""
This file uses various methods including thresholding and canny edge detection to identify circular edges in the image.
This is used to determine pipe layers.

"""

import cv2
import numpy as np
import math
import copy

# Rescaling image 
def rescaleFrame(frame, scale):
    """
    Function to rescale an image based on a scale ratio
    :param frame: image to rescale
    :param scale: positive float representing scale
    :returns: rescaled image
    """
    width = int(frame.shape[1] * scale)
    height = int(frame.shape[0] * scale)
    dimensions = (width, height)

    return cv2.resize(frame, dimensions, interpolation=cv2.INTER_CUBIC)


def detect_edges(binary_image_to_check, outer_radius, inner_radius, center_row, center_col, direction):
    """
    Function to detect edge radii by moving in different directions from the centre of the pipe. Transitions between
    black and white pixels are used to determine edges.

    :param binary_image_to_check: Binary image of pipe cross section
    :param outer_radius: INTEGER of the outer radius of the pipe cross section
    :param inner_radius: INTEGER of inner radius of the pipe cross section
    :param center_row: The INTEGER y coordinate of the center of the pipe image
    :param center_col: The INTEGER  x corrdinate of the center of the pipe image
    :param direction: STRING for the direction to go in from the center
        One of: ["up", "down", "left", "right", "up left", "up right", "down left", "down right"]
    :returns layer_edges: LIST containing INTEGERs of all detected raii
    """
    layer_edges = []
    last_pixel = 255
    for i in range(outer_radius):
        if direction == "up": 
            if binary_image_to_check[center_row - i, center_col] != last_pixel:
                if last_pixel == 255:
                    layer_edges.append(i)
                last_pixel = binary_image_to_check[center_row - i, center_col]
        elif direction == "down": 
            if binary_image_to_check[center_row + i, center_col] != last_pixel:
                if last_pixel == 255:
                    layer_edges.append(i)
                last_pixel = binary_image_to_check[center_row + i, center_col]
        elif direction == "left":
            if binary_image_to_check[center_row, center_col -i] != last_pixel:
                if last_pixel == 255:
                    layer_edges.append(i)
                last_pixel = binary_image_to_check[center_row, center_col -i]
        elif direction == "right":
            if binary_image_to_check[center_row, center_col +i] != last_pixel:
                if last_pixel == 255:
                    layer_edges.append(i)
                last_pixel = binary_image_to_check[center_row, center_col +i]
        elif direction == "up left":
            if binary_image_to_check[center_row -i, center_col -i] != last_pixel:
                if last_pixel == 255:
                    layer_edges.append(int(i*math.sqrt(2)))
                last_pixel = binary_image_to_check[center_row -i, center_col -i]
        elif direction == "up right":
            if binary_image_to_check[center_row -i, center_col +i] != last_pixel:
                if last_pixel == 255:
                    layer_edges.append(int(i*math.sqrt(2)))
                last_pixel = binary_image_to_check[center_row -i, center_col +i]
        elif direction == "down left":
            if binary_image_to_check[center_row +i, center_col -i] != last_pixel:
                if last_pixel == 255:
                    layer_edges.append(int(i*math.sqrt(2)))
                last_pixel = binary_image_to_check[center_row +i, center_col -i]
        elif direction == "down right":
            if binary_image_to_check[center_row +i, center_col +i] != last_pixel:
                if last_pixel == 255:
                    layer_edges.append(int(i*math.sqrt(2)))
                last_pixel = binary_image_to_check[center_row +i, center_col +i]

    return layer_edges


def check_circles_aux(binary_image_to_check, radius_to_check, center_row, center_col, real_radii, check_count):
    """
    Recursive function that checks whether a detected edge is a "real" edge or if it is just noise.
    
    :param binary_image_to_check: Binary image of pipe cross section
    :param radius_to_check: INTEGER of the radius to check legitimacy
    :param center_row: The INTEGER y coordinate of the center of the pipe image
    :param center_col: The INTEGER  x corrdinate of the center of the pipe image
    :param real_radii: LIST to append the INTEGER radius if it is legitimate
    :param check_count: INTEGER Recursive count

    """
    # Count the number of white pixels vs black pixels for the detect radius
    rows, cols = binary_image_to_check.shape
    white_pixels = 0
    black_pixels = 0
    for i in range(rows):
        for j in range(cols):
            if math.sqrt(abs(center_row-i)**2 + abs(center_col-j)**2) == radius_to_check-check_count:
                if binary_image_to_check[i, j] == 0:
                    black_pixels += 1
                else:
                    white_pixels += 1
    if white_pixels > 0.6*black_pixels:
        real_radii.append(radius_to_check-check_count)
    else:
        # Recursively call to check surrouding radii
        if check_count>0 and check_count< 4:
            check_circles_aux(binary_image_to_check, radius_to_check, center_row, center_col, real_radii, check_count+1)
        elif check_count<0 and check_count> -4:
            check_circles_aux(binary_image_to_check, radius_to_check, center_row, center_col, real_radii, check_count-1)
        elif check_count == 0:
            length = len(real_radii)
            check_circles_aux(binary_image_to_check, radius_to_check, center_row, center_col, real_radii, check_count+1)
            if len(real_radii) == length:
                check_circles_aux(binary_image_to_check, radius_to_check, center_row, center_col, real_radii, check_count-1)
    

def check_circles(binary_image_to_check, radii_to_check, center_row, center_col):
    """
    Function that checks whether a set of detected edges/radii are "real" or if it is just noise.

    :param binary_image_to_check: Binary image of pipe cross section
    :param radii_to_check: LIST of INTEGERs of the radii to check legitimacy
    :param center_row: The INTEGER y coordinate of the center of the pipe image
    :param center_col: The INTEGER  x corrdinate of the center of the pipe image
    :returns real_radii: a LIST containing INTEGERs of all "real" and verified edges
    """
    real_radii = []
    for k in range(len(radii_to_check)):
        check_circles_aux(binary_image_to_check, radii_to_check[k], center_row, center_col, real_radii, 0)
    return real_radii


def average_nearby(radii, threshold):
    """
    Given a LIST of integers, it averages and combines nearby integers.

    :param radii: A LIST of INTEGERs containing detected radii values from varius functions
    :param threshold: INTEGER threshold for merging. Example: threshold=5 means all integers within 5 of each other gets merged
    :returns circle_radii: A LIST of remaining integers after merging and averaging
    """
    sorted_radii = sorted(radii)
    circle_radii = []
    total = sorted_radii[0]
    count = 1
    if len(sorted_radii) > 1:
        for i in range(len(sorted_radii)-1):
            if sorted_radii[i+1] - sorted_radii[i] <= threshold:
                total += sorted_radii[i+1]
                count += 1
            else:
                circle_radii.append(int(total/count))
                total = sorted_radii[i+1]
                count = 1
            
            if i == len(sorted_radii)-2: # Last element
                circle_radii.append(int(total/count))
    else:
        circle_radii.append(sorted_radii[0])
    return circle_radii


def count_radius(binary_image_to_check, center_row, center_col, outer_radius, threshold):
    """
    Function to detect circle radii given canny edge detection

    :param binary_image_to_check: Binary image output from canny edge detection
    :param center_row: The INTEGER y coordinate of the center of the pipe image
    :param center_col: The INTEGER  x corrdinate of the center of the pipe image
    :param outer_radius: INTEGER of the outer radius of the pipe cross section
    :param threshold: FLOAT threshold for a valid circle. Example: threshold=0.2 means the circle only
    needs to be 20% complete to be detected.

    :returns circle_radii: A LIST of all detected circle radii INTEGERs
    """
    rows, cols = binary_image_to_check.shape
    radius_count = [0]*(min(rows, cols)//2)
    for i in range(rows):
        for j in range(cols):
            if binary_image_to_check[i, j] == 255:
                radius_count[int(math.sqrt(abs(center_row-i)**2 + abs(center_col-j)**2))] += 1

    possible_radii = []
    for i in range(len(radius_count)):
        if i <= outer_radius+2 and radius_count[i] > threshold*2*3.1415*i:
            possible_radii.append(i)
    
    circle_radii = average_nearby(possible_radii, 5)
    return circle_radii


def mask_image(image):
    """
    Function to mask the pipe cross section by removing background and inside of the pipe.
    
    :param image: OpenCV BGR image of the pipe cross section
    :returns image_masked: Masked image of the pipe
    :returns x1: INTEGER x coordinate of the center of the outer circle of the pipe
    :returns y1: INTEGER y coordinate of the center of the outer circle of the pipe
    :returns outer_radius: outer radius INTEGER of pipe
    :returns x2: INTEGER x coordinate of the center of the inner circle of the pipe
    :returns y2: INTEGER y coordinate of the center of the inner circle of the pipe
    :returns inner_radius: inner radius INTEGER of pipe
    """

    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

    # Converting blurred image into binary image
    thresh = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)[1]
    # cv2.imshow('thresh',thresh)

    # Perform morpholgical operations (patch up small holes etc)
    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (9, 9))
    close = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE, kernel, iterations=1)
    # cv2.imshow('close0', close)
    opening = cv2.morphologyEx(close, cv2.MORPH_OPEN, kernel, iterations=2)
    close = cv2.morphologyEx(opening, cv2.MORPH_CLOSE, kernel, iterations=1)
    # cv2.imshow('close', close)

    # Find the edges of the morphed/closed binary image
    canny = cv2.Canny(close, 0, 0)
    # cv2.imshow("canny", canny)

    # Use the canny edges to list out the contours
    contours, hierarchy = cv2.findContours(canny, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)

    # Getting the contours of the outer and inner outlines of the pipe
    contours = sorted(contours, key=cv2.contourArea, reverse=True)
    outer_contour = contours[0]
    inner_contour = contours[2]
    circle_mask = np.zeros(image.shape, dtype = "uint8")

    # Creating the outer circular mask
    (x1,y1),outer_radius = cv2.minEnclosingCircle(outer_contour)
    # Creating the circular mask
    outer_mask = np.zeros(image.shape, dtype = "uint8")
    cv2.circle(outer_mask, (int(x1), int(y1)), int(outer_radius), (255, 255, 255), 2)
    canny = cv2.Canny(outer_mask, 0, 0)
    contours, hierarchy = cv2.findContours(canny, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
    biggest_contour = max(contours, key=cv2.contourArea)
    cv2.fillPoly(circle_mask, pts=[biggest_contour], color=[255, 255, 255])

    # Creating the inner circular mask
    (x2,y2),inner_radius = cv2.minEnclosingCircle(inner_contour)
    inner_mask = np.zeros(image.shape, dtype = "uint8")
    cv2.circle(inner_mask, (int(x2), int(y2)), int(inner_radius), (255, 255, 255), 2)
    canny = cv2.Canny(inner_mask, 0, 0)
    contours, hierarchy = cv2.findContours(canny, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
    biggest_contour = max(contours, key=cv2.contourArea)
    cv2.fillPoly(circle_mask, pts=[biggest_contour], color=[0, 0, 0])

    # Mask with original image
    image_masked = cv2.bitwise_and(image, circle_mask)
    # cv2.imshow("image_masked", image_masked)

    return image_masked, x1, y1, outer_radius, x2, y2, inner_radius


def get_layers_from_binary(image_masked, center_row, center_col, outer_radius, inner_radius):
    """
    Function to detect the layer edges using the thresholding method.

    :param image_masked: Masked image of the pipe cross section. (background and inside removed)
    :param center_row: The INTEGER y coordinate of the center of the pipe image
    :param center_col: The INTEGER  x corrdinate of the center of the pipe image
    :param outer_radius: INTEGER of the outer radius of the pipe cross section
    :param inner_radius: INTEGER of inner radius of the pipe cross section

    :returns layer_edges: LIST of all the detected layer edges
    :returns binary_image: Binary image used to detect the edges
    :returns  outer_radius: Detected outer radius of the pipe cross section
    """

    # Converting masked image into binary image again
    gray = cv2.cvtColor(image_masked, cv2.COLOR_BGR2GRAY)

    threshold_values = list(reversed([120, 125, 130, 135, 140, 145, 150, 155]))
    layer_edges = []
    for t in range(len(threshold_values)):
        thresh2 = cv2.threshold(gray, threshold_values[t], 255, cv2.THRESH_BINARY_INV)[1]
        # cv2.imshow('thresh2',thresh2)

        # Perform morpholgical operations (increase gaps between layers)
        kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (3, 3))
        close2 = cv2.morphologyEx(thresh2, cv2.MORPH_CLOSE, kernel, iterations=1)
        if t == 0:
            binary_image = close2

        # detect the layer edges going in different directions
        layer_edges_radii_up = detect_edges(close2, outer_radius, inner_radius, center_row, center_col, "up")
        layer_edges_radii_down = detect_edges(close2, outer_radius, inner_radius, center_row, center_col, "down")
        layer_edges_radii_left = detect_edges(close2, outer_radius, inner_radius, center_row, center_col, "left")
        layer_edges_radii_right = detect_edges(close2, outer_radius, inner_radius, center_row, center_col, "right")
        layer_edges_radii_up_left = detect_edges(close2, outer_radius, inner_radius, center_row, center_col, "up left")
        layer_edges_radii_up_right  = detect_edges(close2, outer_radius, inner_radius, center_row, center_col, "up right")
        layer_edges_radii_down_left = detect_edges(close2, outer_radius, inner_radius, center_row, center_col, "down left")
        layer_edges_radii_down_right = detect_edges(close2, outer_radius, inner_radius, center_row, center_col, "down right")

        # choose the one with most detected edges
        layer_edges_copy = copy.deepcopy(layer_edges)
        layer_edges = [layer_edges_copy, layer_edges_radii_up, layer_edges_radii_down, layer_edges_radii_left, layer_edges_radii_right, layer_edges_radii_up_left, layer_edges_radii_up_right, layer_edges_radii_down_left, layer_edges_radii_down_right]
        layer_edges = layer_edges[np.argmax([len(layer_edges_copy), len(layer_edges_radii_up), len(layer_edges_radii_down), len(layer_edges_radii_left), len(layer_edges_radii_right), len(layer_edges_radii_up_left), len(layer_edges_radii_up_right), len(layer_edges_radii_down_left), len(layer_edges_radii_down_right)])]
        direction_index = np.argmax([len(layer_edges_copy), len(layer_edges_radii_up), len(layer_edges_radii_down), len(layer_edges_radii_left), len(layer_edges_radii_right), len(layer_edges_radii_up_left), len(layer_edges_radii_up_right), len(layer_edges_radii_down_left), len(layer_edges_radii_down_right)])
        if direction_index != 0:
            direction = direction_index
            binary_image = copy.deepcopy(close2)

    # Find the outer diameter of pipe cross section
    if direction == 1: # up
        # print("detect direction = up")
        for i in range(outer_radius+5, -1, -1): # the +10 is so that we end up outside the pipe
            if binary_image[center_row - i, center_col] == 0:
                layer_edges.append(i)
                outer_radius = i
                break
    elif direction == 2: #down
        # print("detect direction = down")
        for i in range(outer_radius+5, -1, -1): # the +10 is so that we end up outside the pipe
            if binary_image[center_row + i, center_col] == 0:
                layer_edges.append(i)
                outer_radius = i
                break
    elif direction == 3: #left
        # print("detect direction = left")
        for i in range(outer_radius+5, -1, -1): # the +10 is so that we end up outside the pipe
            if binary_image[center_row, center_col-i] == 0:
                layer_edges.append(i)
                outer_radius = i
                break
    elif direction == 4: #right
        # print("detect direction = right")
        for i in range(outer_radius+5, -1, -1): # the +10 is so that we end up outside the pipe
            if binary_image[center_row, center_col+i] == 0:
                layer_edges.append(i)
                outer_radius = i
                break
    elif direction == 5: #up left
        # print("detect direction = up left")
        for i in range(outer_radius+5, -1, -1): # the +10 is so that we end up outside the pipe
            if binary_image[center_row-i, center_col-i] == 0:
                layer_edges.append(int((i)*math.sqrt(2)))
                outer_radius = int((i)*math.sqrt(2))
                break
    elif direction == 6: #up right
        # print("detect direction = up right")
        for i in range(outer_radius+5, -1, -1): # the +10 is so that we end up outside the pipe
            if binary_image[center_row-i, center_col+i] == 0:
                layer_edges.append(int((i)*math.sqrt(2)))
                outer_radius = int((i)*math.sqrt(2))
                break
    elif direction == 7: #down left
        # print("detect direction = down left")
        for i in range(outer_radius+5, -1, -1): # the +10 is so that we end up outside the pipe
            if binary_image[center_row+i, center_col-i] == 0:
                layer_edges.append(int((i)*math.sqrt(2)))
                outer_radius = int((i)*math.sqrt(2))
                break
    elif direction == 8: #down right
        # print("detect direction = down right")
        for i in range(outer_radius+5, -1, -1): # the +10 is so that we end up outside the pipe
            if binary_image[center_row+i, center_col+i] == 0:
                layer_edges.append(int((i)*math.sqrt(2)))
                outer_radius = int((i)*math.sqrt(2))
                break
    return layer_edges, binary_image, outer_radius


def get_layers_from_canny(image_masked, center_row, center_col, outer_radius, blur, detection_threshold):
    """
    Function to detect the layer edges using the Canny Edges method.

    :param image_masked: Masked image of the pipe cross section. (background and inside removed)
    :param center_row: The INTEGER y coordinate of the center of the pipe image
    :param center_col: The INTEGER  x corrdinate of the center of the pipe image
    :param outer_radius: INTEGER of the outer radius of the pipe cross section
    :param blur: INTEGER ODD blur value used in GaussianBlur
    :param detection_threshold: FLOAT threshold for a valid circle. Example: detection_threshold=0.2 means the circle only
    needs to be 20% complete to be detected.

    :returns canny_circles: LIST of all the detected layer edges using Canny Edges method
    :returns canny_edges: Image of all the edges detected by Canny Edge Detector
    """
    kernel = np.array([[0, -1, 0],
                    [-1, 5,-1],
                    [0, -1, 0]], )
    image_sharp = cv2.filter2D(src=image_masked, ddepth=-1, kernel=kernel)
    # cv2.imshow('sharpened image', image_sharp)

    # Apply Bilateral Blurring (to reduce noise while keeping edges sharp)
    # blurred = cv2.bilateralFilter(image_sharp, 11, 75, 75)
    blurred = cv2.GaussianBlur(image_sharp,(blur,blur),0)

    canny_edges = cv2.Canny(blurred, 45, 135)

    canny_circles = count_radius(canny_edges, center_row, center_col, outer_radius, detection_threshold)
    return canny_circles, canny_edges
######################################### HELPER FUNCTIONS END ###########################################################


def layer_segmentation(image, real_frame_dis):
    """
    Function to segment the layers of the pipe by detecting the layer edges

    :param image_file_path: openCV image
    :returns image_masked: Image with all the layer edges drawn
    :returns radii: A LIST containing all detect layer edge radius
    """
    # Resize image
    image = cv2.resize(image, (768,768), interpolation=cv2.INTER_CUBIC)

    # Extract only the pipe cross section from image
    image_masked, x1, y1, outer_radius, x2, y2, inner_radius = mask_image(image)

    center_row = int(y2)
    center_col = int(x2)
    outer_radius = int(outer_radius)
    inner_radius = int(inner_radius)

    # Get Layer radii from Binary Method
    layer_edges , binary_image, outer_radius = get_layers_from_binary(image_masked, center_row, center_col, outer_radius, inner_radius)
    radii = check_circles(binary_image, layer_edges, center_row, center_col)

    # Creating the outline mask
    outer_mask = np.zeros(image.shape, dtype = "uint8")
    cv2.circle(outer_mask, (center_col, center_row), int(outer_radius), (255, 255, 255), 2)
    canny = cv2.Canny(outer_mask, 0, 0)
    contours, hierarchy = cv2.findContours(canny, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
    biggest_contour = max(contours, key=cv2.contourArea)
    cv2.fillPoly(outer_mask, pts=[biggest_contour], color=[255, 255, 255])
    image_masked = cv2.bitwise_and(image_masked, outer_mask)

    # Get Layer radii from Canny Method
    canny_circles, canny_image = get_layers_from_canny(image_masked, center_row, center_col, outer_radius, blur=11, detection_threshold=0.16)
    

    # Combine and average results
    radii = average_nearby(radii + canny_circles, 3)

    img_width = image.shape[1] # image width
    each_pixel_dis = real_frame_dis/img_width
    real_radii = [x*each_pixel_dis for x in radii]

    for i in range(len(radii)):
        cv2.circle(image_masked,(center_col, center_row),radii[i],(0,0,255),2)
        layer_string = "edge {:.0f}={:.2f}mm".format(i+1, real_radii[i]*1000)
        # cv2.putText(image_masked, layer_string, tuple([round(center_col+radii[i]*math.cos(math.pi/4)), round(center_row-radii[i]*math.sin(math.pi/4))]), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 0, 0), 1)
        # print("Layer Edge %i Radius: %.2f mm" % (i+1, radii[i]*length_per_pixel*1000))
    

    
    # cv2.imshow('original image', image)
    # cv2.imshow("binary image", binary_image)
    # print("Layers detected by Binary Method: " + str(radii))
    # cv2.imshow('canny_edges', canny_image)
    # print("Layers detected by Canny Method: " + str(canny_circles))
    # cv2.imshow("layers", image_masked)
    # print("Combined All Layers: " + str(radii))

    # cv2.waitKey(0)
    # cv2.destroyAllWindows()

    return image_masked, real_radii