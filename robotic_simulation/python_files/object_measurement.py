"""
This script uses Thresholding and Canny Edge Detection method to detect measure the width/length of the pipe in an image.

"""

import cv2
import numpy as np
import math

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


def toBinary(image):
    """
    Function to convert an image to binary version using OTSU thresholding

    :param image: The OpenCV cv::Mat BGR image
    :returns: Binary image
    """
    # Convert image to grayscale
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    # Apply Bilateral Blurring (to reduce noise while keeping edges sharp)
    blur = cv2.bilateralFilter(gray, 5, 75, 75)
    # Converting blurred image into binary image
    thresh = cv2.threshold(blur, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)[1]
    # cv2.imshow('thresh',thresh)

    return thresh


def morphClose(binary):
    """
    Function to apply morphological operations to a binary image.

    :param binary: Binary image
    :returns: Binary image with holes closed
    """
    # Perform morpholgical operations
    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (7, 7))
    opening = cv2.morphologyEx(binary, cv2.MORPH_OPEN, kernel, iterations=1)
    close = cv2.morphologyEx(opening, cv2.MORPH_CLOSE, kernel, iterations=1)
    # cv2.imshow('close', close)

    return close


def findBestContour(closed_binary):
    """
    Function to find the contours in the closed binary image that best fit the shape of a pipe (rectangle)

    :param closed_bianry: Binary image that has undergone closing and opening morphological operations
    :returns: large contour that best fits a rectangular shape
    """
    # Find the edges of the morphed/closed binary image
    canny = cv2.Canny(closed_binary, 0, 0)
    #cv2.imshow("canny", canny)

    # Use the canny edges to list out the contours
    contours, hierarchy = cv2.findContours(canny, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)

    # Finding the contour that fits the best into a rectangle
    for c in contours:
        if cv2.contourArea(c) > 15000:
            best_fit_contour = c
            break

    min_area = math.inf
    for c in contours:
        if cv2.contourArea(c) > 15000:
            # Fit contours into a rectangular box with minimum area
            min_rect = cv2.minAreaRect(c)

            # calculate box points
            box = cv2.boxPoints(min_rect)
            box = np.intp(box)

            # calculate area of rectangle
            x = math.sqrt(((box[0][0] - box[1][0]) ** 2) + ((box[0][1] - box[1][1])) ** 2)
            y = math.sqrt(((box[0][0] - box[3][0]) ** 2) + ((box[0][1] - box[3][1])) ** 2)
            box_area = x*y

            if (box_area - cv2.contourArea(c))/box_area < min_area:
                best_fit_contour = c
                min_area = (box_area - cv2.contourArea(c))/box_area

    return best_fit_contour


def extractPipe(image, best_fit_contour):
    """
    Function to uses a rectangular mask around the best fit contour on the image to extract the general area of the pipe.
    Then thresholding and canny edge detection is used again to finely extract the pipe itself from the surrounding background.

    :param image: original image
    :param best_fit_contour: The best fit contour in the image from Function findBestContour()
    :returns masked_image: The extract image of the pipe
    :returns biggest_contour: The contour of the pipe
    """
    # Fit contour into a rectangular box with minimum area
    min_rect = cv2.minAreaRect(best_fit_contour)

    # calculate box points
    box = cv2.boxPoints(min_rect)
    box = np.intp(box)

    # bit mask and extract object
    mask = np.zeros(image.shape, dtype = "uint8")
    cv2.fillPoly(mask, pts=[box], color=(255,255,255))
    masked_image = cv2.bitwise_and(image, mask)
    # cv2.imshow("masked_image", masked_image)

    # Applying thresholding again
    gray_masked = cv2.cvtColor(masked_image, cv2.COLOR_BGR2GRAY)
    blur_masked = cv2.GaussianBlur(gray_masked, (5, 5), cv2.BORDER_DEFAULT)
    threshold2 = cv2.threshold(blur_masked, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)[1]
    # cv2.imshow('thresh2', threshold2)

    # Perform morpholgical operations
    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (21, 21))
    opening2 = cv2.morphologyEx(threshold2, cv2.MORPH_OPEN, kernel, iterations=1)
    close2 = cv2.morphologyEx(opening2, cv2.MORPH_CLOSE, kernel, iterations=1)
    # cv2.imshow('close2', close2)

    # Find the edges of the morphed/closed binary image
    canny2 = cv2.Canny(close2, 0, 0)
    # cv2.imshow("canny2", canny2)

    # Use the canny edges to list out the contours
    contours2, hierarchy2 = cv2.findContours(canny2, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)

    # Finding the contour with the largest area
    biggest_contour = contours2[0]
    for c in contours2:
        if cv2.contourArea(c) > cv2.contourArea(biggest_contour):
            biggest_contour = c

    # bit mask and extract object
    mask2 = np.zeros(image.shape, dtype = "uint8")
    cv2.fillPoly(mask2, pts=[biggest_contour], color=(255,255,255))
    masked_image = cv2.bitwise_and(image, mask2)
    
    return masked_image, biggest_contour


def drawBox(image, contour, each_pixel_dis, cam_dis):
    """
    Function that draws the dimension box around the pipe on an image and returns the length and diameter of the pipe.

    :param image: pipe image to draw on
    :param contour: contour of the pipe in the image
    :param each_pixel_dis: the real life distance of each pixel in the image
    :returns image: The image with dimension box drawn
    :returns obj_length: length of pipe in metres
    :returns obj_width: diameter/width of pipe in metres
    
    """
    # Fit large contours into a rectangular box with minimum area
    min_rect = cv2.minAreaRect(contour)
    # display box
    box2 = cv2.boxPoints(min_rect)
    box2 = np.intp(box2)
    cv2.drawContours(image, [box2], 0, (0, 0, 255))

    # Calculate and display length and width of box/object in the image
    x = math.sqrt(((box2[0][0] - box2[1][0]) ** 2) + ((box2[0][1] - box2[1][1])) ** 2)
    y = math.sqrt(((box2[0][0] - box2[3][0]) ** 2) + ((box2[0][1] - box2[3][1])) ** 2)
    # print(box2)
    if x >= y: # long edge is 0-1 and 2-3
        length = x
        width = y
        long_edge_is_01 = True
    else: # long edge is 0-3 and 1-2
        length = y
        width = x
        long_edge_is_01 = False
    
    real_length = correct_length(box2, long_edge_is_01, width, image.shape[1], each_pixel_dis, cam_dis)

    # calculate object width and length
    obj_length = each_pixel_dis*real_length
    obj_width = each_pixel_dis*width 

    dimensions_string = "angle_corrected_length={:.4f}m, width={:.4f}m".format(obj_length, obj_width)
    cv2.putText(image, dimensions_string, tuple(box2[0] + [10, -10]), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 255, 0), 1)
    cv2.drawContours(image, [box2], 0, (0, 0, 255))
    # cv2.imshow("dimensions", image)

    return image, obj_length, obj_width

def correct_length(box_coords, long_side_is_01, diameter, sqr_img_size, each_pixel_dis, cam_dis):
    """
    Function to calculate the length of the pipe given the diameter of the pipe accounting for
    the angle of the camera.
    
    """
    if long_side_is_01:
        left_mid_point = [(box_coords[0][0] + box_coords[3][0])//2, (box_coords[0][1] + box_coords[3][1])//2]
        right_mid_point = [(box_coords[1][0] + box_coords[2][0])//2, (box_coords[1][1] + box_coords[2][1])//2]
    else:
        left_mid_point = [(box_coords[0][0] + box_coords[1][0])//2, (box_coords[0][1] + box_coords[1][1])//2]
        right_mid_point = [(box_coords[2][0] + box_coords[3][0])//2, (box_coords[2][1] + box_coords[3][1])//2]

    if left_mid_point[0] != right_mid_point[0]:
        m = (left_mid_point[1]-right_mid_point[1])/(left_mid_point[0]-right_mid_point[0])
        y = sqr_img_size//2
        x1 = left_mid_point[0]
        y1 = left_mid_point[1]
        pipe_centre = [(y-y1+m*x1)/m, y]
    else:
        pipe_centre = [left_mid_point[0], sqr_img_size//2]

    left_dis_to_centre = math.sqrt(((left_mid_point[0] - pipe_centre[0]) ** 2) + ((left_mid_point[1] - pipe_centre[1])) ** 2)
    right_dis_to_centre = math.sqrt(((right_mid_point[0] - pipe_centre[0]) ** 2) + ((right_mid_point[1] - pipe_centre[1])) ** 2)
    z = cam_dis//each_pixel_dis
    left_angle = np.arctan(left_dis_to_centre/z)
    right_angle = np.arctan(right_dis_to_centre/z)

    real_left = left_dis_to_centre - diameter*np.tan(left_angle)
    real_right = right_dis_to_centre - diameter*np.tan(right_angle)
    return real_left + real_right



def object_measurement(image, real_frame_width, cam_dis):
    """
    Function to measure the dimensions of a pipe in an image
    :param image: openCV image
    :param real_frame_width: Width of the image frame in real life in meters
    :returns output_image: Image with dimensions drawn
    :returns obj_length: length of pipe in metres
    :returns obj_width: diameter/width of pipe in metres

    """
    # Rescale image so that the bigger dimension is 768
    scale = 768/max(image.shape[0],image.shape[1])
    image = rescaleFrame(image, scale)

    img_width = image.shape[1] # image width
    each_pixel_dis = real_frame_width/img_width

    # Performing operations to find the dimensions of the pipe. Refer to each of the function documentations
    thresh = toBinary(image)
    close = morphClose(thresh)
    best_fit_contour = findBestContour(close)
    mask_image, contour = extractPipe(image, best_fit_contour)
    output_image, obj_length, obj_width = drawBox(mask_image, contour, each_pixel_dis, cam_dis)

    # # Showing results
    # print('The length of the pipe is: %.4f m'%(obj_length))
    # print('The diameter of the pipe is: %.4f m'%(obj_width))
    # cv2.imshow('original image', image)
    # cv2.imshow("dimensions", output_image)
    # cv2.waitKey(0)
    # cv2.destroyAllWindows()

    return output_image, obj_length, obj_width






