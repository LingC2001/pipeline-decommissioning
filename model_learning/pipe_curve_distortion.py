"""
This script contains functions to offset image pixels to introduce artificial curvature.

"""

import cv2
import numpy as np
import math

def pipe_curve_distortion(file_path, images_to_distort, save_path, sine_curve=False, full_sine_curve=False, x4_curve=False, parabolic_curve=False, half_parabolic_curve=False):
    """
    Function to call the curve function to distort a list of images of pipes with different mathematical offset functions.

    :param file_path: file path STRING from current directory to the directory containing all the images. 
        Example: "dir1/imgs/"
    :param images_to_distort: LIST containing all the file names of images in the file path directory. 
        Example: ["1.jpg", "2.png"]
    :param save_path: file path STRING from current directory to the directory to save the augmented images. 
        Example: "dir1/saved_imgs/"
    :param sine_curve: BOOLEAN value on whether to apply the half sine curve offset
    :param full_sine_curve: BOOLEAN value on whether to apply the full sine curve offset
    :param x4_curve: BOOLEAN value on whether to apply the x^4 curve offset
    :param parabolic_curve: BOOLEAN value on whether to apply the x^2 curve offset
    :param half_parabolic_curve: BOOLEAN value on whether to apply the half parabolic curve offset

    """
    # Iterating over all images to distort
    for i in range(len(images_to_distort)):
        # Reading the image
        img = cv2.imread(file_path + images_to_distort[i])

        # Cropping image to square, and resizing to (1024, 1024)
        crop_size = min(img.shape[0], img.shape[1])
        start_row = (img.shape[0]-crop_size)//2
        end_row = (img.shape[0]-crop_size)//2 + crop_size
        start_col = (img.shape[1]-crop_size)//2
        end_col = (img.shape[1]-crop_size)//2 + crop_size
        img = img[start_row:end_row, start_col:end_col]

        img = cv2.resize(img, (768, 768), interpolation=cv2.INTER_CUBIC)

        # Applying selected curve operations
        if sine_curve:
            curve(img, curve_type="sine_curve", save_file_name=save_path+images_to_distort[i][:len(images_to_distort[i])-4], fill_black = True)
        if full_sine_curve:
            curve(img, curve_type="full_sine_curve", save_file_name=save_path+images_to_distort[i][:len(images_to_distort[i])-4], fill_black = True)
        if x4_curve:
            curve(img, curve_type="x4_curve", save_file_name= save_path+images_to_distort[i][:len(images_to_distort[i])-4], fill_black = True)
        if parabolic_curve:
            curve(img, curve_type="parabolic_curve", save_file_name= save_path+images_to_distort[i][:len(images_to_distort[i])-4], fill_black = True)
        if half_parabolic_curve:
            curve(img, curve_type="half_parabolic_curve", save_file_name= save_path+images_to_distort[i][:len(images_to_distort[i])-4], fill_black = True)


def curve(img, curve_type, save_file_name, fill_black = False):
    """
    Function to apply an mathematical curve offset to an image.

    :param img: OpenCV cv::MAT image
    :param curve_type: STRING of the type of curve offset to be used. 
        Can be one of: ["sine_curve", "full_sine_curve", "x4_curve", "parabolic_curve", "half_parabolic_curve"]
    :param save_file_name: STRING save path + image name (excluding file extension)
        Example: "dir1/save_imgs/img1"
    :param fill_black: BOOLEAN value on whether to wrap the pixels around the image. 
        Pixels that get offset out of the image wrap around to the other side and fills in the black.

    """
    rows, cols = img.shape[:2]

    # Setting amplitude levels
    if curve_type == "sine_curve":
        amplitude_levels = [0.1, 0.15, 0.2, 0.3]
    elif curve_type == "full_sine_curve":
        amplitude_levels = [0.1, 0.15, 0.2, 0.3]
    elif curve_type == "x4_curve":
        amplitude_levels = [0.15, 0.3, 0.5, 0.7]
    elif curve_type == "parabolic_curve":
        amplitude_levels = [0.2, 0.4, 0.6, 0.8]
    elif curve_type == "half_parabolic_curve":
        amplitude_levels = [0.1, 0.3, 0.5, 0.7]

    # Creating distorted images for each amplitude level
    for k in range(len(amplitude_levels)):
        img_output = np.zeros(img.shape, dtype=img.dtype) # Blank output image

        # Iterating over each pixexl of the image to distort
        for i in range(rows):
            for j in range(cols):
                # Calculating curve offset at each pixel
                if curve_type == "sine_curve":
                    offset_x = int(amplitude_levels[k]*cols * math.sin(2 * 3.14 * i / (2*rows)))
                    offset_y = 0
                elif curve_type == "full_sine_curve":
                    offset_x = int(amplitude_levels[k]*cols * math.sin(2 * 3.14 * i / (rows)))
                    offset_y = 0
                elif curve_type == "x4_curve":
                    offset_x = int(amplitude_levels[k]/rows**3 * i**4)
                    offset_y = 0
                elif curve_type == "parabolic_curve":
                    offset_x = int(amplitude_levels[k]/cols * i * (cols-i))
                    offset_y = 0
                elif curve_type == "half_parabolic_curve":
                    offset_x = int(amplitude_levels[k]/cols * i**2)
                    offset_y = 0

                # applying offset
                if i+offset_y >= rows or i+offset_y < 0 or j+offset_x >= cols or j+offset_x <0:
                    if fill_black:
                        img_output[i,j] = img[(i+offset_y)%rows,(j+offset_x)%cols]
                    else:
                        img_output[i,j] = 0
                else:
                    img_output[i,j] = img[(i+offset_y)%rows,(j+offset_x)%cols]

        # Showing and writing image
        if curve_type == "sine_curve":
            # cv2.imshow('Sine Curve '+str(k+1), img_output)
            cv2.imwrite(save_file_name+'_sine_curve_'+str(k+1)+".jpg", img_output)
        elif curve_type == "full_sine_curve":
            # cv2.imshow('Full Sine Curve '+str(k+1), img_output)
            cv2.imwrite(save_file_name+'_full_sine_curve_'+str(k+1)+".jpg", img_output)
        elif curve_type == "x4_curve":
            # cv2.imshow('Half x^4 Curve '+str(k+1), img_output)
            cv2.imwrite(save_file_name+'_half_x4_curve_'+str(k+1)+".jpg", img_output)
        elif curve_type == "parabolic_curve":
            # cv2.imshow('Parabolic Curve '+str(k+1), img_output)
            cv2.imwrite(save_file_name+'_parabolic_curve_'+str(k+1)+".jpg", img_output)
        elif curve_type == "half_parabolic_curve":
            # cv2.imshow('Half Parabolic Curve '+str(k+1), img_output)
            cv2.imwrite(save_file_name+'_half_parabolic_curve_'+str(k+1)+".jpg", img_output)