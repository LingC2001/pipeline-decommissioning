"File is used to convert videos to individual frames of images"

import cv2
import numpy as np

# make sure to change the input_path/ output_path relative to your current dir 

# Create a VideoCapture object and read from input file
# If the input is the camera, pass 0 instead of the video file name
input_path = 'long_side.mp4'
output_path = "video_frames/"
cap = cv2.VideoCapture(input_path)
# Check if camera opened successfully
if (cap.isOpened()== False): 
    print("Error opening video stream or file")

# Read until video is completed
frame_id = 0
img_id = 0
while(cap.isOpened()):
    # Capture frame-by-frame
    ret, frame = cap.read()
    if ret == True:
        if frame_id % 20 == 0:
            cv2.imwrite(output_path + input_path[0:-4]+str(img_id)+".png",frame)
            img_id += 1
        frame_id += 1
    # Break the loop
    else: 
        break

# When everything done, release the video capture object
cap.release()
# Closes all the frames
cv2.destroyAllWindows()