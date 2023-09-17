import cv2
import numpy as np

cap = cv2.VideoCapture('demo_vid.mp4')

# Check if camera opened successfully
if (cap.isOpened()== False): 
  print("Error opening video stream or file")

out = cv2.VideoWriter('demo_vid_resized.mp4', -1, 20.0, (512,512))

while True:
    ret, image = cap.read()
    if ret == True:
        # Crop image to square
        crop_size = min(image.shape[0], image.shape[1])
        start_row = (image.shape[0]-crop_size)//2
        end_row = (image.shape[0]-crop_size)//2 + crop_size
        start_col = (image.shape[1]-crop_size)//2
        end_col = (image.shape[1]-crop_size)//2 + crop_size
        image = image[start_row:end_row, start_col:end_col]
        
        # Rescale image
        image = cv2.resize(image, (512,512),fx=0,fy=0, interpolation=cv2.INTER_CUBIC)
        out.write(image)
    else:
        break
    
cap.release()
out.release()

cv2.destroyAllWindows()


