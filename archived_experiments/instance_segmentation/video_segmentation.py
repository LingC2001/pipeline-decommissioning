import torch
import cv2
import numpy as np
import matplotlib.pyplot as plt
import pickle
import sys
import os 

sys.path.append("detectron2-main/")
from detectron2.engine import DefaultTrainer, DefaultPredictor

GPU_indx = 0
device = torch.device(GPU_indx if torch.cuda.is_available() else 'cpu')

cfg_save_path = "IS_cfg.pickle"
with open(cfg_save_path, 'rb') as f:
    cfg = pickle.load(f)

cfg.MODEL.WEIGHTS = os.path.join(cfg.OUTPUT_DIR, "model_final.pth")
cfg.MODEL.ROI_HEADS.SCORE_THRESH_TEST = 0.4

model = DefaultPredictor(cfg)

def model_output(img):

    # img = cv2.imread("train_images/p22.jpg")
    # img = cv2.imread("train_images/Screenshot_20220831-182505_Gallery.png")
    img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    outputs = model(img)
    labels = outputs["instances"].pred_classes.cpu().numpy()
    masks = outputs["instances"].pred_masks.cpu().numpy()

    return masks, labels

# Create a VideoCapture object and read from input file
# If the input is the camera, pass 0 instead of the video file name
cap = cv2.VideoCapture('section_side.mp4')
out = cv2.VideoWriter('output_side.mp4', -1, 30.0, (512,512))


# Check if camera opened successfully
if (cap.isOpened()== False): 
  print("Error opening video stream or file")

label_map = {
    1: (255, 0, 0),
    2: (0, 255, 0),
    3: (0, 0, 255),
    4: (255, 255, 0),
    5: (0, 255, 255)
}

# Read until video is completed
while(cap.isOpened()):
  # Capture frame-by-frame
  ret, frame = cap.read()
  if ret == True:

    masks, labels = model_output(frame)
    segments = np.zeros_like(frame)
    for k in range(masks.shape[0]):
      for i in range(segments.shape[0]):
          for j in range(segments.shape[1]):
              if masks[k][i][j] != False:
                  segments[i][j] = label_map[labels[k]]
    out.write(segments)
  # Break the loop
  else: 
    break

# When everything done, release the video capture object
cap.release()
out.release()
# Closes all the frames
cv2.destroyAllWindows()