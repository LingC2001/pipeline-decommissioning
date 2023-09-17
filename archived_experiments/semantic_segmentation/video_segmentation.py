import torch
import cv2
import numpy as np
import matplotlib.pyplot as plt
GPU_indx = 0
device = torch.device(GPU_indx if torch.cuda.is_available() else 'cpu')
model = torch.load("Models/underfit.pth")
model2 = torch.load("Models/overfit.pth")

def model_output(img):

    # img = cv2.imread("train_images/p22.jpg")
    # img = cv2.imread("train_images/Screenshot_20220831-182505_Gallery.png")
    img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    img_tensor = torch.unsqueeze(torch.moveaxis(torch.from_numpy(img), 2,0), dim=0).float().to(device)

    output = model(img_tensor)
    output = torch.moveaxis(torch.squeeze(output, dim=0), 0, 2).cpu().detach().numpy()
    carcass, polymer_external_sheath, polymer_fluid_barrier, pressure_armour, tensile_armour = [output[:,:, 0], output[:,:, 1], output[:,:, 2], output[:,:, 3], output[:,:, 4]]

    output2 = model2(img_tensor)
    output2 = torch.moveaxis(torch.squeeze(output2, dim=0), 0, 2).cpu().detach().numpy()
    carcass2, polymer_external_sheath2, polymer_fluid_barrier2, pressure_armour2, tensile_armour2 = [output2[:,:, 0], output2[:,:, 1], output2[:,:, 2], output2[:,:, 3], output2[:,:, 4]]

    return polymer_external_sheath, tensile_armour2

# Create a VideoCapture object and read from input file
# If the input is the camera, pass 0 instead of the video file name
cap = cv2.VideoCapture('input.mp4')
out = cv2.VideoWriter('output.mp4', -1, 20.0, (512,512))


# Check if camera opened successfully
if (cap.isOpened()== False): 
  print("Error opening video stream or file")

# Read until video is completed
while(cap.isOpened()):
  # Capture frame-by-frame
  ret, frame = cap.read()
  if ret == True:

    polymer_external_sheath, tensile_armour = model_output(frame)
    segments = np.zeros_like(frame)
    for i in range(segments.shape[0]):
        for j in range(segments.shape[1]):
            if tensile_armour[i][j] == 1:
                segments[i][j] = [0,0,255]
            elif polymer_external_sheath[i][j] == 1:
                segments[i][j] = [255,0,0]
    out.write(segments)
  # Break the loop
  else: 
    break

# When everything done, release the video capture object
cap.release()
out.release()
# Closes all the frames
cv2.destroyAllWindows()