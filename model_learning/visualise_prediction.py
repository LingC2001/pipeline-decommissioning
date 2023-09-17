from pytorch_grad_cam import GradCAM
from pytorch_grad_cam.utils.image import show_cam_on_image
from torchvision.models import resnet50
import torch
import torch.nn as nn
import cv2

def visualise_prediction(image, model_file):
    """
    Function to visualise the model using GradCAM
    :param image: The image to pass to the model
    :param model_file: The file path of the model 
    :returns img_cam: The class activation mapping of the model
    
    """
    #Set device to GPU_indx if GPU is avaliable
    GPU_indx = 0
    device = torch.device(GPU_indx if torch.cuda.is_available() else 'cpu')

    model_external = resnet50(pretrained=True)
    num_ftrs = model_external.fc.in_features
    model_external.fc = nn.Linear(num_ftrs, 2)
    model_external = model_external.to(device)

    check_point_external = torch.load(model_file)
    model_external.load_state_dict(check_point_external['model_state_dict'])

    target_layers = [model_external.layer4[-1]]

    # Construct the CAM object once, and then re-use it on many images:
    cam = GradCAM(model=model_external, target_layers=target_layers, use_cuda=True)
    
    img_tensor = torch.unsqueeze(torch.moveaxis(torch.from_numpy(image), 2,0), dim=0).to(device)
    targets = None

    # You can also pass aug_smooth=True and eigen_smooth=True, to apply smoothing.
    grayscale_cam = cam(input_tensor=img_tensor, targets=targets, eigen_smooth=True, aug_smooth=True)

    # In this example grayscale_cam has only one image in the batch:
    grayscale_cam = grayscale_cam[0, :]
    img_cam = show_cam_on_image(cv2.cvtColor(image, cv2.COLOR_BGR2RGB), grayscale_cam, use_rgb=True)

    return cv2.resize(img_cam, (768, 768), interpolation=cv2.INTER_CUBIC)
