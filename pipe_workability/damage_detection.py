"""
This script called trained Resnet model to determine pipe workability by detecting damage

"""

import torch
from torchvision.models import resnet50
import numpy as np
import torchvision.transforms as T
import torch.nn as nn

def softmax(x):
    """Compute softmax values for each sets of scores in x."""
    e_x = np.exp(x - np.max(x))
    return e_x / e_x.sum() 

def damage_detection_predictions(model_path, image):
    """
    Function to use a trained Resnet50 Model to detect whether a pipe is damaged or not

    :param model_path: The file path to the model
    :param image: The openCV image to classify
    :returns predicted_label: 0=undamaged, 1=damaged
    :returns output: The output for the predictions
    
    """
    #Set device to GPU_indx if GPU is avaliable
    GPU_indx = 0
    device = torch.device(GPU_indx if torch.cuda.is_available() else 'cpu')

    # Loading trained model
    model = resnet50(pretrained=True)
    num_ftrs = model.fc.in_features
    model.fc = nn.Linear(num_ftrs, 2)
    model = model.to(device)

    check_point = torch.load(model_path)
    model.load_state_dict(check_point['model_state_dict'])

    # Get model prediction of image
    img_tensor = torch.unsqueeze(torch.moveaxis(torch.from_numpy(image), 2,0), dim=0).to(device)
    output = model(img_tensor).cpu().detach().numpy()
    predicted_label = np.argmax(output, 1)[0]

    return predicted_label, output

def damage_detection(external_image, external_model, internal_image, internal_model):
    """
    Function: given an outer image of a pipe and a cross sectional image of a pipe, classifies the
    type of damage the pipe has

    :param external_image: image of the exterior of the pipe
    :param external_model: model file path for external damage
    :param internal_image: image of the cross sectional view of the pipe
    :param internal_model: model file path for internal damage
    :returns damage_type: Type of damage
    :returns external_prob: The probability of having external damage vs no external damage
    :returns internal_prob: The probability of having internal damage vs no internal damage
    
    
    """

    external_output, external_prob = damage_detection_predictions(external_model, external_image)
    internal_output, internal_prob = damage_detection_predictions(internal_model, internal_image)

    print(external_prob, internal_prob)
    print(external_output, internal_output)

    # conveerting to probability using softmax
    external_prob = softmax(np.squeeze(external_prob,0))
    internal_prob = softmax(np.squeeze(internal_prob,0))

    if external_output == 0 and internal_output == 0: # Not damaged
        damage_type = "Not Damaged"
    elif external_output == 0 and internal_output == 1: # Only internally damaged
        damage_type = "Only Internally Damaged"
    elif external_output == 1 and internal_output == 0: # Only externally damaged
        damage_type = "Only Externally Damaged"
    else: # Both internally and externally damaged 
        damage_type = "Both Internally and Externally Damaged"
        
    return damage_type, external_prob, internal_prob




