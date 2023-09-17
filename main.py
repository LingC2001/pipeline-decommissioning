"""
This is the main file used to test all the image analysis methods on images together.
"""

import cv2
import numpy as np
from data_interpretation.object_measurement import object_measurement
from pipe_workability.damage_detection import damage_detection
from pipe_workability.resize_for_model import resize_for_model
from model_learning.visualise_prediction import visualise_prediction
from data_interpretation.layer_segmentation import layer_segmentation


def main():
    image_long = cv2.imread("relevant_images/p8.jpg")
    image_cross = cv2.imread("relevant_images/51.JPG")
    cam_height = 2.2
    real_frame_width = 1.08
    # image_long = cv2.imread("robotic_simulation/images/side.jpg")
    # image_cross = cv2.imread("robotic_simulation/images/section.jpg")
    # # real_frame_width = (1100/5)*512*11.276/1000/1000
    # cam_height = 1.1
    # real_frame_width = 2*cam_height*np.tan(30*np.pi/180)
    

    #Showing original images
    cv2.imshow("Original external image", cv2.resize(image_long, (768,768), interpolation=cv2.INTER_CUBIC))
    cv2.imshow("Original cross section image", cv2.resize(image_cross, (768,768), interpolation=cv2.INTER_CUBIC))


    # TODO: Determine pipe workability using both images
    damage_type, external_prob, internal_prob = damage_detection(resize_for_model(image_long), "pipe_workability/Models/ResNet50_External.pt", resize_for_model(image_cross), "pipe_workability/Models/ResNet50_Internal.pt")
    external_CAM = visualise_prediction(resize_for_model(image_long), "pipe_workability/Models/ResNet50_External.pt")
    internal_CAM = visualise_prediction(resize_for_model(image_cross), "pipe_workability/Models/ResNet50_Internal.pt")
    
    # Showing results
    print("Detected Damage Type: " + damage_type)
    print("Externally damaged probability: " + str(external_prob))
    print("Internally damaged probability: " + str(internal_prob))
    cv2.imshow("external CAM", cv2.cvtColor(external_CAM, cv2.COLOR_RGB2BGR))
    cv2.imshow("internal CAM", cv2.cvtColor(internal_CAM, cv2.COLOR_RGB2BGR))


    # TODO: Determine pipe dimensions
    dimensions, obj_length, obj_width = object_measurement(image_long, real_frame_width, cam_height, background=False, very_long=False, is_horizontal = None)

    # Showing results
    print('The length of the pipe is: %.2f mm'%(obj_length))
    print('The diameter of the pipe is: %.2f mm'%(obj_width))
    cv2.imshow("dimensions", dimensions)


    # TODO: Layer segmentation
    cam_dis = 0.2
    real_frame_dis = 2*cam_dis*np.tan(30*np.pi/180)
    image_layers , radii = layer_segmentation(image_cross, real_frame_dis)

    # Showing results
    cv2.imshow("Pipe Layers", image_layers)
    print("Layers Radii: " + str(radii))

    
    cv2.waitKey(0)
    cv2.destroyAllWindows()


if __name__ == "__main__":
    main()