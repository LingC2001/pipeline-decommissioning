# PIPELINE-DECOMMISIONING

## The purpose of the project
This repo contains work done performed together with Monash Software and Civil Engineering team to:

1. Digitalisation of Pipe Images to Determine Workability
    - Detect damage and deformation from pipe images
    - Determine suitability in undergoing separation procedures

2. Image Data Interpretation for Process Optimisation
    - Derive important information from pipe images
    - Determine important parameters of the separation process

3. Information Translation for Automation Systems
    - Design interface to relay optimised parameters to the autonomous separation system
  
4. Integration with Kuka robotic arm
   - Perform cutting/proof of concept with robotic arm path

## General approach
1. check pipe intergrity/damage **pipe_workability** - ResNet50 Pre-trained Model
    - **model_learning** -> check if the model learns the correct feature of pipes
2. identify the pipe's width and length - Thresholding & Canny Edge Detection
3. detect, identify and segment the different layers of pipes - Semantic Segmentation **semantic_segmentation**
4. simulate the whole process by using RoboDK **robotic_sumulation**

## How to set up
1. Clone the repo from [GitHub](https://github.com/LingC2001/pipeline-decommissioning).
2. Follow [this](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/) link to create a virtual environment.

## Contact people
Chunyang Chen (Faculty of IT) Chunyang.Chen@monash.edu

Zijin Chen (Faculty of IT) zche0135@student.monash.edu

Victor Chang (Faculty of Engineering) victor.chang@monash.edu

Jenny Zhou (Faculty of Engineering) jenny.zhou@monash.edu

Ling Chen (Faculty of Engineering) lche0147@student.monash.edu

## Useful rescourses
Image Augmentation [Keras](https://machinelearningmastery.com/how-to-configure-image-data-augmentation-when-training-deep-learning-neural-networks/)

Image Augmentation [PyTorch](https://pytorch.org/vision/main/transforms.html)

Canny Edge Detection [YouTube](https://www.youtube.com/watch?v=gmrbZOpPeno)

Canny Edge Detection [Web Page](canny-edge-detection-step-by-step-in-python-computer-vision-b49c3a2d8123)

Image Distortion [Web Page](https://subscription.packtpub.com/book/application-development/9781785283932/1/ch01lvl1sec16/image-warping)
