
pipe - v1 2022-09-19 4:25am
==============================

This dataset was exported via roboflow.com on September 18, 2022 at 6:33 PM GMT

Roboflow is an end-to-end computer vision platform that helps you
* collaborate with your team on computer vision projects
* collect & organize images
* understand unstructured image data
* annotate, and create datasets
* export, train, and deploy computer vision models
* use active learning to improve your dataset over time

It includes 79 images.
Pipe-layers are annotated in COCO format.

The following pre-processing was applied to each image:
* Auto-orientation of pixel data (with EXIF-orientation stripping)
* Resize to 512x512 (Stretch)

The following augmentation was applied to create 3 versions of each source image:
* 50% probability of horizontal flip
* 50% probability of vertical flip
* Equal probability of one of the following 90-degree rotations: none, clockwise, counter-clockwise, upside-down
* Randomly crop between 0 and 29 percent of the image
* Random rotation of between -45 and +45 degrees
* Random shear of between -23° to +23° horizontally and -23° to +23° vertically
* Random brigthness adjustment of between -25 and +25 percent
* Random exposure adjustment of between -6 and +6 percent
* Random Gaussian blur of between 0 and 1 pixels
* Salt and pepper noise was applied to 2 percent of pixels


