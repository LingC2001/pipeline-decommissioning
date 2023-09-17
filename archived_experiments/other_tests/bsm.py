import cv2 as cv
import numpy as np
import os
from PIL import Image

# Normalise color hue and save them in a new directory - from color to grayscale
def normalise_hue(img_dir,grayscale_dir,width,height):
    for file in sorted(os.listdir(img_dir)):
        filepath = ''
        if file.endswith('.jpg'):
            filepath = img_dir + '/' + file
            image = cv.imread(filepath)
            # im = Image.open(filepath) # Can be many different formats.
            # pix = im.load()
            # print(im.size)
            if image is None:
                break

            #change color to gray by using (r+g+b)*0.33
            for i in range(height):
                for j in range(width):
                    if image[i][j][0] != 0 or image[i][j][1] != 0 or image[i][j][2] != 0:
                        gray_value =sum(image[i][j])*0.33
                        image[i][j] =  [gray_value, gray_value,gray_value] 

            save_path = grayscale_dir + '/'+file
            cv.imwrite(save_path, image)


# apply background model 
def bg_learn(grayscale_dir):
    backSub = cv.createBackgroundSubtractorMOG2()
    for file in sorted(os.listdir(grayscale_dir)):
        filepath = ''
        if file.endswith('.jpg'):
            filepath = grayscale_dir + '/' + file
            frame = cv.imread(filepath)
            if frame is None:
                print(file, None)
                break

            fgMask = backSub.apply(frame, learningRate=1)

    return backSub


# apply foreground extraction
def bg_generate(grayscale_dir,bg_dir,backSub):
    for file in sorted(os.listdir(grayscale_dir)):
        filepath = ''
        if file.endswith('.jpg'):
            filepath = grayscale_dir + '/' + file
            frame = cv.imread(filepath)
            if frame is None:
                print(file, None)
                break

            fgMask = backSub.apply(frame, learningRate=0)

            fgMask = np.expand_dims(fgMask,axis=2)
            fgMask[fgMask<=255/4] = 0
            fgMask[fgMask>255/4] = 1
            fgExtract = frame*fgMask

            save_path = bg_dir + '/'+file
            cv.imwrite(save_path, fgExtract)
            print(save_path,'done')



if __name__ == "__main__":
    # set path according to your local file path
    img_dir = "C:/Users/Hp/OneDrive - Monash University/Desktop/test"
    grayscale_dir = "C:/Users/Hp/OneDrive - Monash University/Desktop/test_gray"
    bg_dir = "C:/Users/Hp/OneDrive - Monash University/Desktop/test_bg"
    # use Pillow -> im.size to check image pixel
    width,height = 3024,3024
    normalise_hue(img_dir,grayscale_dir,width,height)

    backSub = bg_learn(grayscale_dir)
    bg_generate(grayscale_dir,bg_dir,backSub)
