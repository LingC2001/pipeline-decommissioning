import labelme

from PIL import Image
import os
import os.path as osp 
import glob
import cv2
import json

from labelme import utils
import argparse
import base64

def labelme2images(input_dir, output_dir, save_img=False):
    """
    new_size_width, new_size_height = new_size
    """
    
    print("Generating dataset")
    
    filenames = glob.glob(osp.join(input_dir, "*.json"))

    for filename in filenames:

        data = json.load(open(filename))
        imageData = data.get("imageData")

        # base name
        base = osp.splitext(osp.basename(filename))[0]

        if not imageData:
            imagePath = os.path.join(os.path.dirname(json_file), data["imagePath"])
            with open(imagePath, "rb") as f:
                imageData = f.read()
                imageData = base64.b64encode(imageData).decode("utf-8")
        img = utils.img_b64_to_arr(imageData)


        label_name_to_value = {"_background_": 0}
        for shape in sorted(data["shapes"], key=lambda x: x["label"]):
            label_name = shape["label"]
            if label_name in label_name_to_value:
                label_value = label_name_to_value[label_name]
            else:
                label_value = len(label_name_to_value)
                label_name_to_value[label_name] = label_value
        lbl, _ = utils.shapes_to_label(
            img.shape, data["shapes"], label_name_to_value
        )
        
        utils.lblsave(osp.join(output_dir, base), lbl)

input_dir = "annotated"
output_dir = "labels"

labelme2images(input_dir, output_dir, save_img=True)