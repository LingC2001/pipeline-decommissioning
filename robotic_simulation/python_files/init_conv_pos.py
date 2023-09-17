# Type help("robodk.robolink") or help("robodk.robomath") for more information
# Press F5 to run the script
# Documentation: https://robodk.com/doc/en/RoboDK-API.html
# Reference:     https://robodk.com/doc/en/PythonAPI/robodk.html
# Note: It is not required to keep a copy of this file, your Python script is saved with your RDK project

# You can also use the new version of the API:
from robolink import *
from robodk import *
RDK = Robolink()

def init_conv_pos():
    MECHANISM_NAME = "Conveyor"
    mechanism = RDK.Item(MECHANISM_NAME, itemtype=ITEM_TYPE_ROBOT)

    if mechanism.Valid():
        mechanism.setJoints([1000])