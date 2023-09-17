from robolink import *
from robodk import *
RDK = Robolink()

def move_conv():
    MECHANISM_NAME = "Conveyor"
    mechanism = RDK.Item(MECHANISM_NAME, itemtype=ITEM_TYPE_ROBOT)
    PART_TRAVEL_MM = 900

    if mechanism.Valid():
        mechanism.MoveJ(mechanism.Joints()+ PART_TRAVEL_MM)

