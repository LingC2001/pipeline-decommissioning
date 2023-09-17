# Import Libraries
import RPi.GPIO as GPIO
import numpy as np
from RpiMotorLib import RpiMotorLib
import matplotlib.pyplot as plt
import time
import serial
import cv2
# stepper hor: DIR/STEP = GPIO05/GPIO06
# stepper left vert: GPIO23/GPIO24
# stepper right vert: GPIO20/GPIO21

############################# DEFINE VARIABLES #############################
# Scan resolution in mm
resolution  = 2
# Declare the clearance/length that can be covered by laser in mm
x_clearance = 220
y_clearance = 220
#############################################################################

# Address of the laser sensor
ADDR = 0x80

# Declaring the serial communication parameters for the sensor
ser = serial.Serial(
    port = '/dev/serial0',
    baudrate = 9600,
    parity = serial.PARITY_NONE,
    stopbits = serial.STOPBITS_ONE,
    bytesize = serial.EIGHTBITS,
    timeout = 1
)

set_addr = [0xFA,0x04,0x01,ADDR,0x81]
ser.write(bytearray(set_addr))
time.sleep(0.05)
while ser.inWaiting()>0:
    x = ser.read()

#char_buff = [ADDR, 0x06, 0x03, 0x77]
#ser.write(bytearray(char_buff))

#define GPIO pins
motor_hor_MS_res = (17,27,22)
motor_hor_step = 5
motor_hor_dir = 6

motor_left_vert_MS_res = (10,9,11)
motor_left_vert_dir = 26
motor_left_vert_step = 16 # False is downwards

motor_right_vert_MS_res = (13,19,26)
motor_right_vert_dir = 24
motor_right_vert_step = 23

# Declare an named instance of class pass GPIO pins numbers
motor_hor = RpiMotorLib.A4988Nema(motor_hor_dir, motor_hor_step, motor_hor_MS_res, "A4988")
motor_left_vert = RpiMotorLib.A4988Nema(motor_left_vert_dir, motor_left_vert_step, motor_left_vert_MS_res, "A4988")
motor_right_vert = RpiMotorLib.A4988Nema(motor_right_vert_dir, motor_right_vert_step, motor_right_vert_MS_res, "A4988")

# Declare origin as the top left corner
x_laser = 0.0
y_laser = 0.0



# Function declaration for checking origin
def origin_check(x,y):
    if x == 0.0 and y == 0.0:
        return True
    else:
        return False

# Function declaration for 'simultaneous' vertical motor movement
def move_vertical(num_of_mm, left_vert, right_vert, directional = False):
    num_of_steps = num_of_mm*25
    steps_remaining = int(num_of_steps)
    while steps_remaining>0:
        left_vert.motor_go(directional, "Full", 1, 0.001, False, 0)
        right_vert.motor_go(directional, "Full", 1, 0.001, False, 0)
        steps_remaining = steps_remaining - 1

def ResizeWithAspectRatio(image, width=None, height=None, inter=cv2.INTER_AREA):
    dim = None
    (h, w) = image.shape[:2]

    if width is None and height is None:
        return image
    if width is None:
        r = height / float(h)
        dim = (int(w * r), height)
    else:
        r = width / float(w)
        dim = (width, int(h * r))

    return cv2.resize(image, dim, interpolation=inter)

def grey_scale(value):
    minimum, maximum = 0, 3
    ratio = (value-minimum)/(maximum-minimum)
    grey = int(ratio*255)
    grey = 255-grey
    print(grey)
    return grey
    
def colour_mapping(value):
    """
    minimum = closest value laser will scan
    maximum = furthest value laser will display colour
    """
    minimum, maximum = 0.05, 0.55
    minimum, maximum = float(minimum), float(maximum)

    if value > maximum:
        return (50, 50, 50)

    ratio = (value-minimum) / (maximum - minimum)
    colour_val = int(ratio*255) 
    return (colour_val, 0, 255-colour_val)

def reset_serial(num=1):
    for i in range(num):
        char_buff = [ADDR, 0x06, 0x02, 0x78]
        ser.write(bytearray(char_buff))
        time.sleep(0.1)
        while ser.inWaiting() > 0:
            data = []
            for i in range(11):
                data.append(ser.read())
            print(data)

            try:
                if chr(data[3][0]) =='E' and chr(data[4][0]) == 'R' and chr(data[5][0]) =='R':
                    print("Out of Range ERROR")
                else:
                    distance = 0.0
                    distance = ((data[3][0])-0x30)*100+((data[4][0])-0x30)*10+((data[5][0])-0x30)*1+((data[7][0])-0x30)*0.1+((data[8][0])-0x30)*0.01+((data[9][0])-0x30)*0.001
                    print("Dumped Distance = " + str(round(distance,3)) + " M from serial")
            except IndexError:
                pass

def soft_clear_serial():
    while ser.inWaiting() > 0:
        data = []
        for i in range(11):
            data.append(ser.read())
        print(data)

        try:
            if chr(data[3][0]) =='E' and chr(data[4][0]) == 'R' and chr(data[5][0]) =='R':
                print("Out of Range ERROR")
            else:
                distance = 0.0
                distance = ((data[3][0])-0x30)*100+((data[4][0])-0x30)*10+((data[5][0])-0x30)*1+((data[7][0])-0x30)*0.1+((data[8][0])-0x30)*0.01+((data[9][0])-0x30)*0.001
                print("Dumped Distance = " + str(round(distance,3)) + " M from serial")
        except IndexError:
            pass

# Read laser data
def read_laser(ADDR, ser, change_ser, dir_of_motion, motor_move):
    # dir_of_motion = True means right
    counter = 0
    reset_counter = 0
    edge_failure = False
    #motor_move = 0
    #ser.reset_input_buffer()
    while counter == 0:
        time.sleep(0.01)
        if reset_counter>300:
            change_ser = True
            char_buff = [ADDR, 0x06, 0x02, 0x78]
            ser.write(bytearray(char_buff))
            time.sleep(0.5)
        if ser.inWaiting() > 0:
            #print(ser.in_waiting)

            data = []
            for i in range (11):
                
                data.append(ser.read())
                #print(len(data))
            try:
                if chr(data[3][0]) == 'E' and chr(data[4][0]) == 'R' and chr(data[5][0])=='R':
                    print("err")
                    distance = 0.0
                    '''
                    motor_hor.motor_go(dir_of_motion, "Full" , 2, .0001, False, 0)
                    motor_move = motor_move + 2
                    reset_serial(1)
                    change_ser = True
                    '''
                    edge_failure = True
                    distance = 0.0
                    counter = 1
                else:   
                    distance = 0.0
                    distance = ((data[3][0])-0x30)*100+((data[4][0])-0x30)*10+((data[5][0])-0x30)*1+((data[7][0])-0x30)*0.1+((data[8][0])-0x30)*0.01+((data[9][0])-0x30)*0.001
                    distance = round(distance,3)

                    if distance>4.0 or distance<0.0:
                        print("Out of Range!")
                        print("Distance = " + str(distance) + " M")
                        reset_serial(3)
                        change_ser = True
                          
                    else:
                        
                        print("Distance = " + str(distance) + " M")
                        reset_counter = 0
                        counter = 1
            except IndexError:
                print("Index Error")
                reset_serial(1)
                change_ser = True 
        reset_counter = reset_counter + 1
    return distance, change_ser, motor_move, ser, edge_failure


# Move the motors slightly initially to remove initdelay
motor_hor.motor_go(False, "Full" , 50, .001, False, .05)
motor_hor.motor_go(True, "Full" , 50, .001, False, 0)

motor_left_vert.motor_go(False, "Full" , 5, .001, False, .05)
motor_left_vert.motor_go(True, "Full" , 5, .001, False, 0)

motor_right_vert.motor_go(False, "Full" , 5, .001, False, .05)
motor_right_vert.motor_go(True, "Full" , 5, .001, False, 0)

# Flags
scan_flag = False
origin_flag = origin_check(x_laser,y_laser)

# Some constants
hor_steps_per_mm = 5
ver_steps_per_mm = 25


num_hor_steps_per_scan = hor_steps_per_mm * resolution
num_vert_steps_per_movement = ver_steps_per_mm * resolution

total_horizontal_steps = x_clearance*hor_steps_per_mm
total_vertical_steps = y_clearance*ver_steps_per_mm

total_hor_dist = x_clearance
hor_dist_per_step = total_hor_dist/total_horizontal_steps
hor_dist_per_scan = hor_dist_per_step*num_hor_steps_per_scan

total_hor_scans = total_horizontal_steps/num_hor_steps_per_scan
total_ver_scans = total_vertical_steps/num_vert_steps_per_movement

# Predeclaring a matrix for recording the values
record_height = np.zeros((int(total_ver_scans),int(total_hor_scans)+1))
'''
plt.imshow(record_height,cmap='autumn', interpolation='nearest', extent=[0,total_hor_dist,220,0])
plt.title("2-D Heatmap")
plt.show()
'''


height_map = np.zeros((int(total_ver_scans),int(total_hor_scans)+1,3), dtype = np.uint8)
print(record_height.shape)
                


# Scanning process
while not scan_flag:
    vert_steps_remaining = total_vertical_steps
    i = 0
    j = 0
    flag = True
    change_ser = True
    motor_move = 0
    last_measurement = 0
    
    


    while vert_steps_remaining>0.0:
        if change_ser:
            change_ser = False
            char_buff = [ADDR, 0x06, 0x03, 0x77]
            ser.write(bytearray(char_buff))
        hor_steps_remaining = total_horizontal_steps
        '''

        Code for performing one scan at initial position of each level

        '''
        motor_hor.motor_go(not flag, "Full" , 2*num_hor_steps_per_scan, .001, False, 0)
        soft_clear_serial()
        record_height[i,j], change_ser, motor_move, ser, edge_failure = read_laser(ADDR,ser, change_ser, flag, motor_move)
        if edge_failure:
            record_height[i,j] = last_measurement
        last_measurement = record_height[i,j]
        height_map[i][j] = colour_mapping(record_height[i,j])
        cv2.imshow("height heat-map", ResizeWithAspectRatio(height_map, width=768))
        cv2.waitKey(1)
        
        while hor_steps_remaining>0.0:
            if change_ser:
                change_ser = False
                char_buff = [ADDR, 0x06, 0x03, 0x77]
                ser.write(bytearray(char_buff))

            if motor_move>0:
                corrector = motor_move/num_hor_steps_per_scan
                num_corrector_scans = np.floor(corrector)
                if num_corrector_scans>0:
                    fixer = int(num_corrector_scans*num_hor_steps_per_scan)
                    motor_hor.motor_go(flag, "Full" , fixer, .001, False, 0)
                    motor_move = motor_move - fixer
                    
            if flag == True:
                motor_hor.motor_go(False, "Full" , num_hor_steps_per_scan-motor_move, .001, False, 0)
                j = j + 1
                x_laser = x_laser + 1
            else:
                motor_hor.motor_go(True, "Full" , num_hor_steps_per_scan-motor_move, .001, False, 0)
                j = j - 1
                x_laser = x_laser - 1
            
            motor_move = 0

            '''
            
            Code for performing one scan
            
            '''
            soft_clear_serial()
            record_height[i,j], change_ser, motor_move, ser, edge_failure = read_laser(ADDR,ser, change_ser, flag, motor_move)
            if edge_failure:
                record_height[i,j] = last_measurement
            last_measurement = record_height[i,j]
            hor_steps_remaining = hor_steps_remaining - num_hor_steps_per_scan
            
            height_map[i][j] = colour_mapping(record_height[i,j])
            cv2.imshow("height heat-map", ResizeWithAspectRatio(height_map, width=768))
            cv2.waitKey(1)
        '''
        if motor_move>0:
            motor_hor.motor_go(not flag, "Full" , motor_move, .001, False, 0)
            motor_move = 0    
        '''
        move_vertical(resolution, motor_left_vert, motor_right_vert, False)

        i = i+1
        y_laser = y_laser - 1
        flag = not flag
        vert_steps_remaining = vert_steps_remaining - num_vert_steps_per_movement

    scan_flag = True

cv2.imshow("height heat-map", ResizeWithAspectRatio(height_map, width=768))
cv2.waitKey(1)

while not origin_check(x_laser,y_laser):
    if x_laser != 0:
        hor_steps_to_origin = int(x_laser/hor_dist_per_step)
        if x_laser>0:
            motor_hor.motor_go(True, "Full" , hor_steps_to_origin, .001, False, 0)
            x_laser = x_laser - hor_steps_to_origin*hor_dist_per_step
        else:
            motor_hor.motor_go(False, "Full" , hor_steps_to_origin, .001, False, 0)
            x_laser = x_laser + hor_steps_to_origin*hor_dist_per_step
        
    if y_laser != 0:
        vert_steps_to_origin = int(np.round(np.abs(y_laser*25)))
        if y_laser<0:
            while np.abs(vert_steps_to_origin>0):
                move_vertical(resolution, motor_left_vert, motor_right_vert, True)
                y_laser = y_laser + 1
                vert_steps_to_origin = vert_steps_to_origin - num_vert_steps_per_movement

print("Saving scan")
cv2.imwrite("scan.png", height_map)

# save height_record into ttxt file
file = open("scan.txt", "w")
for i in range(record_height.shape[0]):
    if i!= 0 and i!= record_height.shape[0]:
        file.write("\n")
    for j in range(record_height.shape[1]):
        file.write(str(round(record_height[i][j], 3)))
        file.write(" ")
file.close()




print("shutting down laser")           
shut_down = [ADDR,0x04,0x02,0x7A]
ser.write(bytearray(shut_down))









