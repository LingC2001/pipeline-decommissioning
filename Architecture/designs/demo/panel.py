import streamlit as st
from demo_server import Server
from demo_client import Client
import subprocess
import time
import random as rnd

IMG_DIR = 'img/'
ROBOT_SIM_DIR = '../../../sim/Matlab Example/simp_Robot_Model_Main_test.m'
_, serverCol, _ = st.columns(3)

with serverCol:
    st.markdown("<h1 style='text-align: center; '>Server</h1>", unsafe_allow_html=True)
    st.image(IMG_DIR + "servers.png", width=240)


server_open_ports = st.multiselect(
    'Specify the listening ports:',
    list(range(65536)))

if st.button('Run Server'):
    subprocess.Popen(['python', 'demo_server.py', str(server_open_ports)], close_fds=True)
    st.write('Server now listening to ports:', str(server_open_ports))
# ----------------------------------------------------------------------------------------
st.markdown("""---""")

r2c1, r2c2, r2c3 = st.columns(3)

with r2c1:
    # get port
    # images directory

    st.markdown("<h2 style='text-align: center; '>Camera</h2>", unsafe_allow_html=True)
    st.image(IMG_DIR + "photo-camera.png")

    camera_port = st.number_input('Specify a port',
     min_value=0,
     max_value=65536,
     key='cp')

    if st.button('Run Camera'):
        subprocess.Popen(['python', 'demo_client.py', str(camera_port), 'camera'], close_fds=True)
        st.write('Camera sending signals to port no.', str(camera_port))

with r2c2:
    # get port

    st.markdown("<h2 style='text-align: center; '>IR Sensor</h2>", unsafe_allow_html=True)
    st.image(IMG_DIR + "sensor.png")

    IR_scanner_port = st.number_input('Specify a port',
     min_value=0,
     max_value=65536,
     key='irp')

    if st.button('Run IR Sensor'):
        subprocess.Popen(['python', 'demo_client.py', str(IR_scanner_port), 'IR'], close_fds=True)
        st.write('IR Sensor sending signals to port no.', str(IR_scanner_port))

with r2c3:
    # get port
    st.markdown("<h2 style='text-align: center; '>Dust Detector</h2>", unsafe_allow_html=True)
    st.image(IMG_DIR + "smoke-detector.png")

    sensor_port = st.number_input('Specify a port',
     min_value=0,
     max_value=65536,
     key='sp')

    # title__ = st.text_input('Movie title', 'Life of Brian')
    # st.write('The current movie title is', title)

    if st.button('Run Dust Detector'):
        subprocess.Popen(['python', 'demo_client.py', str(sensor_port), 'dust'], close_fds=True)
        st.write('Dust Detector sending signals to port no.', str(sensor_port))

st.markdown("""---""")

r3c1, r3c2, r3c3 = st.columns(3)

with r3c1:
    st.markdown("<h2 style='text-align: center; '>Server</h2>", unsafe_allow_html=True)
    st.image(IMG_DIR + "servers.png", width=200)
        

with r3c2:
    st.markdown("<h3 style='text-align: center; '>exchange</h3>", unsafe_allow_html=True)
    st.image(IMG_DIR + "exchange.png", width=200)

    robot_port = st.number_input('Specify a port for robot status updates',
     min_value=0,
     max_value=65536,
     key='rpp')

    if st.button('Send Commands'):
        ROBOT_PORT = 11111
        # Open a port on robot'
        subprocess.Popen(['python', 'demo_server.py', str([ROBOT_PORT])], close_fds=True)
        # Send commands from the central server to robot
        subprocess.Popen(['python', 'demo_client.py', str(ROBOT_PORT), 'coords'], close_fds=True)
        # Send status updates from the robot
        subprocess.Popen(['python', 'demo_client.py', str(robot_port), 'robot_st'], close_fds=True)
        # Move the robot
        subprocess.Popen(['matlab', '-nosplash','-nodesktop' , '-r', f"""run('{ROBOT_SIM_DIR}');exit;"""], close_fds=True)
        # Delay
        time.sleep(rnd.randint(2,4))
        # Send status updates from the robot
        subprocess.Popen(['python', 'demo_client.py', str(robot_port), 'robot_end'], close_fds=True)
    
with r3c3:
    st.markdown("<h2 style='text-align: center; '>Robot</h2>", unsafe_allow_html=True)
    st.image(IMG_DIR + "robotic-arm.png", width=220)
