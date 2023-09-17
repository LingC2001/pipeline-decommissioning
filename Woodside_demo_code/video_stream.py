from flask import Flask, render_template, Response
import cv2
from long_camera_generate_frames import generate_frames


app = Flask(__name__)
camera = cv2.VideoCapture(0)

@app.route('/')

def index():
    return render_template('index.html')

@app.route('/video')

def video():
    return Response(generate_frames(camera), mimetype='multipart/x-mixed-replace; boundary=frame')

if __name__ == "__main__":
    app.run(host="0.0.0.0")