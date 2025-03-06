import cv2
import numpy as np

# Path to your AVI file
avi_path = "SWI_LifeTime_data/LifeTime_best_det.avi"
npy_output = "photon_decay_data.npy"

# Open video file
cap = cv2.VideoCapture(avi_path)

# Check if the file opened successfully
if not cap.isOpened():
    print("Error: Could not open AVI file.")
    exit()

# Read frames and store them in a list
frames = []
while True:
    ret, frame = cap.read()
    if not ret:
        break  # Stop when there are no more frames
    gray_frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)  # Convert to grayscale
    frames.append(gray_frame)

# Convert list to NumPy array (Frames, Height, Width)
frames_array = np.array(frames)

# Save as .npy file
np.save(npy_output, frames_array)

# Release resources
cap.release()

print(f"Saved {len(frames)} frames to {npy_output}")
