import numpy as np
import cv2
from matplotlib import pyplot as plt
from pathlib import Path
import os
import shutil
import skvideo.io
import time

def rotate_and_crop(img, img_ref, template, cropped_width, cropped_height):
    """
    Rotate img to match img_ref, then match template (a consistent, agreed upon anchor)
    in the rotated img and crop a width x height image based on template location.

    Parameters
    ----------
    img : ndarray(height, width) or ndarray(height, width, 3)
        Image to be processed
    img_ref : ndarray(height, width) or ndarray(height, width, 3)
        Reference image for alignment. Should be consistent across whole experiment and image shape should match subject img
    template : ndarray(height, width)
        Anchor template to set as reference point. Should be consistent across whole experiment. Template shape should be smaller than img shape
    cropped_width : int
        Final cropped image width
    cropped_height : int
        Final cropped image height

    Returns
    -------
    final_img : ndarray(height, width)
        Rotated and cropped img
    2d_rotation : ndarray(2,3)
        Calculated transformation matrix for 2d rotation. For use in cv2.warpAffine(img, 2d_rotation, (width,height)).
    top_left : (x,y)
        x,y coordinates of top left corner for cropping
    bot_right : (x,y)
        x,y coordinates of bot right corner for copping
    """
    # Image shape
    height, width = img.shape[:2]

    assert img.shape == img_ref.shape # img and ref_img shape should match

    # Extract orb key points in both images
    orb_detector = cv2.ORB_create(4000)
    kp1, d1 = orb_detector.detectAndCompute(img, None)
    kp2, d2 = orb_detector.detectAndCompute(img_ref, None)

    # Find key points that match in both images
    matcher = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck = True)

    # Match the two sets of descriptors.
    matches = list(matcher.match(d1, d2))

    # Sort matches on the basis of their Hamming distance.
    matches.sort(key = lambda x: x.distance)

    # Take the top 90 % matches forward.
    matches = matches[:int(len(matches)*90)]
    no_of_matches = len(matches)

    # Define empty matrices of shape no_of_matches * 2.
    p1 = np.zeros((no_of_matches, 2))
    p2 = np.zeros((no_of_matches, 2))

    # Match point coordinate matrix
    for i in range(len(matches)):
        p1[i, :] = kp1[matches[i].queryIdx].pt
        p2[i, :] = kp2[matches[i].trainIdx].pt

    # Find the homography matrix.
    homography, mask = cv2.estimateAffinePartial2D(p1, p2, cv2.LMEDS)

    # Remove the scaling factor, we are only interested in rotation.
    s = np.sqrt(homography[0,0]**2 + homography[0,1]**2)
    homography = homography/s

    # Remove translation, we are only interested in rotation.
    homography[0,2]=0
    homography[1,2]=0

    # Calculate rotation angle
    rot_angle = np.arccos(homography[0,0])

    # Perform rotation
    transformed_img = cv2.warpAffine(img,
                        homography, (width, height))

    # Apply template Matching
    res = cv2.matchTemplate(transformed_img[:60,:60],template,cv2.TM_CCOEFF_NORMED)
    min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(res)

    top_left = max_loc
    print(top_left)
    bottom_right = (top_left[0] + cropped_width, top_left[1] + cropped_height)

    final_img = transformed_img[top_left[1]:bottom_right[1],top_left[0]:bottom_right[0]]

    return (final_img,homography,top_left,bottom_right)

def preprocess_raw_video(room_name, session_name, data_path, config):
    """
    Preprocessing pipeline for video data

    Parameters
    ----------
    room_name : string
        'b151', 'b149f', or 'b149b'
    session_name : string
        name of session (e.g. 211105_29837)
    config
        config dict containing general experiment information and settings

    Returns
    -------

    I/O
    ---
    Processed videos saved in session_path with raw/ replaced with processed/
    """
    root = config['project_root']

    assert room_name in ('b151', 'b149f', 'b149b')

    # Set paths
    data_path = os.path.join(data_path,room_name)
    cameras_path = os.path.join(data_path, 'cameras')

    # Camera data paths for all cameras in the specified room
    camera_paths = {}
    for dir in os.listdir(cameras_path)[:]:
        for f in os.listdir(os.path.join(cameras_path, dir))[:]:
            fpath = os.path.join(cameras_path,dir,f)
            camera_paths[dir] = fpath

    # Preprocess each camera based on settings in config
    for camera in camera_paths.keys():
        if(camera in config[room_name]['cameras'].keys()):
            # ---------------
            # Rotate and Crop
            # ---------------
            # Load configuration
            cam_settings = config[room_name]['cameras'][camera]
            fps = cam_settings['fps']
            anchor_path = cam_settings['anchor']
            reference_path = cam_settings['reference']
            crop_w = cam_settings['crop_width']
            crop_h = cam_settings['crop_height']

            # Get anchor and reference image
            anchor_img = cv2.imread(os.path.join(root,anchor_path))
            ref_img = cv2.imread(os.path.join(root,reference_path))

            # Read video
            cap = cv2.VideoCapture(camera_paths[camera])
            num_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
            cap.set(cv2.CAP_PROP_POS_FRAMES, 50) # The first few frames have more artifacts, skip them here.
            ret, frame = cap.read()
            # Compute rigid homography on first frame. Use this homography transformation for entire video
            # 5000x5000 is a place holder to retain as much of the frame as possible after cropping.
            transformed_img, rotation_transform, top_left, _ = rotate_and_crop(frame, ref_img, anchor_img, 5000, 5000)
            print("Cropping reference (B151 bottom):",top_left)
            cap.set(cv2.CAP_PROP_POS_FRAMES, 0) # Reset read pointer to frame 0

            out_dir = os.path.dirname(camera_paths['bottom'].replace('raw','processed'))
            Path(out_dir).mkdir(parents=True, exist_ok=True)

            video_out = {}

            for (width,height) in zip(crop_w,crop_h):
                width = min(width,transformed_img.shape[1]) # Maximum width
                height = min(height,transformed_img.shape[0]) # Maximum height
                out_name = os.path.join(out_dir,'{}x{}.mp4'.format(width, height))
                #use the h.264 codec
                #set the constant rate factor to 0, which is lossless
                #the slower the better compression, in princple})
                video_out[(width,height)] = skvideo.io.FFmpegWriter(out_name, inputdict={'-r':str(fps),'-s':'{}x{}'.format(width,height)},outputdict={'-r':str(fps),'-c:v': 'libx264','-crf': '0','-preset':'ultrafast', '-s':'{}x{}'.format(width,height), '-pix_fmt': 'rgb24'})

            start_time = time.time()
            for i in range(num_frames):
                if(i%int(num_frames/20) == 0 and i!=0):
                    print("Processed {}/{} frames | {} | FPS: {}".format(i,num_frames,time.strftime("%H:%M:%S"),i/(time.time()-start_time)))

                ret, frame = cap.read()

                for (width,height) in zip(crop_w,crop_h):
                    # Rotate and Crop frame
                    width = min(width,transformed_img.shape[1]) # Maximum width
                    height = min(height,transformed_img.shape[0]) # Maximum height
                    bottom_right = (top_left[0] + width, top_left[1] + height)

                    # Perform rotation
                    rotated_frame = cv2.warpAffine(frame, rotation_transform, frame.shape[:2][::-1])
                    # Crop based on anchor position and configured crop size
                    cropped_frame = rotated_frame[top_left[1]:bottom_right[1], top_left[0]:bottom_right[0]]
                    # Write to corresponding video file
                    video_out[(width,height)].writeFrame(cropped_frame)

            # Cleanup
            for (width,height) in zip(crop_w,crop_h):
                width = min(width,transformed_img.shape[1]) # Maximum width
                height = min(height,transformed_img.shape[0]) # Maximum height
                video_out[(width,height)].close()
            cap.release()
            cv2.destroyAllWindows()
        else: # Skip cameras that don't have configurations set
            print('{} {} camera has no configuration. Please set camera configuration in config.json'.format(room_name, camera))
            Path(os.path.dirname(camera_paths[camera].replace('raw','processed'))).mkdir(parents=True, exist_ok=True)
            shutil.copyfile(camera_paths[camera], camera_paths[camera].replace('raw','processed'))
