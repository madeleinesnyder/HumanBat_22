import deeplabcut
import os
import numpy as np

print(os.listdir())
print(deeplabcut.__version__)
vid_path = os.environ.get('VIDEO_PATH')

# Project config path (point to desired trained net project)
project_config_path = 'HumanBat/HumanBat_DLC/dlc_projects/B151_220104_IMPLANT_LIGHTS_BOTTOM-Madeleine-2022-01-13/config.yaml'

# Video file
vid_name = [fname for fname in os.listdir(vid_path) if fname.endswith('.mp4')][0]
vid_path = os.path.join(vid_path, vid_name)
print(vid_path)

print("\nAnalyzing video (Inference)\n")
deeplabcut.analyze_videos(project_config_path,
                          [vid_path],
                          videotype='.mp4',
                          dynamic=(True,.5,10))
                          #use_shelve=True) use_shelve is not accepted as argument for analyze_videos for some reason...

print("\nCreating video with detections\n")
deeplabcut.create_video_with_all_detections(project_config_path,
                                            [vid_path],
                                            videotype='.mp4')

print("\nConvert detections to tracklets\n")
deeplabcut.convert_detections2tracklets(project_config_path,
                                        [vid_path],
                                        videotype='.mp4',
                                        shuffle=1,
                                        trainingsetindex=0,
                                        track_method='ellipse')

print("\nStitch tracklets together\n")
deeplabcut.stitch_tracklets(project_config_path,
                            [vid_path],
                            videotype='.mp4',
                            shuffle=1,
                            trainingsetindex=0,
                            track_method='ellipse')

print("\nFilter stitched tracklets (predictions)\n")
deeplabcut.filterpredictions(project_config_path,
                             [vid_path],
                             track_method='ellipse')

print("\nPlot trajectories\n")
deeplabcut.plot_trajectories(project_config_path,
                             [vid_path],
                             filtered=True,
                             track_method='ellipse')

# Complete flag file
np.save(os.path.join(os.environ.get('VIDEO_PATH'),'done.npy'),[])
