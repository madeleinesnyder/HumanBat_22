import matlab.engine
import io
import sys
import os
import pathlib
import numpy as np
import matplotlib.pyplot as plt
import cv2
import pickle
import gc
import scipy.signal

def test_ephys_noise(ephys_extracted_path, out_path):
    """
    ephys_extracted_path : string
        Path to extracted_data/ folder from Julie's extract_logger_data script

    out_path : string
        Test result output directory
    """
    # Start matlab engine
    eng = matlab.engine.start_matlab()

    # Add LoggerDataProcessing Path
    # Path is absolute
    matlab_script_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),'tests')
    matlab_script_path = eng.genpath(matlab_script_path)
    q_bats_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),'qbats')
    q_bats_path = eng.genpath(q_bats_path)
    eng.addpath(matlab_script_path, nargout=0)
    eng.addpath(q_bats_path, nargout=0)

    print("Calculating Ephys Power Spectrum.... ")
    res = eng.test_ephys_noise(ephys_extracted_path,out_path,nargout=1)

    return res

def visualize_synchrony(motu_path, arduino_path, ephys_path, camera_path, out_path, config):
    motu_data = np.load(motu_path, allow_pickle=True)['motu']
    motu_fs = motu_data['fs']
    event_motu = motu_data['data'][2,:]
    motu_data = None
    gc.collect()

    arduino_data = np.load(arduino_path, allow_pickle=True)
    arduino_data = arduino_data.tolist()['arduino']

    # TODO: Ephys Data
    ephys = np.load(ephys_path, allow_pickle=True).tolist()
    ephys_data = ephys['data']
    ephys_fs = ephys['fs']
    ephys = None
    gc.collect()

    state_mapping = {'RESET':0, 'READY':4, 'STIM':3, 'OPEN':2, 'REW':1}
    state_history = np.array(arduino_data['state'])

    camera_paths = {}
    for dir in os.listdir(camera_path)[:]:
        for f in os.listdir(os.path.join(camera_path, dir))[:]:
            fpath = os.path.join(camera_path,dir,f)
            camera_paths[dir] = fpath

    fps = 50
    T = 1/fps
    frame_step_size = 10
    shift = 2 #sec
    delta_t = frame_step_size*T

    motu_step_size = T * motu_fs
    ms_step_size = T * 1000
    ephy_step_size = T * motu_fs

    last_ready_enter_index = np.where(state_history == 'READY')[0][-2]
    t0_ms = int(arduino_data['corrected_state_enter_ms'][last_ready_enter_index] - shift*1000) # ms

    duration = 75000 # ms
    states = arduino_data['state'][last_ready_enter_index:]
    state_duration = np.diff(np.hstack([arduino_data['corrected_state_enter_ms'][last_ready_enter_index:]]))
    state_waveform = []
    for i in range(len(states)-1):
        state_waveform += (int(state_duration[i])*int(motu_fs/1000))*[state_mapping[states[i]]]
    state_waveform = np.array(state_waveform)[:duration*int(motu_fs/1000)]
    state_waveform = np.hstack([[0]*shift*motu_fs,state_waveform])
    duration = np.min([len(state_waveform), duration*int(motu_fs/1000)])

    t0_motu = int(t0_ms*motu_fs/1000)
    t0_camera = int(t0_ms*(fps/1000))
    t0_ephy = int(t0_ms*(ephys_fs/1000))
    print('Synchrony check segment:')
    print('t0_camera', t0_camera)
    print('t0_motu',t0_motu)
    print('t0_ms',t0_ms)
    print('t0_ephy',t0_ephy)
    print('duration (samples at motu_fs)', duration)
    #print(t0_ephy+int(duration*(ephys_fs/motu_fs)))
    print(ephys_data.shape)
    print(int(duration*(ephys_fs/motu_fs)))
    print(t0_ephy)
    print(ephys_fs)
    ephys_data = np.mean(ephys_data[:,t0_ephy:t0_ephy+int(duration*(ephys_fs/motu_fs))],axis=0)
    print(ephys_data.shape)
    print(duration)
    ephys_data = scipy.signal.resample(ephys_data, duration)

    fourcc = cv2.VideoWriter_fourcc('m', 'p', '4', 'v')
    writer = cv2.VideoWriter(out_path,fourcc,fps = int(fps/frame_step_size), frameSize=(1200, 800))

    caps = {}
    for camera in camera_paths.keys():
        caps[camera] = cv2.VideoCapture(camera_paths[camera])
        print(int(caps[camera].get(cv2.CAP_PROP_FRAME_COUNT)))
        caps[camera].set(cv2.CAP_PROP_POS_FRAMES, t0_camera)

    fig = plt.figure(figsize=(12,8), dpi=100)
    ax1 = plt.subplot(4,3,3)
    ax2 = plt.subplot(4,3,6,sharex=ax1)
    ax3 = plt.subplot(4,3,9,sharex=ax1)
    ax4 = plt.subplot(4,3,12)
    ax5 = plt.subplot(2,3,1)
    ax6 = plt.subplot(2,3,2)
    ax7 = plt.subplot(2,3,4)
    ax8 = plt.subplot(2,3,5)
    vid_axs = [ax5,ax6,ax7,ax8]
    im_axs = {}
    ax1.plot(np.abs(event_motu[t0_motu:t0_motu+duration]))
    ax2.plot(state_waveform)
    ax1_vline = ax1.axvline(0, 0,1, color='red')
    ax2_vline = ax2.axvline(0, 0,5, color='red')
    ax2.set_ylim([0,5])

    ax4_text = ax4.text(0.5, 0.5, 'Camera Frame: {}'.format(0), horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, fontsize='large')
    ax4.axis('off')

    ax3.plot(ephys_data)
    ax3_vline = ax3.axvline(0,0,2024,color='red')
    ax3.axis('off')

    for i in range(int(duration*(1/motu_fs)*50)):
        frames = {}
        for camera in caps.keys():
            #print(camera)
            ret, frames[camera] = caps[camera].read()
            #print(frames[camera].shape)

        if(i%frame_step_size == 0):
            if i%(frame_step_size*10) == 0:
                print(i)

            j = 0
            for key in frames.keys():
                ax = vid_axs[j]
                title = key
                if(i==0):
                    im_axs[key] = ax.imshow(frames[key])
                    ax.axis('off')
                else:
                    im_axs[key].set_data(frames[key])
                j += 1

            ax1_vline.set_xdata([i*motu_step_size])
            ax2_vline.set_xdata([i*motu_step_size])
            ax3_vline.set_xdata([i*motu_step_size])

            ax4_text.set_text('Camera Frame: {}'.format(i+t0_camera))

            fig.canvas.draw()
            # convert canvas to image
            img = np.frombuffer(fig.canvas.tostring_rgb(), dtype = np.uint8)
            img = img.reshape(fig.canvas.get_width_height()[::-1] + (3, ))
            img = cv2.cvtColor(img, cv2.COLOR_RGB2BGR)
            writer.write(img)
    for camera in camera_paths.keys():
        caps[camera].release()
    cv2.destroyAllWindows()
    writer.release()

    return True
