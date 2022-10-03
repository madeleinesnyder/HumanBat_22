import matlab.engine
import io
import sys
import os
import pathlib
import mat73
import re
import scipy.signal
import scipy.io
import numpy as np
from HumanBat.qbats.utils import PyMatlab

def concat_audio(data_path):
    """
    Python wrapper for audioconcat matlab script.

    Parameters
    ----------
    data_path : string
    
    TODO
    """
    assert os.path.exists(data_path), "{} does not exist!".format(data_path)

    folder_name = pathlib.PurePath(data_path).name
    # Start matlab engine
    eng = matlab.engine.start_matlab()

    # Add Path to audio matlab scripts
    # Path is absolute
    AudioDataProcessingPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'audio')
    AudioDataProcessingPath = eng.genpath(AudioDataProcessingPath)

    matlab_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), data_path)

    eng.addpath(AudioDataProcessingPath, nargout=0)

    print("Running audio concatination script... ")
    eng.HumanBat_audio_concat(matlab_path, nargout=0)

    print("Audio Concatination complete!")


def downsample_audio(data_path):
    """
    Python wrapper for downsampling the concatinated audio data
    
    Parameters
    ----------
    data_path : string 
    
    """
    assert os.path.exists(data_path), "{} does not exist!".format(data_path)

    folder_name = pathlib.PurePath(data_path).name
    # Start matlab engine
    eng = matlab.engine.start_matlab()

    # Add Path to audio matlab scripts
    # Path is absolute
    AudioDataProcessingPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'audio')
    AudioDataProcessingPath = eng.genpath(AudioDataProcessingPath)

    matlab_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), data_path)

    eng.addpath(AudioDataProcessingPath, nargout=0)

    print("Running audio downsampling script... ")
    eng.HumanBat_downsample_concatinated_audio(matlab_path, nargout=0)

    print("Audio downsampling complete!")
    
"""
def load_extracted_data(extracted_ephys_path, fs, ds_factor):
    #print(os.listdir(extracted_ephys_path))
    fnames = os.listdir(extracted_ephys_path)
    r = re.compile('(.*CSC)\d*.mat')
    fname = list(filter(r.match, fnames))[0]
    fname = re.search('(.*)CSC\d*.mat', fname).group(1)

    event_data = scipy.io.loadmat(os.path.join(extracted_ephys_path, fname+'EVENTS.mat'))
    event_names = event_data['event_types_and_details']
    event_ts = event_data['event_timestamps_usec']
    ttl_ts = []
    for i in range(len(event_names)):
        name = ''.join(event_names[i][0])
        print(name)
        if('Started recording' in name):
            t0 = event_ts[i]
        if('Digital in rising edge' in name):
            ttl_ts.append(event_ts[i])
    print(t0)
    ttl_ts = np.array(ttl_ts) - t0
    print(ttl_ts)
    ttl_t0 = int(ttl_ts[0]*(fs/ds_factor)/(1e6))
    print(ttl_t0)

    ephys_data = []
    for i in range(16):
        print("Loading ephys channel {}...".format(i))
        ch_data = mat73.loadmat(os.path.join(extracted_ephys_path, fname+'CSC{}.mat'.format(i)))['AD_count_int16']

        # Resample
        if(ds_factor>1):
            ephys_ds = scipy.signal.resample_poly(ch_data, 1, ds_factor)
        else:
            ephys_ds = ch_data
        print(ephys_ds.shape)
        ephys_data.append(ephys_ds)

    return {'data': np.vstack(ephys_data)[:,ttl_t0:], 'ttl_ts': ttl_ts, 'fs':fs/ds_factor}"""

def format_extracted_logger_data(extracted_ephys_path):
    pylab = PyMatlab(os.path.dirname(os.path.realpath(__file__)))
    print("Formatting extracted logger data... ")
    pylab.eng.format_extracted_logger_data(extracted_ephys_path, nargout=0)

def extracted2binary(extracted_ephys_path):

    # Calls matlab script to extract ephys to a binary file

    # Start matlab engine
    eng = matlab.engine.start_matlab()

    # Add LoggerDataProcessing Path
    # Path is absolute
    LoggerDataProcessingPath = os.path.dirname(os.path.realpath(__file__))
    LoggerDataProcessingPath = eng.genpath(LoggerDataProcessingPath)
    eng.addpath(LoggerDataProcessingPath, nargout=0)

    print("Running extraction2binary script... ")
    eng.format_extracted_logger_data(extracted_ephys_path, nargout=0)

    return

def run_kilosort2(extracted_ephys_path):

    # Runs kilosort2

    # Start matlab engine
    eng = matlab.engine.start_matlab()

    # Add LoggerDataProcessing Path
    # Path is absolute
    LoggerDataProcessingPath = os.path.dirname(os.path.realpath(__file__))
    LoggerDataProcessingPath = eng.genpath(LoggerDataProcessingPath)
    eng.addpath(LoggerDataProcessingPath, nargout=0)

    print("Running kilosort2 script... ")
    eng.mymaster_kilosort('',extracted_ephys_path,nargout=0)

    return
