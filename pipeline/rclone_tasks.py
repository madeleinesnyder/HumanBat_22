import luigi
import os
from pathlib import Path
#import h5py
import json
import sys
import numpy as np
#import hdf5storage
import re
import gc
import pickle
from distutils import dir_util
import getopt, sys
import luigi.tools.deps_tree as deps_tree
#from pipeline.audio_tasks import *
#from pipeline.video_tasks import *
from shared_utils import rclone
#from HumanBat import process_motu, process_arduino, process_ephys, process_video, b151_tests, DockerTask
#from process_b149f import *
#from HumanBat.qbats.utils import savemat73

class PullServerData(luigi.Task):
    """
    Pulls data from NAS servers for specified date

    Pulls both raw/ and processed/ for specified date
    """
    bat_id = luigi.Parameter()
    date = luigi.Parameter() # YYMMDD
    data_path = luigi.Parameter()

    #server_path = luigi.Parameter() # Path to server root directory containing data from each day
    #local_path = luigi.Parameter() # Path to local root directory containing data from each day

    def output(self):
        script_dir = os.path.dirname(__file__)
        with open(os.path.join(script_dir,'config.json'),'r') as f:
            config = json.load(f)

        server_path = config['remote_data_path']
        self.server_path = os.path.join(server_path,'{}/raw/{}'.format(self.bat_id, self.date))

        local_path = config['local_data_path']
        self.local_path = os.path.join(local_path,'{}/raw/{}'.format(self.bat_id, self.date))

        rclone.check(self.server_path, self.local_path)

        return luigi.LocalTarget(self.local_path)

    def run(self):
        print("\n\nPulling raw data for bat {} from {}\n".format(self.bat_id, self.date))
        rclone.copy(self.server_path, self.local_path)
        print("\n\nPulling processed data for bat {} from {}\n".format(self.bat_id, self.date))
        rclone.copy(self.server_path.replace('raw','processed'), self.local_path.replace('raw','processed'))

class PushServerData(luigi.Task):
    """
    Pushes data to NAS servers for specified date(s)

    Pushes processed/ folder only (rclone copy command, never deletes any files so it is safe)
    """
    bat_id = luigi.Parameter()
    date = luigi.Parameter() # YYMMDD
    data_path = luigi.Parameter()

    #server_path = luigi.Parameter() # Path to server root directory containing data from each day
    #local_path = luigi.Parameter() # Path to local root directory containing data from each day

    def requires(self):
        return [PullServerData(self.bat_id, self.date, self.data_path),
                ConcatAudio(self.bat_id, self.date, self.data_path),
                StackVideos(self.bat_id, self.date, self.data_path)]

    def output(self):
        script_dir = os.path.dirname(__file__)
        with open(os.path.join(script_dir,'config.json'),'r') as f:
            config = json.load(f)

        server_path = config['remote_data_path']
        self.server_path = os.path.join(server_path,'{}/processed/{}'.format(self.bat_id, self.date))

        local_path = config['local_data_path']
        self.local_path = os.path.join(local_path,'{}/processed/{}'.format(self.bat_id, self.date))

        rclone.check(self.local_path,self.server_path)

        return None

    def run(self):
        rclone.copy(self.local_path, self.server_path)

"""
class B151CheckDataIntegrity(luigi.Task):
    data_path = luigi.Parameter()
    task_complete = False

    with open('./config/config.json', 'r') as f:
        config = json.load(f)

    def output(self):
        return None

    def complete(self):
        return True#self.task_complete

    def run(self):
        ephys_path = os.path.join(self.data_path, 'b151/ephys')
        logger_dirs = os.listdir(ephys_path)

        assert len(logger_dirs) == 1, "Data Integrity Error: Multiple logger folders found in {}".ephys_path

        # EVENTLOG.CSV exists (extracted by EVENT_FILE_READER.EXE from NEURALYNX)
        assert os.path.exists(os.path.join(ephys_path, logger_dirs[0], 'EVENTLOG.CSV')), "Data Integrity Error: EVENTLOG.CSV not found. Did you remember to extract it?"

        # TODO: Check number of NEUR____.DAT files match expected number in EVENTLOG.CSV

        # Check arduino logs exist
        session_name = Path(self.data_path).name
        assert os.path.exists(os.path.join(self.data_path, 'b151/{}_logs.txt'.format(session_name)))

        self.task_complete = True

class B151ExtractMotuData(luigi.Task):
    data_path = luigi.Parameter()
    print(data_path)

    with open('./config/config.json', 'r') as f:
        config = json.load(f)

    def requires(self):
        return B151CheckDataIntegrity(self.data_path)

    def output(self):
        self.out_path = os.path.join(os.path.join(self.data_path, 'b151/motu')).replace('raw','processed')
        Path(self.out_path).mkdir(parents=True, exist_ok=True)
        self.in_path = os.path.join(os.path.join(self.data_path, 'b151/motu'))

        return luigi.LocalTarget(os.path.join(self.out_path,'motu.pkl'))

    def run(self):
        raw_motu_data = process_motu.load_motu_data(self.in_path)
        print(raw_motu_data.shape)
        fs = self.config['b151']['motu']['fs']

        motu_data, ttl_indices = process_motu.slice_valid_motu_data(raw_motu_data, fs)
        mat_data = {'motu': {'data': motu_data, 'ttl_indices': ttl_indices, 'fs': fs}}
        #hdf5storage.savemat(self.output().path, mat_data, format='7.3')
        with open(self.output().path, 'wb') as f:
            pickle.dump(mat_data, f, protocol=pickle.HIGHEST_PROTOCOL)

        savemat73(mat_data, self.output().path.replace('.pkl','.mat'))
        #np.save(os.path.join(self.out_path, 'motu_data.npy'), motu_data, allow_pickle=False)

class B151ExtractArduinoData(luigi.Task):
    data_path = luigi.Parameter()

    with open('./config/config.json', 'r') as f:
        config = json.load(f)

    def requires(self):
        return B151CheckDataIntegrity(self.data_path)

    def output(self):
        self.out_path = os.path.join(os.path.join(self.data_path, 'b151/arduino')).replace('raw','processed')
        Path(self.out_path).mkdir(parents=True, exist_ok=True)
        self.in_path = os.path.join(os.path.join(self.data_path, 'b151'))
        return luigi.LocalTarget(os.path.join(self.out_path,'arduino.npy'))

    def run(self):
        session_name = Path(self.data_path).name
        arduino_data = process_arduino.parse_b151_arduino_logs(os.path.join(self.in_path, '{}_logs.txt'.format(session_name)))

        arduino_data = {'arduino': arduino_data}
        #hdf5storage.savemat(self.output().path, arduino_data, format='7.3')
        np.save(self.output().path, arduino_data)

class B151ExtractEphysData(luigi.Task):
    data_path = luigi.Parameter()

    with open('./config/config.json', 'r') as f:
        config = json.load(f)

    def requires(self):
        return B151CheckDataIntegrity(self.data_path)

    def output(self):
        # Get logger directory path
        self.in_path = os.path.join(os.path.join(self.data_path, 'b151/ephys'))
        dirs = [a for a in os.listdir(self.in_path) if os.path.isdir(os.path.join(self.in_path,a))]

        r = re.compile("(.*\d\d)")
        logger_dirs = list(filter(r.match, dirs))
        print("logger directories found: {}".format(logger_dirs))
        assert len(logger_dirs) == 1, "There should only be 1 logger folder! Found {}".format(logger_dirs)
        self.in_path = os.path.join(self.in_path, logger_dirs[0])

        # Create output path
        self.out_path = os.path.dirname(self.in_path.replace('raw','processed'))
        Path(self.out_path).mkdir(parents=True, exist_ok=True)

        return luigi.LocalTarget(os.path.join(self.out_path,'extracted_data'))

    def run(self):
        fs = self.config['b151']['ephys']['fs']

        # Extract logger data
        process_ephys.extract(self.in_path)

        # Copy extracted_data/ to processed folder
        if(os.path.isdir(os.path.join(self.out_path, 'extracted_data'))):
            print("Overwriting existing {}".format(os.path.join(self.out_path, 'extracted_data')))
        dir_util.copy_tree(os.path.join(self.in_path, 'extracted_data'),self.output().path)

class B151DownsampleEphysData(luigi.Task):
    data_path = luigi.Parameter()

    def requires(self):
        return (B151CheckDataIntegrity(self.data_path), B151ExtractEphysData(self.data_path))

    def output(self):
        # Get logger directory path
        self.in_path = os.path.join(os.path.join(self.data_path, 'b151/ephys'))

        # Create output path
        self.out_path = self.in_path.replace('raw','processed')

        return luigi.LocalTarget(os.path.join(self.out_path,'logger_data.mat'))

    def run(self):
        ephys_ds = process_ephys.format_extracted_logger_data(os.path.join(self.out_path))


class B151KilosortEphysData(luigi.Task):
    data_path = luigi.Parameter()
    resources = {"gpu": 1} # Uses GPU resources. Prevents multiple tasks from simultaneously using GPU

    def requires(self):
        return (B151CheckDataIntegrity(self.data_path), B151ExtractEphysData(self.data_path), B151DownsampleEphysData(self.data_path))

    def output(self):
        # Get logger directory path
        self.in_path = os.path.join(self.data_path.replace('raw','processed'),'b151/ephys')

        # Create output path
        self.out_path = os.path.join(self.data_path.replace('raw','processed'),'b151/ephys/params.py')

        return luigi.LocalTarget(self.out_path)

    def run(self):
        # Run kilosort2
        process_ephys.run_kilosort2(self.in_path)

class B151ExtractCameraData(luigi.Task):
    data_path = luigi.Parameter()

    with open('./config/config.json', 'r') as f:
        config = json.load(f)

    def requires(self):
        return (B151CheckDataIntegrity(self.data_path))

    def output(self):
        # Get camera data path
        self.in_path = os.path.join(os.path.join(self.data_path, 'b151/cameras'))

        # Create output path
        self.out_path = self.in_path.replace('raw','processed')

        return luigi.LocalTarget(self.out_path)

    def run(self):
        # Extract camera data
        room_name = 'b151'
        session_name = Path(self.data_path).name
        process_video.preprocess_raw_video(room_name, session_name, self.data_path, self.config)

class B151BottomCameraDLC(DockerTask.DockerTask):
    data_path = luigi.Parameter()
    resources = {"gpu": 1} # Uses GPU resources. Prevents multiple tasks from simultaneously using GPU

    with open('./config/config.json', 'r') as f:
        config = json.load(f)

    @property
    def environment(self):
        return {'NVIDIA_VISIBLE_DEVICES': 'all', 'VIDEO_PATH':os.path.join(self.data_path.replace('raw','processed'), 'b151/cameras/bottom')}

    @property
    def image(self):
        return 'deeplabcut/deeplabcut:latest-core'
#ipython /app/dlc_analyze_videos.py
    @property
    def command(self):
        return 'ipython HumanBat/dlc_analyze_videos.py'

    @property
    def container_options(self):
        return {"working_dir": "/app"}

    @property
    def use_gpu(self):
        return True

    @property
    def auto_remove(self):
        return True

    def output(self):
        return luigi.LocalTarget(os.path.join(self.data_path.replace('raw', 'processed'),'b151/cameras/done.npy'))

    def requires(self):
        return B151ExtractCameraData(self.data_path)

    @property
    def binds(self):
        '''
        Override this to mount local volumes, in addition to the /tmp/luigi
        which gets defined by default. This should return a list of strings.
        e.g. ['/hostpath1:/containerpath1', '/hostpath2:/containerpath2']
        '''
        return ['/home/batlab/Desktop/HumanBat:/app:z']


class B151EphysPowerSpectrum(luigi.Task):
    data_path = luigi.Parameter()

    with open('./config/config.json', 'r') as f:
        config = json.load(f)
        print(config)

    def requires(self):
        return B151ExtractEphysData(self.data_path), B151CheckDataIntegrity(self.data_path)

    def output(self):
        p = os.path.join(os.path.join(self.data_path, 'b151/ephys/'))
        p = self.data_path.replace('raw','processed')

        self.in_path = os.path.join(p,'b151/ephys/extracted_data')
        self.out_path = os.path.join(p,'b151/tests/ephys_power_spectrum.jpg')

        Path(os.path.dirname(self.out_path)).mkdir(parents=True, exist_ok=True)

        return luigi.LocalTarget(self.out_path)

    def run(self):
        # Extract camera data
        res = b151_tests.test_ephys_noise(self.in_path, self.out_path)
        assert res == True, 'Ephys test failed, see {} for test result'.format(self.out_path)

class B151VisualizeSynchronyTest(luigi.Task):
    data_path = luigi.Parameter()

    with open('./config/config.json', 'r') as f:
        config = json.load(f)

    def requires(self):
        return (B151ExtractEphysData(self.data_path),\
                B151ExtractCameraData(self.data_path),\
                B151ExtractArduinoData(self.data_path),\
                B151ExtractMotuData(self.data_path),\
                B151DownsampleEphysData(self.data_path),\
                B151CheckDataIntegrity(self.data_path))

    def output(self):
        p = os.path.join(os.path.join(self.data_path, 'b151')).replace('raw','processed')

        self.ephys_path = os.path.join(p,'ephys/ephys_ds.npy')
        self.motu_path = os.path.join(p,'motu/motu.pkl')
        self.arduino_path = os.path.join(p,'arduino/arduino.npy')
        self.camera_path = os.path.join(p,'cameras/')
        self.out_path = os.path.join(p,'tests/visualize_synchrony.mp4')

        Path(os.path.dirname(self.out_path)).mkdir(parents=True, exist_ok=True)

        return luigi.LocalTarget(self.out_path)

    def run(self):
        # Visualize Synchrony
        res = b151_tests.visualize_synchrony(self.motu_path,self.arduino_path,self.ephys_path,self.camera_path,self.out_path,self.config)
        assert res == True, 'Visualize Synchrony test failed, see {} for test result'.format(self.out_path)
"""
