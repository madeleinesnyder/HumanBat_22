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
import glob
from distutils import dir_util
import getopt, sys
import luigi.tools.deps_tree as deps_tree
import process_ephys, process_video, process_cortex, process_ciholas
from utils import savemat73


class B149fCheckDataIntegrity(luigi.Task):
    data_path = luigi.Parameter()
    task_complete = False

    with open('./config.json', 'r') as f:
        config = json.load(f)

    def output(self):
        return luigi.LocalTarget(os.path.join(self.data_path,'good_integrity.npy'))


    def run(self):
        ephys_path = os.path.join(self.data_path, 'ephys')

        # EVENTLOG.CSV exists (extracted by EVENT_FILE_READER.EXE from NEURALYNX)
        dirs = [a for a in os.listdir(ephys_path) if os.path.isdir(os.path.join(ephys_path,a))]

        r = re.compile("(.*\d\d)")
        logger_dirs = list(filter(r.match, dirs))
        print("logger directories found: {}".format(logger_dirs))
        assert len(logger_dirs) == 1, "There should only be 1 logger folder! Found {}".format(logger_dirs)
        assert os.path.exists(os.path.join(ephys_path, logger_dirs[0], 'EVENTLOG.CSV')), "Data Integrity Error: EVENTLOG.CSV not found. Did you remember to extract it?"

        # TODO: Check number of NEUR____.DAT files match expected number in EVENTLOG.CSV

        # Checking ciholas data file is proper
        process_ciholas.check_ciholas_integrity(os.path.join(self.data_path,'ciholas'), 9000000)

        self.task_complete = True
        np.save(os.path.join(self.data_path,'good_integrity.npy'), np.array([1]))

class B149fExtractCiholasData(luigi.Task):
    data_path = luigi.Parameter()

    with open('./config.json', 'r') as f:
        config = json.load(f)

    def requires(self):
        return B149fCheckDataIntegrity(self.data_path)

    def output(self):
        self.out_path = os.path.join(os.path.join(self.data_path, 'ciholas')).replace('raw','processed')
        Path(self.out_path).mkdir(parents=True, exist_ok=True)
        self.in_path = os.path.join(os.path.join(self.data_path, 'ciholas'))
        self.in_file = glob.glob(os.path.join(self.in_path, '*_cdp_*.txt'))[0]
        print(self.in_file)
        return luigi.LocalTarget(os.path.join(self.out_path,'extracted_{}.mat'.format(Path(self.in_file).stem)))

    def run(self):
        tag_1_sn = self.config['b149f']['ciholas']['tag_1']['serial_number']
        tag_2_sn = self.config['b149f']['ciholas']['tag_2']['serial_number']
        tag_SNs = [tag_1_sn, tag_2_sn]
        process_ciholas.extract_ciholas(self.in_file, tag_SNs)
        dir_util.copy_tree(os.path.join(self.in_path),self.out_path)

class B149fExtractCortexData(luigi.Task):
    data_path = luigi.Parameter()

    with open('./config.json', 'r') as f:
        config = json.load(f)

    def requires(self):
        return B149fCheckDataIntegrity(self.data_path)

    def output(self):
        self.out_path = os.path.join(os.path.join(self.data_path, 'cortex')).replace('raw','processed')
        Path(self.out_path).mkdir(parents=True, exist_ok=True)
        self.in_path = os.path.join(os.path.join(self.data_path, 'cortex'))

        return luigi.LocalTarget(os.path.join(self.out_path,'done.npy'))

    def run(self):
        session_name = Path(self.data_path).name
        process_cortex.extract_cortex_c3d(self.in_path)
        dir_util.copy_tree(os.path.join(self.in_path,'Generated_C3D_files/processed'),self.out_path)
        np.save(os.path.join(self.out_path,'done.npy'), np.array([1]))

class B149fExtractEphysData(luigi.Task):
    data_path = luigi.Parameter()

    with open('./config.json', 'r') as f:
        config = json.load(f)

    def requires(self):
        return B149fCheckDataIntegrity(self.data_path)

    def output(self):
        # Get logger directory path
        self.in_path = os.path.join(os.path.join(self.data_path, 'ephys'))
        dirs = [a for a in os.listdir(self.in_path) if os.path.isdir(os.path.join(self.in_path,a))]

        r = re.compile("(.*\d\d)")
        logger_dirs = list(filter(r.match, dirs))
        self.in_path = os.path.join(self.in_path, logger_dirs[0])

        # Create output path
        self.out_path = os.path.dirname(self.in_path.replace('raw','processed'))
        Path(self.out_path).mkdir(parents=True, exist_ok=True)

        return luigi.LocalTarget(os.path.join(self.out_path,'extracted_data'))

    def run(self):
        fs = self.config['b149f']['ephys']['fs']

        # Extract logger data
        process_ephys.extract(self.in_path)

        # Copy extracted_data/ to processed folder
        if(os.path.isdir(os.path.join(self.out_path, 'extracted_data'))):
            print("Overwriting existing {}".format(os.path.join(self.out_path, 'extracted_data')))
        dir_util.copy_tree(os.path.join(self.in_path, 'extracted_data'),self.output().path)

class B149fEphysPowerSpectrum(luigi.Task):
    data_path = luigi.Parameter()

    with open('./config.json', 'r') as f:
        config = json.load(f)
        print(config)

    def requires(self):
        return B149fExtractEphysData(self.data_path), B149fCheckDataIntegrity(self.data_path)

    def output(self):
        p = os.path.join(os.path.join(self.data_path, 'ephys'))
        p = self.data_path.replace('raw','processed')

        self.in_path = os.path.join(p,'ephys/extracted_data')
        self.out_path = os.path.join(p,'tests/ephys_power_spectrum.jpg')

        Path(os.path.dirname(self.out_path)).mkdir(parents=True, exist_ok=True)

        return luigi.LocalTarget(self.out_path)

    def run(self):
        # Calculate power spectrum
        res = b151_tests.test_ephys_noise(self.in_path, self.out_path) # b151_tests.test_ephys_noise is identical to b149f
        assert res == True, 'Ephys test failed, see {} for test result'.format(self.out_path)


class B149fDownsampleEphysData(luigi.Task):
    data_path = luigi.Parameter()

    with open('./config.json', 'r') as f:
        config = json.load(f)

    def requires(self):
        return (B149fCheckDataIntegrity(self.data_path), B149fExtractEphysData(self.data_path))

    def output(self):
        # Get logger directory path
        self.in_path = os.path.join(os.path.join(self.data_path, 'ephys'))

        # Create output path
        self.out_path = self.in_path.replace('raw','processed')

        return luigi.LocalTarget(os.path.join(self.out_path,'logger_data.mat'))

    def run(self):
        ephys_ds = process_ephys.format_extracted_logger_data(os.path.join(self.out_path))

class B149fKilosortEphysData(luigi.Task):
    data_path = luigi.Parameter()
    resources = {"gpu": 1} # Uses GPU resources. Prevents multiple tasks from simultaneously using GPU

    with open('./config.json', 'r') as f:
        config = json.load(f)

    def requires(self):
        return (B149fCheckDataIntegrity(self.data_path), B149fDownsampleEphysData(self.data_path))

    def output(self):
        # Get logger directory path
        self.in_path = os.path.join(self.data_path.replace('raw','processed'),'ephys')

        # Create output path
        self.out_path = os.path.join(self.data_path.replace('raw','processed'),'ephys/params.py')

        return luigi.LocalTarget(self.out_path)

    def run(self):
        # Sort with kilosort2
        process_ephys.run_kilosort2(self.in_path)
        
        
class B149fConcatAudio(luigi.Task):
    data_path = luigi.Parameter()

    with open('./config.json', 'r') as f:
        config = json.load(f)

    #def requires(self):
    #    return (B149fCheckDataIntegrity(self.data_path), B149fExtractEphysData(self.data_path))

    def output(self):
        # Get logger directory path
        self.in_path = os.path.join(os.path.join(self.data_path, 'audio'))

        # Create output path
        self.out_path = self.in_path.replace('raw','processed')

        return luigi.LocalTarget(os.path.join(self.out_path,'logger_data.mat'))

    def run(self):
        ephys_ds = process_audio.audio_concat(os.path.join(self.out_path))

