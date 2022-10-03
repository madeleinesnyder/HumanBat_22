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
from pipeline.rclone_tasks import *
from pipeline.ephys_tasks import *
from shared_utils.utils import *

class DownsampleEphysData(luigi.Task):
    bat_id = luigi.Parameter()
    date = luigi.Parameter() # YYMMDD
    data_path = luigi.Parameter()
    resources = {"gpu":1} # Uses GPU resorces, prevents multiple tasks from simultaneously using GPU


    def requires(self):
        return [PullServerData(self.bat_id, self.date, self.data_path),ExtractEphysData(self.bat_id,self.date,self.data_path)]

    def output(self):
        # Get logger directory path
        print(self.data_path)
        self.in_path = os.path.join(os.path.join(self.data_path, 'ephys'))

        # Get all directories (not files) in path
        dirs = [a for a in os.listdir(self.in_path) if os.path.isdir(os.path.join(self.in_path,a))]
        r = re.compile("(.*\d\d)")
        self.logger_dirs = list(filter(r.match, dirs))

        # Create output path
        self.out_path = os.path.dirname(self.in_path.replace('raw','processed'))
        Path(self.out_path).mkdir(parents=True, exist_ok=True)

        return luigi.LocalTarget(os.path.join(self.out_path,'logger_data.mat'))

    def run(self):

        HumanBatPath = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        pyMatlab = PyMatlab(HumanBatPath)
        print(HumanBatPath)

        # Downsample ephys data to kilosort
        for logger in self.logger_dirs: # Extract all loggers
            print("Downsampling ephys {} data... ".format(logger))
            print(os.path.join(self.out_path, logger))
            pyMatlab.eng.format_extracted_logger_data(os.path.join(self.out_path, logger),nargout=0)
