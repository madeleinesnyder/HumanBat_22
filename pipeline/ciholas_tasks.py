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
from shared_utils.utils import PyMatlab
import mat73
import glob

class ExtractCiholasData(luigi.Task):
    bat_id = luigi.Parameter()
    date = luigi.Parameter() # YYMMDD
    data_path = luigi.Parameter()

    def requires(self):
        return [PullServerData(self.bat_id, self.date, self.data_path)]

    def output(self):
        # Get logger directory path
        self.in_path = os.path.join(os.path.join(self.data_path, 'ciholas'))
        self.in_file = glob.glob(os.path.join(self.in_path, '*_cdp_*.txt'))[0]
        print(self.in_file)

        # Create output path
        self.out_path = os.path.dirname(self.in_path.replace('raw','processed'))
        Path(self.out_path).mkdir(parents=True, exist_ok=True)
        print(self.out_path)

        return luigi.LocalTarget(os.path.join(self.out_path,'doneciholas.npy'))

    def run(self):
        # Extract logger data
        HumanBatPath = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        pyMatlab = PyMatlab(HumanBatPath)
        print(HumanBatPath)

        print("Checking Ciholas integrity")
        expected_num_lines = 9000000
        #pyMatlab.eng.HumanBat_checkCiholasIntegrity(self.in_file, expected_num_lines, nargout = 0);
        print("Ciholas looks good!")

        print("Extracting ciholas data... ")
        print(self.in_path)
        tagSNs = [17106917,17106934,17107055,17107419,17107431,17107443,17106951]
        pyMatlab.eng.ExtractCdp_AF_v0(self.in_file, 'Include', tagSNs, nargout=0)

        # Copy extracted_data/ to processed folder
        dir_util.copy_tree(os.path.join(self.in_path),self.out_path)
