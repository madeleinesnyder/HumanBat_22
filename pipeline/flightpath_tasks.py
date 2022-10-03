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

class ExtractFlightPaths(luigi.Task):
    bat_id = luigi.Parameter()
    date = luigi.Parameter() # YYMMDD
    data_path = luigi.Parameter()

    def requires(self):
        return [PullServerData(self.bat_id, self.date, self.data_path),
                ExtractCortexData(self.bat_id, self.date, self.data_path)]

    def output(self):
        # Get trackfile directory path
        self.in_path = os.path.join(os.path.join(self.data_path, 'cortex'))
        self.in_file = os.path.dirname(self.in_path.replace('raw','processed'))

        # Create output path
        self.out_path = os.path.dirname(self.in_path.replace('raw','processed'))
        Path(self.out_path).mkdir(parents=True, exist_ok=True)

        return luigi.LocalTarget(os.path.join(self.out_path,'done.npy'))

    def run(self):
        # Extract logger data
        HumanBatPath = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        pyMatlab = PyMatlab(HumanBatPath)
        print(HumanBatPath)

        print("Characterizing FlightPaths... ")
        print(self.in_path)
        #pyMatlab.eng.HumanBat_characterizeFlights(self.in_path nargout=0)

        print("Extracting and Clustering FlightPaths")
        print(self.in_path)
        #pyMatlab.eng.HumanBat_plot_trajectories(self.in_path, nargout=0)
