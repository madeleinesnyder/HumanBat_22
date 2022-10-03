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

def check_ciholas_integrity(data_path, expected_num_lines):
    # Start matlab engine
    eng = matlab.engine.start_matlab()

    # Add q-bats path
    # Path is absolute
    QbatsPath = os.path.join(os.path.dirname(os.path.realpath(__file__)))
    QbatsPath = eng.genpath(QbatsPath)
    eng.addpath(QbatsPath, nargout=0)

    print("Checking Ciholas integrity")

    eng.HumanBat_checkCiholasIntegrity(data_path, expected_num_lines, nargout = 0);

    print("Ciholas looks good!")


def extract_ciholas(data_path, tag_SNs):
    """
    Python wrapper for extracting C3D data files from cortex

    Parameters
    ----------
    data_path : string
        Path to ciholas folder containing *_cdp_*.txt file
    tag_SNs:
        array of serial numbers of ciholas tags used
    """
    assert os.path.exists(data_path), "{} does not exist!".format(data_path)

    folder_name = pathlib.PurePath(data_path).name

    # Start matlab engine
    eng = matlab.engine.start_matlab()

    # Add q-bats path
    # Path is absolute
    QbatsPath = os.path.join(os.path.dirname(os.path.realpath(__file__)))
    QbatsPath = eng.genpath(QbatsPath)
    eng.addpath(QbatsPath, nargout=0)

    print("Running Angelo's Ciholas extraction script... ")
    print(data_path)

    eng.ExtractCdp_AF_v0(data_path,'Include', tag_SNs, nargout=0)

    print("Ciholas extraction complete!")
