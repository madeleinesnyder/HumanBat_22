import os
import numpy as np
import re
from pathlib import Path

def parse_b151_arduino_logs(log_file_path):
    arduino_ttl = []
    arduino_state = []
    arduino_state_enter_ms = []
    arduino_state_exit_ms = []
    with open(log_file_path, 'r') as f:
        for line in f.readlines():
            ttl_ts_matches= re.findall('(\d*)!', line) # Match ttl timestamps
            if(not ttl_ts_matches == []):
                arduino_ttl += ttl_ts_matches

            state_enter_matches=re.search('(\w+)_ENTER:(\d+)|', line)
            if(state_enter_matches.group(1) != None):
                arduino_state.append(state_enter_matches.group(1))
                arduino_state_enter_ms.append(state_enter_matches.group(2))

            state_exit_matches=re.search('.*(\w+)_EXIT:(\d+)|', line)
            if(state_exit_matches.group(1) != None):
                arduino_state_exit_ms.append(state_exit_matches.group(2))
    arduino_ttl = np.array(arduino_ttl, dtype=np.uint64)
    arduino_state_enter_ms = np.array(arduino_state_enter_ms, dtype=np.int64) - arduino_ttl[0]
    arduino_state_exit_ms = np.array(arduino_state_exit_ms, dtype=np.int64) - arduino_ttl[0]
    arduino_ttl_ms = arduino_ttl - arduino_ttl[0]

    # Arduino clock drift correction
    # Basic correction by interpolating ppms (parts per millisecond)
    # between first and last TTL
    known_ttl_interval = 1000 # 1 sec ttl
    num_ttl = len(arduino_ttl_ms)
    true_ttl_ts = (num_ttl-1) * known_ttl_interval # What the final ttl timestamp should be theoretically
    meas_ttl_ts = arduino_ttl_ms[-1]
    print("Arduino clock drift correction: true_ttl_ts: {}, meas_ttl_ts: {}, diff: {}".format(true_ttl_ts, meas_ttl_ts, meas_ttl_ts-true_ttl_ts))
    diff = meas_ttl_ts - true_ttl_ts
    ppms = diff/meas_ttl_ts # Clock drift per ms

    corrected_ttl_ms = (1-ppms)*arduino_ttl_ms
    corrected_state_enter_ms = (1-ppms)*arduino_state_enter_ms
    corrected_state_exit_ms = (1-ppms)*arduino_state_exit_ms


    return {'corrected_ttl_ms':corrected_ttl_ms,\
            'state':arduino_state,\
            'corrected_state_enter_ms':corrected_state_enter_ms,\
            'corrected_state_exit_ms':corrected_state_exit_ms}
