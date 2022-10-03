import numpy as np
import os
import scipy.signal

def load_motu_data(motu_path):
    motu_data = []
    num_files = len(os.listdir(motu_path))
    print("loading audio data...")
    for i in range(num_files):
        try: # Try block handles the occasional error where the last data block is not saved properly.
            data = np.load(os.path.join(motu_path,'audio_{}.npy'.format(i)))
            motu_data.append(data)
        except ValueError as e:
            print(e)
        #data_ds = scipy.signal.decimate(data, 5)

    motu_data = np.hstack(motu_data)
    return motu_data

def get_motu_ttl_indices(motu_data, motu_fs, channel_map = {'audio':0,"ttl":1,"events":2}):
    """
    Get TTL indices from MOTU ttl channel

    Parameters
    ----------
    motu_data : ndarray (num_channels, num_samples)

    motu_fs : int
        MOTU sampling frequency

    channel_map : dict
        Mapping of each MOTU channel.
    """
    ttl_channel = motu_data[channel_map['ttl'],:]
    peak_x = np.array(scipy.signal.find_peaks(ttl_channel,distance=motu_fs/1.5, height=np.max(ttl_channel)/2)[0],dtype=np.uint64)
    return peak_x

def slice_valid_motu_data(motu_data, motu_fs):
    """
    Slice MOTU data based on first and last TTLs

    Parameters
    ----------
    motu_data : ndarray (num_channels, num_samples)

    motu_fs : int
        MOTU sampling frequency

    Returns
    -------
    (
        ndarray(num_channels, num_samples)
            Sliced MOTU data
        ,
        ndarray(num_ttls)
            TTL indices
    )
    """
    motu_ttl_indices = get_motu_ttl_indices(motu_data, motu_fs)

    motu_first_ttl_index = motu_ttl_indices[0]
    motu_last_ttl_index = motu_ttl_indices[-1]

    valid_motu_data = motu_data[:,motu_first_ttl_index:motu_last_ttl_index]
    shifted_ttl_indices = get_motu_ttl_indices(valid_motu_data, motu_fs)

    return valid_motu_data, shifted_ttl_indices
