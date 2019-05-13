'''
collections of frequently used functions
'''
import numpy as np
import h5py



def read_raw_segnal(fast5, start = None, end = None):
    '''
    dependency: 
        h5py
        
    take the fast5 filename as the input, return a np array of raw signal
    input:
        start, end: <int>
            index (0-based) for truncating the output signal
    '''
    with h5py.File(fast5, 'r') as h5_f:
        
        try:
            read_key = list(h5_f["Raw/Reads/"].keys())[0]
            signal = h5_f["Raw/Reads/"][read_key]["Signal"]
            
            return(signal[int(start):int(end)])
        except:
            exit("Error: Fail to read {}".format(fast5))
            return(None)

def mean_normalization(signal):
    '''
    dependency: numpy as np
    input: signal dataset from fast5
    output: nomalised signal <np.array>
    '''
    signal = np.array(signal)
    mad_scale = np.median(abs(signal - np.median(signal)))

    norm_signal = (signal - np.median(signal)) / mad_scale


    return(norm_signal)