'''
collections of frequently used functions
'''
import numpy as np
import h5py,os


def read_raw_signal(fast5, start = None, end = None):
    '''
	dependency: h5py
		
	take the fast5 filename as the input, return a np array of raw signal
	input:
		start, end: <int>
			index (0-based) for truncating the output signal

    '''
    if start and end:
        start = int(start)
        end = int(end)
    with h5py.File(fast5, 'r') as h5_f: 
        
        read_key = list(h5_f["Raw/Reads/"].keys())[0]
        signal = list(h5_f["Raw/Reads/"][read_key]["Signal"])                                                                                                       
        if not start and not end:
            return(signal[::1])
        elif start < end:
            return(signal[start:end])
        else:
            print("Warning!! The signal is fetched backwards")
            return(signal[start:end:-1])
        try:
            read_key = list(h5_f["Raw/Reads/"].keys())[0]
            signal = list(h5_f["Raw/Reads/"][read_key]["Signal"])                                                                                                       
            if start < end:
                return(signal[start:end])
            else:
                return(signal[start:end:-1])
        except:
            exit("Error: Fail to read {}".format(fast5))


def normalization(signal, method = "z_score"):
	'''
	dependency: numpy as np
	input: signal dataset from fast5
	output: nomalised signal <np.array>
	'''
	signal = np.array(signal)
	
	if method == "z_score":
		return((signal - np.mean(signal)) / np.std(signal))
	
	elif method == "median_shift":    
	  	mad_scale = np.median(abs(signal - np.median(signal)))
	  	norm_signal = (signal - np.median(signal)) / mad_scale
	  	return(norm_signal)
	  
	else:
	  	exit("ERROR during normalization. The normalization method is not recognized. The accepted method: z_score,  meadian-shift")




def reverse_complement(seq):
    '''
	input: <str>
		queried seq
	return: <str>
		reverse_complement seq
	'''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    for element in seq:
        if element not in complement.keys():
            exit("Error in reverse_complement!!Bad bases in the seq.")
    letters = [complement[base] for base in seq]
    return ''.join(letters)[::-1]

def get_junction_signal_by_pos(fast5, junction_pos = None, window = 20, start_pos = None , end_pos = None):

    '''
	fast5: signal file with tombo output (a transcript reference needed)
	junction_pos: 0-based position of the splicing site of each transcript reference
	window: signals assigned to a window around the junction_pos will be fetched.
	'''
    if start_pos and start_pos:
        junction_bases_start, junction_bases_end = start_pos, end_pos
    elif junction_pos and window:
        junction_bases_start = junction_pos - int(np.floor(window/2))
        junction_bases_end = junction_pos + int(np.ceil(window/2))

    with h5py.File(fast5, 'r') as h5_f:
        
        path = "Analyses/"
       
        subpath = list(h5_f[path].keys())
        for i in subpath:
            if "RawGenomeCorrected" in i:
                path = path + i +'/BaseCalled_template/'


        try:
            h5_f[path]["Alignment"]
        except:
            print("Alignment doesn't exist in fast5!!")
            exit(0)  


        mapped_start = h5_f[path]["Alignment"].attrs['mapped_start']
        mapped_end = h5_f[path]["Alignment"].attrs['mapped_end']
        assert mapped_end - mapped_start > window,\
                                "window size larger than the read"

        strand = h5_f[path]["Alignment"].attrs['mapped_strand']
        Events = h5_f[path]["Events"]
        read_start_rel_to_raw = Events.attrs['read_start_rel_to_raw']
        
        #from transcript related pos to read specific pos

        
        if junction_bases_start >=0 and junction_bases_end >= 0:
            junction_bases_start -= mapped_start
            junction_bases_end -= mapped_start

        # transcript ref from the reverse strand of genome ref    
        elif junction_bases_start <=0 and junction_bases_end <= 0:
            junction_bases_start = mapped_end - abs(end_pos)  - mapped_start
            junction_bases_end = mapped_end - abs(start_pos) - mapped_start


        if junction_bases_start < 0:
            print("Warning: Read discarded!!junction pos is to close to the start of the mapped region")
            return 0
            junction_bases_start = 0
            junction_bases_end = window

        if junction_bases_end > mapped_end - mapped_start:
            print("Warning: Read discarded!!junction pos is to close to the end of the mapped region")
            return 0
            junction_bases_end = mapped_end - mapped_start
            junction_bases_start = mapped_end - mapped_start - window
        
        if strand == '-':
            Events = np.array(Events)[::-1]
            junction_signal_start = read_start_rel_to_raw +\
                                    Events[junction_bases_end][2]
        
            junction_signal_end = read_start_rel_to_raw + \
                Events[junction_bases_start+1][2] + Events[junction_bases_start][3]
        
        if strand == "+":
            junction_signal_start = read_start_rel_to_raw +\
                                    Events[junction_bases_start][2]
        
            junction_signal_end = read_start_rel_to_raw + \
                Events[junction_bases_end-1][2] + Events[junction_bases_start][3]
        
    print("Corresponding signal start and end:")
    print(junction_signal_start, junction_signal_end)
    return read_raw_signal(fast5, junction_signal_start, junction_signal_end)


def get_mapped_info_from_fast5(fast5,window = 20):
    """ extract:
            mapped_start
            mapped_end
            mapped_strand
            as dist from the tombo processed fast5
    """
    
    mapping_info = {}
    with h5py.File(fast5, 'r') as h5_f:
        path = "Analyses/"
        


        
        
        subpath = list(h5_f[path].keys())
        for i in subpath:
            if "RawGenomeCorrected" in i:
                path = path + i +'/BaseCalled_template/'
        
        try:
            h5_f[path]["Alignment"]
        except:
            print("Alignment doesn't exist in fast5!!")
            exit(0)
        
        mapped_start = h5_f[path]["Alignment"].attrs['mapped_start']
        mapped_end = h5_f[path]["Alignment"].attrs['mapped_end']
        assert mapped_end - mapped_start > window,\
                                "window size larger than the read"

        strand = h5_f[path]["Alignment"].attrs['mapped_strand']

        
        mapping_info["start"]=mapped_start
        mapping_info["end"]=mapped_end
        mapping_info["strand"]=strand
    
    return mapping_info

