'''
collections of frequently used functions
'''

import numpy as np
import h5py
import os
from collections import defaultdict
import scrappy
import matplotlib.pyplot as plt
import sys


def read_raw_signal(fast5, start = None, end = None):
	'''
	take the fast5 filename as the input, return a np array of raw signal
	Args:
		start, end: <int>
			index (0-based) for truncating the output signal
	Returns: 
		raw signal <numpy array>
	'''
	start = int(start) if start else 0
	end = int(end) if end else None

	if end and end <= start:
		print("InputError: Invalid start and end position when fetching \
			when fetching the raw signals")
		sys.exit(0)
	
	with h5py.File(fast5, 'r') as h5_f: 
		read_key = list(h5_f["Raw/Reads/"].keys())[0]
		signal = list(h5_f["Raw/Reads/"][read_key]["Signal"])        
		return(signal[start:end])

#def normalization(signal, method = "z_score"):
#	'''
#	Args:
#		signal dataset from fast5
#	Returns:
#		 nomalised signal <np.array>
#	'''
#	signal = np.array(signal)
#	
#	if method == "z_score":
#		return((signal - np.mean(signal)) / np.std(signal))
#
#	elif method == "median_shift":    
#	  	mad_scale = np.median(abs(signal - np.median(signal)))
#	  	norm_signal = (signal - np.median(signal)) / mad_scale
#	  	return(norm_signal)
#
#	else:
#		print("ERROR during normalization. The normalization method is  \
#			not recognized. The accepted method: z_score,  meadian-shift")
#		sys.exit(0)
#



	
def reverse_complement(seq):
	'''
	Args: <str>
		queried seq
	Returns: <str>
		reverse_complement seq
	'''
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
	for element in seq:
		if element not in complement.keys():
			print("Error in reverse_complement!!Bad bases in the seq.")
			sys.exit(0)
	letters = [complement[base] for base in seq]
	return ''.join(letters)[::-1]

# TOMBO related functions
def get_junction_signal_by_pos(fast5, junction_pos = None, 
		window = None, start_pos = None , end_pos = None, verbose = False):

	'''
	Args:
		fast5: 
			signal file with tombo output (a transcript reference needed)
		junction_pos: 
			0-based position of the splicing site of each transcript reference
		window: 
			signals assigned to a window around the junction_pos will be fetched.
		start_pos:
			start position (0-based) of the junction sequence
		end_pos:
			end position (0-based) of the junction sequence
		verbose: bool, verbose mode.
	Returns: 
		read_raw_signal: <numpy array>

	'''
	
	if start_pos and start_pos:
		junction_bases_start, junction_bases_end = start_pos, end_pos
		window = end_pos - start_pos
	elif junction_pos and window:
		junction_bases_start = junction_pos - int(np.floor(window/2))
		junction_bases_end = junction_pos + int(np.ceil(window/2))
	else:
		print("Missing valid arguments in function get_junction_signal_by_pos")
		sys.exit(0)

	
	fast5object = Fast5Class(fast5)
	fast5object.get_alignment()
	if not fast5object.mapped:
		print("Alignment doesn't exist in fast5!!")
		sys.exit(0)  
		
	assert fast5object.mapped['mapped_end'] - \
		fast5object.mapped['mapped_start'] > window, \
		 "window size larger than the read"
   
		
	# adjust junction seq start pos inside aligned region
	if junction_bases_start >=0 and junction_bases_end >= 0:
		junction_bases_start -= fast5object.mapped['mapped_start']
		junction_bases_end -= fast5object.mapped['mapped_start']
	# the positions are recorded in negative when
	
	
	## transcript sequence is consistent to reverse strand of genome ref    
	#elif junction_bases_start <=0 and junction_bases_end <= 0:
	#	junction_bases_start = fast5object.mapped['mapped_end'] - \
	#		abs(end_pos)  - fast5object.mapped['mapped_start']
	#	junction_bases_end = fast5object.mapped['mapped_end'] - \
	#		abs(start_pos) - fast5object.mapped['mapped_start']

	if junction_bases_start < 0 or \
		junction_bases_end > fast5object.mapped['mapped_end'] - \
		fast5object.mapped['mapped_start']:
		print("Warning: Read discarded! \
		Junction pos is too close to the start of the mapped region")
		return []
	
	#print(junction_bases_start,junction_bases_end)


	if fast5object.mapped['mapped_strand'] == '-':
		
		junction_signal_start = fast5object.read_start_rel_to_raw +\
								fast5object.events[::-1][junction_bases_end][2]
	
		junction_signal_end = fast5object.read_start_rel_to_raw + \
			fast5object.events[::-1][junction_bases_start+1][2] + \
			fast5object.events[::-1][junction_bases_start][3]
	
	if fast5object.mapped['mapped_strand'] == "+":
		junction_signal_start = fast5object.read_start_rel_to_raw +\
			fast5object.events[junction_bases_start][2]
	
		junction_signal_end = fast5object.read_start_rel_to_raw + \
			fast5object.events[junction_bases_end-1][2] + \
			fast5object.events[junction_bases_start][3]
	if verbose:
		print("Corresponding signal start and end:")
		print(junction_signal_start, junction_signal_end)
	
	return fast5object.get_signal(junction_signal_start, junction_signal_end)
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
			sys.exit(0)
		
		mapped_start = h5_f[path]["Alignment"].attrs['mapped_start']
		mapped_end = h5_f[path]["Alignment"].attrs['mapped_end']
		assert mapped_end - mapped_start > window,\
								"window size larger than the read"

		strand = h5_f[path]["Alignment"].attrs['mapped_strand']

		
		mapping_info["start"]=mapped_start
		mapping_info["end"]=mapped_end
		mapping_info["strand"]=strand
	
	return mapping_info

# scrappie squiggle
def sequence_to_squiggle(seq, trim = 0, model = 'squiggle_r94'):
    '''input:
        seq:
			<str> sequence to be converted to squiggle
		trim:
			<int> the number of events that will be ignored in each side.
		model:
			<str> scrappy model name:	{'squiggle_r94',
										'squiggle_r94_rna',
										'squiggle_r10'}
        
		output:
			numpy array: [[mean, std, dwell]...]
    '''
    simulated_seq = scrappy.sequence_to_squiggle(seq,model = model, \
		rescale =True).data(as_numpy = True, sloika = False)
    if trim:
        return simulated_seq[trim:-trim]
    else:
        return simulated_seq

def expect_squiggle_dict(seqs, trim = 0, model = 'squiggle_r94'):
	'''
	read squiggle data from scrappie, ready for dtw
	Args:
		seqs: list of sequence motifs <list of strings>
	Returns: 
		python dictionary of simulated squiggle by scrappie squiggle
	Raises:

	'''
	
	if seqs:
		expect_squiggle_dic = defaultdict(list)
		for seq in seqs:
			squiggle = sequence_to_squiggle(seq = seq, trim = trim, model = model)
			for mean, std, dwell_time in squiggle:
				expect_squiggle_dic[seq] += [[mean, std]] *int(round(dwell_time))
	else:
		print("No valid input detected when generating expect squiggle")
		sys.exit(0)
	return expect_squiggle_dic

def parse_candidate_file(filename):
	'''
	Args:
		filename: Candidate file generated by Candidates_from_gtf.py
	Returns:
		candidate_list: list of <candidate class>
	'''
	class candidate(object):
		def __init__(self, sequences, start, end):
			self.sequences = sequences
			self.start = int(start)
			self.end = int(end)
			
	with open(filename, 'r') as f:
		candidate_list = []
		for line in f:
			line = line.strip().split(',')
			try:
				candidate_list.append(candidate(line[:-3], line[-3], line[-2]))
			except:
				sys.exit(0)
	return candidate_list

def plot_dtw_alignment( long_seq, short_seq, dtw_path, dtw_score = None, \
	show_sd = True, figure_name = "Untitled", \
	figure_title = "Untitled",**plot_args):
	'''
	Args:
		figure_name: <string>
			the figure name that will be saved as
		figure_title: "string"
			figure title
		long_seq: <list or np.array>
			the long sequence in the dtw alignment
		short_seq: 2D np.array [[mean, sd]..]
			the short sequence in the dtw alignment
		dtw_path: <list or np.array>
			the best path of the alignment
		dtw_score: <INT>
			alignment score
		show_sd: True or False
			whether or not plot the sd in dash line
		**plot_args:
			arguments for matplotlib.pyplot
	Returns:
		DTW alignment visulisation.
	'''

	plt.figure(**plot_args)
	plt.plot(long_seq)
	path = np.array(dtw_path)
	plt.plot(path[:,1]-1, short_seq[[path[:,0]-1]][:,0],'g')
	if show_sd:
		plt.plot(path[:,1]-1, short_seq[[path[:,0]-1]][:,0]
						 + short_seq[[path[:,0]-1]][:,1],'g--')
		plt.plot(path[:,1]-1, short_seq[[path[:,0]-1]][:,0]
						 - short_seq[[path[:,0]-1]][:,1],'g--')

	add_info = "\nDist: {:.2f}, path length: {},  \
	Adjusted dist: {:.2f}".format(dtw_score,len(path), \
		dtw_score/len(path)) if dtw_score else ""
	
	plt.title(figure_title + add_info, fontsize=20)
	
	
	i = -1
	while os.path.exists(figure_name + ".png"):
		i += 1
	if i >-1:
		figure_name += str(i)
		
	plt.savefig(figure_name + ".png")
	plt.close()


class Fast5Class(object):
	def __init__(self, filename):
		self.filename = filename
		

	def get_read_id(self):
		with h5py.File(self.filename, 'r') as h5_f: 
			read_key = list(h5_f["Raw/Reads/"].keys())[0]
			read_id = h5_f["Raw/Reads/"][read_key].attrs['read_id']
		return read_id
		
	def get_signal(self, start = None, end = None, normalization = True,\
	 rm_outlier = True):

		start = int(start) if start else 0
		end = int(end) if end else -1

		if end <= start:
			print("InputError: Invalid start and end position when fetching \
				when fetching the raw signals")
			sys.exit(0)
		


		with h5py.File(self.filename, 'r') as h5_f: 
			read_key = list(h5_f["Raw/Reads/"].keys())[0]
			signal = list(h5_f["Raw/Reads/"][read_key]["Signal"]) 

		if normalization:
			self.get_alignment()
			signal = (signal - self.norm_shift)/self.norm_scale
		if rm_outlier:
			signal = signal[start:end]
			signal = self.remove_outlier(signal, thresh = 3)
			return(signal)
 
		return(signal[start:end])

	def get_alignment(self, output=None):
		with h5py.File(self.filename, 'r') as h5_f:
			path = "Analyses/"
			subpath = list(h5_f[path].keys())
			
			for i in subpath:
				if "RawGenomeCorrected" in i:
					path = path + i +'/BaseCalled_template/'

			try:
				self.mapped = dict(h5_f[path]["Alignment"].attrs)
				self.events = np.array(h5_f[path]["Events"])
				self.read_start_rel_to_raw = \
					h5_f[path]["Events"].attrs['read_start_rel_to_raw']
				self.norm_shift = h5_f[path].attrs['shift']
				self.norm_scale = h5_f[path].attrs['scale']*1.4826
			except:
				print("Alignment doesn't exist in fast5!!")
				self.mapped = False
			
			if output == "mapping_info":
				return self.mapped
			elif output == "events":
				return self.events
			elif output == "norm_params":
				return self.norm_shift, self.norm_scale
			else:
				return None
	def remove_outlier(self, normalised_signal, thresh = 3):
		'''
		remove the data points in signal whose obsolute value reaches the thresh.

		Args:
			<list>/<np.array>Normalised signal
			<int> threshold for outliers
		Returns:
			<np.array>Normalised signal with outlier removed
		
		'''
		normalised_signal = np.array(normalised_signal)
		return normalised_signal[np.abs(normalised_signal < thresh)]	