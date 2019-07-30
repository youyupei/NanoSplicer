'''
For looking at the structure inside the fast5 file
python version on spartan:
    Python/2.7.13-intel-2017.u2

Usage:
script.py [HDF5_FILENAME]

'''
import h5py, sys, os, getopt
from collections import defaultdict
import pickle


def parse_arg(argv):
    """ 
    parsing flagged command line argument
    return inputfile then outputfile
    """
    try: 
        opts, args = getopt.getopt(argv,"hp:o:v",
                    ["help","path=","output=","verbose"])
    except getopt.GetoptError:
        print('Usage: {} -p <inputpath> -o <outputfile>'.format(sys.argv[0]))
        sys.exit()
    
    #default setting:
    verbose = False
    
    for opt, arg in opts:
        if opt in ('-h', "--help"):
            print('Usage: {} -p <inputpath> -o <outputfile>'.format(sys.argv[0]))
            sys.exit()
        elif opt in ("-p", "--path"):
            inputpath = arg
        elif opt in ("-o", "--output"):
            outputfile = arg
        elif opt in ("-v", "--verbose"):
            verbose=True   
    
    return inputpath, outputfile, verbose

def print_hdf5_structure(hdf5_f):
    # print the structure of a opened hdf5 file
    print("Attributes:")
    print(hdf5_f.attrs.keys())
    for key in hdf5_f.keys():
        print(key)
        print(hdf5_f[key])
        print_hdf5_structure(hdf5_f[key])

# get read id given the fast5 filename
def read_id(filename):
    hdf5_f = h5py.File(filename, 'r')
    try:
        read_id = hdf5_f["Raw/Reads/"].values()[0].attrs["read_id"]
    except:
        print("Error: Fail to read id from {}".format(filename))
        return(None)
    
    hdf5_f.close()
    return(read_id)

def read_id_dic_generator(PATH):
    # walk through all file in PATH
    # read the Read ID and return a dictionary of key: readID value: file_path
    
    ReadId_to_SquiggleFile_dict = defaultdict(list)
    try:
        count = 0
        for root, dirs, files in os.walk(PATH, topdown=False):
            fast5s = [os.path.join(root, name) for name in files]
            
            # extract the read ID          
            for fast5 in fast5s:
                ReadId_to_SquiggleFile_dict[read_id(fast5)].append(fast5)
                count += 1
                print("{} files finished".format(count))
        print("Dictionary creation successful")
        return ReadId_to_SquiggleFile_dict

    except:
        print("Error: Creating Dictionary Failed!!!!!")
        exit()
    
    
def pickle_dump(py_object, output_filename):
    # Pickle the object to <output.pickle> 
    # using the highest protocol available.
    with open(output_filename, 'wb') as f:
        pickle.dump(py_object, f, pickle.HIGHEST_PROTOCOL)

def pickle_load(pickle_file):
    #load pickle object
    with open(pickle_file, "rb") as f:
        try:
            return pickle.load(f)
        except EOFError:
            print("pickle.load failed!!!!!!!")        

if __name__ == "__main__": 
    inputpath, outputfile, verbose = parse_arg(sys.argv[1:])
    read_id_dic = read_id_dic_generator(inputpath)
    pickle_dump(read_id_dic, outputfile)