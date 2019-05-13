'''
 Plot the squiggle from fast5 file
'''
import sys, h5py
import numpy as np
import matplotlib.pyplot as plt

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

def read_raw_segnal(fast5, start = None, end = None):
    '''
    take the fast5 filename as the input, return a np array of raw signal
    input:
        start, end: <int>
            index (0-based) for truncating the output signal
    '''
    with h5py.File(fast5, 'r') as h5_f:
        
        try:
            read_key = list(h5_f["Raw/Reads/"].keys())[0]
            signal = h5_f["Raw/Reads/"][read_key]["Signal"]
            
            return(signal[start:end])
        except:
            print("Error: Fail to read {}".format(fast5))
            return(None)

def mean_normalization(signal):
    '''
    input: signal dataset from fast5
    output: nomalised signal <np.array>
    '''
    signal = np.array(signal)
    mad_scale = np.median(abs(signal - np.median(signal)))

    norm_signal = (signal - np.median(signal)) / mad_scale


    return(norm_signal)

def main():
    argv = sys.argv
    if len(argv) < 3:
        exit("Usage: python {} <fast5 filename> <save fig as> <optianl: start> <optional: end>.format(argv[0])")
    filename = argv[1]
    
    if len(argv) == 5:
        signal = read_raw_segnal(filename, int(argv[3]), int(argv[4]))
    else:
        signal = read_raw_segnal(filename)

    # Normalisation
    signal = mean_normalization(signal)




    fig, ax = plt.subplots()
    fig.set_size_inches(20, 9.5)
    ax.plot(np.array(signal))

    ax.set_ylabel("Current levle", fontsize = 25)
    ax.set_xlabel("index", fontsize = 25)
    ax.grid()

    fig.savefig(argv[2])


if __name__ == "__main__":
    main()
    
