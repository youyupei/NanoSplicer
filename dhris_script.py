import sys
from collections import defaultdict
import scrappy
import numpy as np
from numpy import exp, pi, sqrt, log

#dtw
def cost(x, *y, dist_type = None, upper = np.inf):
    '''
    input:
        two float when dist_type = "manhattan"
        three float: x, y_mean, y_std when dist_type = "z_score"
    '''
    def manhattan(a, b, upper = upper):
        return abs(a - b)
    
    def z_score(a, b_mean, b_std, upper = upper):
        diff = min(abs(a-b_mean), upper)
        return diff/b_std

    def log_likelihood(a, b_mean, b_std, upper = upper):
        diff = min(abs(a-b_mean), upper)
        z = diff/b_std
        return -(log(1/sqrt(2*pi)) - z**2/2)
        #return 0.6931472+log(1+z**2)
    

    y_len = len(y)
    if y_len not in (1,2):
        exit("FATAL!: unexpected input in distance matrics.")
    if dist_type == None:
        dist_type = "manhattan" if y_len == 1 else "log_likelihood" 
    if dist_type not in ("manhattan", "z_score", "log_likelihood"):
        exit("FATAL: can't recognise the distance matrics ['manhattan', 'z_score', 'log_likelihood'],")

    if dist_type == "manhattan":
        return manhattan(x, y[0])

    if dist_type == "z_score":
        return z_score(x, y[0], y[1])

    if dist_type == "log_likelihood":
        return log_likelihood(x, y[0], y[1])


def dtw_local_alignment_max_sum(long, short, radius = None, \
            dist_type = None, upper = np.inf):


    
    short_len = len(short)
    long_len = len(long)
    cum_matrix = np.zeros((short_len + 1, long_len + 1))
    cum_matrix[1:, 0] = np.inf
    
    pre_step_matrix = np.zeros((short_len + 1, long_len + 1), dtype = int)
    '''
    matrix for recording the best path. Each of the cells contains one of three
    possible integer: 0 ,1 ,2 indecating the corresponding value in cum_matrix
    (cum_matrix[i,j]) came from the cum_matrix[i-1,j], cum_matrix[i - 1,j - 1],
     and cum_matrix[i, j - 1] respectively.
    '''
    
    for i in range(1, short_len + 1):
        for j in range(1, long_len + 1):

            pre_values = (cum_matrix[i-1,j], 
                        cum_matrix[i - 1,j - 1],
                        cum_matrix[i, j - 1])
            
            
            #pre_step_matrix[i, j] = np.argmin(pre_values)
            pre_step_matrix[i, j] = np.random.choice(\
                            np.where(pre_values==min(pre_values))[0])
            
            cum_matrix[i, j] = min(pre_values) +\
                        cost(long[j -1], *short[i - 1], dist_type = dist_type)
    best_score = min(cum_matrix[-1,:])
    best_path = []
  
    traced_short_index = short_len
    #traced_long_index = np.argmin(cum_matrix[-1:])
    traced_long_index = np.random.choice(\
                np.where(cum_matrix[-1,:]==min(cum_matrix[-1,:]))[0])
    


    #plt.plot(cum_matrix[-1:][0])
    #plt.savefig('path_score.png')

    while True:
        best_path.append([traced_short_index, traced_long_index])
        pre_step = pre_step_matrix[traced_short_index, traced_long_index]

        if traced_short_index == 1:
            break

        if pre_step in (0, 1):
            traced_short_index -= 1
        
        if pre_step in (1, 2):
            traced_long_index -= 1
    # best_path: 0-based coordinate on the (i+1)*(j+1) matrix
    best_path = np.array(best_path)
    return best_path[::-1], best_score/len(best_path[::-1]),cum_matrix

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


def main():
    args = sys.argv
    dtw_local_alignment = dtw_local_alignment_max_sum
    if len(args) < 3:
        print("Usage: python3 {} <squiggle.csv> <candidate.csv>".format(args[0]))
        sys.exit(0)
    squiggle_file, candidate_file = args[1:3]
    output_file = "dtw_score.csv"

    with open(squiggle_file, 'r') as squiggle_f, open(candidate_file,'r') as candidate_f, \
                open(output_file, 'w') as output_f:
        count = 0
        while True:
            squiggle = next(squiggle_f, "end_of_file")
            if squiggle != "end_of_file":
                candidates = next(candidate_f).strip().split(",")[0:2]
                model_dic = expect_squiggle_dict(candidates)
                dtw_long = np.array([float(x) for x in squiggle.strip().split(",")], float)

                for key in candidates:
                    dtw_short = np.array(model_dic[key],float)
                    path , score = dtw_local_alignment(
                        dtw_long, dtw_short, dist_type = "log_likelihood")[0:2]
                    output_f.write(str(score) + ',')

                output_f.write('\n')
                count += 1
                print('{} squiggle finished'.format(count))


            else:
                break

if __name__ == "__main__":
    main()