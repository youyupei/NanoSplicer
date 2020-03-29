'''
apply dynamic time warping to 2 time series
'''

import numpy as np
from numpy import exp, pi, sqrt, log
from scipy.stats import poisson

import matplotlib.pyplot as plt


# require fastdtw to be installed 
# "pip install fastdtw"
#from fastdtw import fastdtw
''' 
    fastdtw(x, y, radius=1, dist=None)      

    return the approximate distance between 2 time series with O(N)
        time and memory complexity
        Parameters
        ----------
        x : array_like
            input array 1
        y : array_like
            input array 2
        radius : int
            size of neighborhood when expanding the path. A higher value will
            increase the accuracy of the calculation but also increase time
            and memory consumption. A radius equal to the size of x and y will
            yield an exact dynamic time warping calculation.
        dist : function or int
            The method for calculating the distance between x[i] and y[j]. If
            dist is an int of value p > 0, then the p-norm will be used. If
            dist is a function then dist(x[i], y[j]) will be used. If dist is
            None then abs(x[i] - y[j]) will be used.
        Returns
        -------
        distance : float
            the approximate distance between the 2 time series
        path : list
            list of indexes for the inputs x and y
'''

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
        
    # f(z)
        return 0.9189385 +  z**2/2
    # f(x|mean,st)
        #return 0.9189385 +  z**2/2 + log(b_std)
    

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

def dtw_local_alignment_max_mean(long, short, radius = None, dist_type = None, \
         upper = np.inf):

    short_len = len(short)
    long_len = len(long)
    mean_matrix = np.zeros((short_len + 1, long_len + 1, 2))
    mean_matrix[1:, 0, 0] = np.inf
    
    pre_step_matrix = np.zeros((short_len + 1, long_len + 1), dtype = int)
    '''
    matrix for recording the best path. Each of the cells contains one of three
    possible integer: 0 ,1 ,2 indecating the corresponding value in mean_matrix
    (mean_matrix[i,j]) came from the mean_matrix[i-1,j], mean_matrix[i - 1,j - 1],
     and mean_matrix[i, j - 1] respectively.
    '''
    
    for i in range(1, short_len + 1):
        for j in range(1, long_len + 1):

            pre_values = np.array((mean_matrix[i-1,j] , 
                        mean_matrix[i - 1,j - 1],
                        mean_matrix[i, j - 1]))
            
            
            #pre_step_matrix[i, j] = np.argmin(pre_values)
            min_index = np.random.choice(\
                            np.where(pre_values[:,0]==min(pre_values[:,0]))[0])
            
            pre_step_matrix[i, j] = min_index
            pre_mean ,pre_p_length = pre_values[min_index]

            
            mean_matrix[i, j] = (pre_mean * pre_p_length +\
                cost(long[j -1], *short[i - 1], dist_type = dist_type))/\
                (pre_p_length + 1), pre_p_length + 1


    best_score = min(mean_matrix[-1,:,0])
    best_path = []
  
    traced_short_index = short_len
    #traced_long_index = np.argmin(mean_matrix[-1:])
    traced_long_index = np.random.choice(\
                np.where(mean_matrix[-1,:,0]==min(mean_matrix[-1,:,0]))[0])
    
    if False:
        plt.plot(mean_matrix[-1, :, 0])
        plt.savefig('path_score.png')

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

    return best_path[::-1], best_score,mean_matrix[:,:,0]




    best_score = min(mean_matrix[-1,:,0])
    best_path = []
  
    traced_short_index = short_len
    #traced_long_index = np.argmin(mean_matrix[-1:])
    traced_long_index = np.random.choice(\
                np.where(mean_matrix[-1,:,0]==min(mean_matrix[-1,:,0]))[0])
    
    if False:
        plt.plot(mean_matrix[-1, :, 0])
        plt.savefig('path_score.png')

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

    return best_path[::-1], best_score,mean_matrix[:,:,0]



def dtw_global_alignment_max_sum(seq1, seq2, radius = None, \
            dist_type = None, upper = np.inf):


    # init accumulated matrix
    seq2_len = len(seq2)
    seq1_len = len(seq1)
    cum_matrix = np.zeros((seq2_len + 1, seq1_len + 1))
    cum_matrix[1:, 0] = np.inf
    cum_matrix[0, 1:] = np.inf
    
    pre_step_matrix = np.zeros((seq2_len + 1, seq1_len + 1), dtype = int)
    '''
    matrix for recording the best path. Each of the cells contains one of three
    possible integer: 0 ,1 ,2 indecating the corresponding value in cum_matrix
    (cum_matrix[i,j]) came from the cum_matrix[i-1,j], cum_matrix[i - 1,j - 1],
     and cum_matrix[i, j - 1] respectively.
    '''
    
    # build matrix
    for i in range(1, seq2_len + 1):
        for j in range(1, seq1_len + 1):

            pre_values = (cum_matrix[i-1,j], 
                        cum_matrix[i - 1,j - 1],
                        cum_matrix[i, j - 1])
            
            
            #pre_step_matrix[i, j] = np.argmin(pre_values)
            pre_step_matrix[i, j] = np.random.choice(\
                            np.where(pre_values==min(pre_values))[0])
            
            cum_matrix[i, j] = min(pre_values) +\
                        cost(seq1[j -1], *seq2[i - 1], dist_type = dist_type)
    best_score = cum_matrix[-1,-1]
    best_path = []
  

    # init trace back index
    traced_seq1_index = seq1_len
    traced_seq2_index = seq2_len

    


    #plt.plot(cum_matrix[-1:][0])
    #plt.savefig('path_score.png')

    while True:
        best_path.append([traced_seq2_index, traced_seq1_index])
        pre_step = pre_step_matrix[traced_seq2_index, traced_seq1_index]

        if traced_seq2_index == 1:
            break

        if pre_step in (0, 1):
            traced_seq2_index -= 1
        
        if pre_step in (1, 2):
            traced_seq1_index -= 1
    # best_path: 0-based coordinate on the (i+1)*(j+1) matrix
    best_path = np.array(best_path)
    return best_path[::-1], best_score/len(best_path[::-1]),cum_matrix

def dtw_local_alignment_max_sum_band(long, short, band_prop = 0.4, \
            dist_type = None, upper = np.inf):

    short_len = len(short)
    long_len = len(long)
    band_len = int(np.ceil(long_len* band_prop))
    band_move = (long_len - band_len)/short_len

    cum_matrix = np.full((short_len + 1, long_len + 1), np.inf)
    cum_matrix[0, 0:band_len+1] = 0
    
    pre_step_matrix = np.zeros((short_len + 1, long_len + 1), dtype = int)
    '''
    matrix for recording the best path. Each of the cells contains one of three
    possible integer: 0 ,1 ,2 indecating the corresponding value in cum_matrix
    (cum_matrix[i,j]) came from the cum_matrix[i-1,j], cum_matrix[i - 1,j - 1],
     and cum_matrix[i, j - 1] respectively.
    '''
    
    for i in range(1, short_len + 1):
        j_start = int((i - 1) * band_move)
        for j in range(j_start, j_start + band_len):

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



def dtw_local_alignment_max_sum_band_flipped(long, short, band_prop = 0.4, \
            dist_type = None, upper = np.inf):
    


    #flip long and short
    long, short = short, long
    short_len = len(short)
    long_len = len(long)
    band_len = int(np.ceil(long_len* band_prop))
    band_move = (long_len - band_len)/short_len

    cum_matrix = np.full((short_len + 1, long_len + 1), np.inf)
    cum_matrix[0, 0:band_len+1] = 0
    
    pre_step_matrix = np.zeros((short_len + 1, long_len + 1), dtype = int)
    '''
    matrix for recording the best path. Each of the cells contains one of three
    possible integer: 0 ,1 ,2 indecating the corresponding value in cum_matrix
    (cum_matrix[i,j]) came from the cum_matrix[i-1,j], cum_matrix[i - 1,j - 1],
     and cum_matrix[i, j - 1] respectively.
    '''
    
    for i in range(1, short_len + 1):

        j_start = int((i - 1) * band_move)
        for j in range(j_start, j_start + band_len):

            pre_values = (cum_matrix[i-1,j], 
                        cum_matrix[i - 1,j - 1],
                        cum_matrix[i, j - 1])
            
            
            #pre_step_matrix[i, j] = np.argmin(pre_values)
            pre_step_matrix[i, j] = np.random.choice(\
                            np.where(pre_values==min(pre_values))[0])

            cum_matrix[i, j] = min(pre_values) +\
                        cost(short[i -1], *long[j - 1], dist_type = dist_type)
    best_score = min(cum_matrix[-1,:])
    best_path = []
  
    traced_short_index = short_len
    #traced_long_index = np.argmin(cum_matrix[-1:])
    traced_long_index = np.random.choice(\
                np.where(cum_matrix[-1,:]==min(cum_matrix[-1,:]))[0])

    while traced_short_index > 0:

        best_path.append([traced_short_index, traced_long_index])
        pre_step = pre_step_matrix[traced_short_index, traced_long_index]

        if pre_step in (0, 1):
            traced_short_index -= 1
        
        if pre_step in (1, 2):
            traced_long_index -= 1
    # best_path: 0-based coordinate on the (i+1)*(j+1) matrix
    best_path = np.array(best_path)
    return best_path[::-1], best_score/len(best_path[::-1]),cum_matrix

def dist_to_likelihood(long, short, path, dist_type = None):
    pre_long = path[0][1]# first long index
    num_of_match = 0
    accu_log_density = 0
    log_likelihood = []
    for short_index, long_index in path:
        if long_index == pre_long:
            num_of_match += 1
            accu_log_density += cost(long[long_index -1], *short[short_index - 1], dist_type = dist_type)
        else:
            log_likelihood.append(accu_log_density/num_of_match)
            pre_long = long_index
            num_of_match = 1
            accu_log_density = cost(long[long_index -1], *short[short_index - 1], dist_type = dist_type)
    log_likelihood.append(accu_log_density/num_of_match)
    return sum(log_likelihood), len(long) - len(log_likelihood)

def dist_to_likelihood_flipped(long, short, path, dist_type = None):
    pre_short = path[0][0]# first short index
    log_density_same_x = []
    log_likelihood = []
    for short_index, long_index in path:
        if short_index == pre_short:
            log_density_same_x.append(cost(short[short_index - 1], *long[long_index -1], dist_type = dist_type))
        else:
            log_likelihood.append(max(log_density_same_x))
            pre_short = short_index
            log_density_same_x = [cost(short[short_index - 1], *long[long_index -1], dist_type = dist_type)]
    log_likelihood.append(max(log_density_same_x))
    return sum(log_likelihood), len(long) - len(log_likelihood)

def dist_to_likelihood_flipped_new_path(long, short, path, dist_type = None):
    pre_short = path[0][0]
    log_density_list = []
    multi_long = [path[0][1]]
    log_likelihood = 0
    new_path = []
    for short_index, long_index in path:
        if short_index == pre_short:
            multi_long.append(long_index)
            log_density_list.append(cost(short[short_index - 1], *long[long_index -1], dist_type = dist_type))
        else:
            log_likelihood += min(log_density_list)
            new_path.append([short_index -1, multi_long[log_density_list.index(min(log_density_list))]])
            multi_long = [long_index]
            pre_short = short_index
            log_density_list = [cost(short[short_index - 1], *long[long_index -1], dist_type = dist_type)]
    log_likelihood += min(log_density_list)
    new_path.append([short_index, multi_long[log_density_list.index(min(log_density_list))]])
    return new_path, log_likelihood



def dist_to_likelihood_flipped_time_serie(long, short, path, dist_type = None):
    num_of_match = 0
    pre_long_index = path[0][1]# first long index
    acc_poi_log_dens = 0
    for short_index, long_index in path:
        print(long[long_index-1],long[pre_long_index-1])
        if (long[long_index-1] == long[pre_long_index-1]).all():
            num_of_match += 1
        else:
            print(num_of_match,long_index,pre_long_index)
            acc_poi_log_dens += poisson.logpmf(num_of_match,long_index-pre_long_index)
            pre_long_index = long_index
            num_of_match = 1
    print(num_of_match,path[-1][1], pre_long_index)
    acc_poi_log_dens += poisson.logpmf(num_of_match,path[-1][1]-pre_long_index + 1)
    return acc_poi_log_dens

def main():
    short = np.array([[1,1],[1,1],[1,1],[3,1],[4,1],[5,1],[5,1],[5,1],[7,1],[5,1],[5,1],[4,1],[5,1]])
    long = np.array([1,1,1,3,4,5,7,5,4,5])


    path , score = dtw_local_alignment_max_sum(long, short)
    path , score = dtw_local_alignment_max_sum_band_flipped(long, short)
    plt.plot(long)
    plt.plot(np.array(path)[:,1]-1, short[[np.array(path)[:,0]-1]][:,0])
    plt.savefig('test.png')
    print(path)
if __name__ == '__main__':
    main()

    
            
