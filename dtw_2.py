'''
apply dynamic time warping to 2 time series
'''

import numpy as np
import matplotlib.pyplot as plt
from math import exp, pi, sqrt, log

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



def dtw_local_alignment(long, short, radius = None, dist_type = None, upper = np.inf):

    def cost(x, *y, dist_type = None, upper = upper):
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

    return best_path, best_score

def main():
    short = np.array([[1,4],[2,9],[3,1],[4,1],[5,2]])
    long = np.array([1,1,1,3,4,5,7,5,4,5])

    if False:
        path , score = dtw_local_alignment(long, short)
        plt.plot(long)
        plt.plot(np.array(path)[:,1]-1, short[[np.array(path)[:,0]-1]][:,0])
        plt.savefig('test.png')
        print(path)
if __name__ == '__main__':
    main()

    
            





