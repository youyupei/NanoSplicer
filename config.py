# CHROMOSOME_NAME (temp)
CHROMOSOME_NAME = "NC_000001.11"
#CHROMOSOME_NAME = "chrIS"

# wether or not output junction squiggles and candidate squiggles matches as csv
SAVE_DATA = False

# wavelet denoising level (0 if wavelet is not applying)
WAVELET_DENOISE = True
# WAVELET_LEVEL = 3 
# np.ceil(np.log(len(squiggle)))

# minimum data point assigned to a k-mer in dtw
UNIFORM_DWELL = 4
MINIMUM_POINT_FOR_DIST_SEG = 4

# prior ratio for sequence pattern
PRIOR_RATIO = 9

# output
PLOT = False
PLOT_LR = False
RESULT = True