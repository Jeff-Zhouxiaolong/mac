The following files are available in the package:

nemd.m: Computes multivariate emd for an input signal. See the comments
in the actual nemd.m to find out the usage of this function

The following .mat files contain dataset which may serve as an input to
nemd.m. Just load any of these files in matlab and send the resulting output 
vector as an input to the function nemd.

1) syn_12channel_inp.mat: contains synthetically generated 12 channel
data set (with combination of 5 tones (sinewaves) and noise added to few channels).

2) syn_16channel_inp.mat: contains synthetically generated 16 channel
data set (combined 6 tones (sinewaves) and noise added to some channels)

3) syn_hex_inp.mat: contains synthetically generated 6 channel data set
(combined 4 sinewaves and noise added to some channels - see the
Multivariate EMD paper and the Supplementary Material for more detaiil)

4) taichi_hex_inp.mat: contains hexavariate real world taichi dataset.
(two 3D recordings from intertial bodysensors combined into a single
hexavariate signal [left wrist and left ankle])

