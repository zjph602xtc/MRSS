% This file will download required files and folders that my toolbox relies
% on to compute some of the dysphonia measures. The user may wish to take
% this step manually. The author of this toolbox has no responsibility for 
% suggesting the downloading of third party software

% Copyright (c) A. Tsanas, 2014
% modified by Tianchen Xu

unzip('voicebox.zip', 'voicebox')
unzip('pack_emd.zip', 'EMD')
unzip('fastdfa.zip', 'DFA')
unzip('rpde.zip', 'RPDE')
unzip('pda_matlab.zip', 'SHRP'); % if there is any problem, search online for "shrp Matlab" and download the set of functions
