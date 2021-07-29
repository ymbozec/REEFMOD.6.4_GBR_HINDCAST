# REEFMOD.6.4_GBR_HINDCAST
Individual-based model of coral populations on Australia's Great Barrier Reef. Model scripts to simulate the reconstruction of coral trajectories between 2008-2020 published in: Bozec, Y.-M., K. Hock, R. A. Mason, M. E. Baird, C. Castro-Sanguino, S. A. Condie, M. Puotinen, A. Thompson, and P. J. Mumby. 2021. Cumulative impacts across Australia’s Great Barrier Reef: A mechanistic evaluation. Ecological Monographs.


%__________________________________________________________________________
%
% REEFMOD-GBR INSTRUCTIONS (09/2020)
%
% Yves-Marie Bozec,
% School of Biological Sciences
% The University of Queensland
% y.bozec@uq.edu.au
%__________________________________________________________________________


(requires MATLAB R2017)


To run the model and reproduce by simulation the 2008-2020 reconstruction of 3,806 reefs of the GBR:

-> simply run the front script 'MAIN_REEFMOD_GBR.m' (~30 min per run)

Number of replicate runs is currently set to 40.

The front script launches other scripts/data files which are in the folders 'settings', 'functions', 'data'.
Parameterization can be changed in f_multiple_reefs.m' and 'settings_GBR'.
PARAMETERS.m gives the default parameterization. 

This produces an output file 'R0_HINDCAST_GBR.mat' (~500MB).
