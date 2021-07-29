# REEFMOD.6.4_GBR_HINDCAST


Model scripts to simulate the reconstruction of coral trajectories between 2008-2020. Work published in: Bozec, Y.-M., K. Hock, R. A. Mason, M. E. Baird, C. Castro-Sanguino, S. A. Condie, M. Puotinen, A. Thompson, and P. J. Mumby. 2021. Cumulative impacts across Australia’s Great Barrier Reef: A mechanistic evaluation. Ecological Monographs.

To run the model and simulate the reconstruction of 3,806 reefs of the GBR, simply run the front script MAIN_REEFMOD_GBR.m.

Number of replicate runs is currently set to 40 (~30 min per run).

The front script launches other scripts/data files which are in the folders 'settings', 'functions', 'data'. Parameterization can be changed in f_multiple_reefs.m and settings_GBR.m. PARAMETERS.m gives the default parameterization.

This produces an output file R0_HINDCAST_GBR.mat (~500MB).

(requires MATLAB R2017)

Yves-Marie Bozec (y.bozec@uq.edu.au), School of Biological Sciences, The University of Queensland, Australia.
