# REEFMOD.6.4 (GBR HINDCAST)

This repository contains the scripts of ReefMod-GBR, a coral individual-based model that simulates coral populations across 2,300 km of Australia's Great Barrier Reef (GBR). This version was used to simulate the reconstruction of coral trajectories across the GBR between 2008-2020, based on temporally- and spatially-explicit forcing of water quality, cyclones, heat stress (coral bleaching) and the simulated population dynamics of the coral-eating (crown-of-thorns) starfish.

Citation: Bozec, Y.-M., K. Hock, R. A. Mason, M. E. Baird, C. Castro-Sanguino, S. A. Condie, M. Puotinen, A. Thompson, and P. J. Mumby. 2021. Cumulative impacts across Australia’s Great Barrier Reef: A mechanistic evaluation. Ecological Monographs.

---------------------

The model was developed in MATLAB (R2017b).

To simulate the 2008-2020 reconstruction of 3,806 reefs of the GBR, simply run the front script MAIN_REEFMOD_GBR.m. This will produce the output file R0_HINDCAST_GBR.mat (~500MB). Number of replicate runs is currently set to 40 (~20 min per run) but can be changed in the front script. 

---------------------

Yves-Marie Bozec (y.bozec@uq.edu.au)

School of Biological Sciences, The University of Queensland, Australia.

