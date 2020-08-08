# T0_Tagged_Tracks_Scripts
This is a repository of scripts separate from the UBOONE_SUITE (the MicroBooNE standard repository that operates within LArSoft, a framework for Liquid Argon Time Projection Chamber experiments).  They tag a specific type of reconstructed particle track, and anode-piercing/cathode-piercing (ACPT) cosmic ray muon track, in the detector. Each of the files fulfills the following purpose(s):

SpaceChargeStudyWithSpacePoints.cpp - This C++ script uses the 'gallery' package on the Fermilab computers to extract, identify, and characterize an ACPT cosmic ray muon track from MicroBooNE artroot files.
mcc83_tracks_top_anode_plots.py - This python script uses pandas and matplotlib to plot the path within the MicroBooNE detector of top-piercing, anode-piercing ACPT cosmic ray muon tracks.
