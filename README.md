# Fo_effector_clustering (FoEC) synopsis
This script will take a folder with genome fasta files, find mimps and mimp terminal inverted repeats and try to identify candidate effectors. These will be clustered into families and then BLASTed against each of the genomes to identify presence (1) or absence (0). These binary patterns will be hierarchically clustered in R to produce a clustering figure using "heatmap.3.R". 

# Usage
Usage: `python FoEC.py -i [infolder] <options>` (N.B. make sure to provide the _absolute_ path to this folder!)
Type `python FoEC.py -h` for a detailed help page including options.

# Dependencies
The pipeline used a number of different 3rd party programs including Augustus /(optional)/
