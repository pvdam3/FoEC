# Fo_effector_clustering (FoEC) synopsis
This script will take a folder with genome fasta files, find mimps and mimp terminal inverted repeats and try to identify candidate effectors. These will be clustered into families and then BLASTed against each of the genomes to identify presence (1) or absence (0). These binary patterns will be hierarchically clustered in R to produce a clustering figure using "heatmap.3.R". 

# Code Example
Usage: python FoEC.py -i [infolder] <options>

Options:
  -h, --help            show this help message and exit
  -i INFOLDER, --in=INFOLDER
                        provide a folder with genome fasta files
  -o OUTFOLDER, --out=OUTFOLDER
                        output folder; default is: /Users/Peter/Programming/py
                        thon/01.Fo_effector_clustering/output/output_16.10.20_
                        16h44m12
  -e EFFECTOR_FASTA, --effector_fasta=EFFECTOR_FASTA
                        Skip the effector prediction pipeline and go straight
                        to blasting and clustering using a self-supplied
                        effector list in fasta format.
  -m, --met_stop_prediction_only
                        Do not use Augustus gene prediction, only identify
                        ORFs by looking for first occurrence of a Methionine
                        (met) to first Stop codon (stop).
  -u, --use_nonredundant_effectors
                        In case multiple ORF prediction methods are used (NOT
                        using -m), the outputs from these methods can be
                        clustered first, then concatenated, then clustered
                        (use this -u option).
                        Default behaviour is not to use this option, meaning
                        taking the output from both methods, concatenate them
                        and then cluster them into gene families.
