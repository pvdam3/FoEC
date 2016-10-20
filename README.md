# FoEC
This script will take a folder with genome fasta files, find mimps and mimp terminal inverted repeats and try to identify candidate effectors. These will be clustered into families and then BLASTed against each of the genomes to identify presence (1) or absence (0). These binary patterns will be hierarchically clustered in R to produce a clustering figure using "heatmap.3.R".

## Steps
The following steps are executed in this pipeline:

1. Miniature Impala (mimp) terminal inverted repeat (TIR) identification based on regular expression of the consensus sequence of this repeat.
  * Test
2. Following item

## Usage
Usage: 
```bash
python FoEC.py -i [infolder] <options>
```
(N.B. make sure to provide the **absolute** path to this folder!)
Type `python FoEC.py -h` for a detailed help page including options.

## Dependencies
The pipeline relies a number of different 3rd party programs and libraries including Augustus *(optional)* and BioPython.

## References
van Dam, P., Fokkens, L., Schmidt, S. M., Linmans, J. H. J., Kistler, H. C., Ma, L.-J., & Rep, M. (2016). Effector profiles distinguish <I>formae speciales </I>of <I>Fusarium oxysporum</I>. Environmental Microbiology. http://doi.org/10.1111/1462-2920.13445