# FoEC
This script will take a folder with genome fasta files, find mimps and mimp terminal inverted repeats and try to identify candidate effectors. These will be clustered into families and then BLASTed against each of the genomes to identify presence (1) or absence (0). These binary patterns will be hierarchically clustered in R to produce a clustering figure using "heatmap.3.R".

## Concept
The following steps are executed in this pipeline (graphical overview available in [pipeline_overview.pdf](pipeline_overview.pdf)):

1. Candidate effector identification in each of the provided genome fasta files:
  * Miniature Impala (mimp) terminal inverted repeat (TIR) identification based on regular expression of the consensus sequence of this repeat.
  * Parsing a sequence (default 2500bp) downstream of this TIR.
  * Finding possible Open Reading Frames (ORFs) within this sequence using two methods; i) translating the sequence in three frames and finding the first Methionine (M) residue followed by a sequence of threshold length and a STOP codon or end of contig, or ii) using [AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/) gene prediction. 
  * These translated sequences are feeded to SignalP to identify potentially secreted proteins.
  * The records that pass this criterium are saved.
2. Duplicate effector candidates are removed:
  * All records are combined into a concatenated fastafile containing all the identified sequences from step 1. This file will probably contain many duplicates.
  * A BLAST database is created from this file and each of the fasta records inside the file are BLASTed against this database. This creates a network of 'gene families', thereby essentially marking redundancy.
  * The longest record from each gene family is extracted and saved to a new fastafile, which should contain far fewer records due to the duplicate removal step.
3. Identifying presence-absence patterns in the genomes:
  * The list of candidate effectors obtained in step 2 is used as a set of query sequences for BLASTing against a database of each of the genome fasta files. This will result (using some threshold values) in a binary presence or absence of each of the individual effector candidates in each of the genomes.
  * These binary values are stored in a table (.txt) with the effectors on one axis and the genomes on the other.
4. Hierarchical clustering of binary effector presence patterns:
  * To discover which genomes are most alike in terms of effector pallette, binary table is imported in an R script, which applies hierarchical clustering on the rows and columns.
  * The resulting matrix and tree are plotted using a script called heatmap.3.R (available on [GitHub](https://gist.github.com/nachocab/3853004))

## Usage
Please first make sure to have all the dependencies installed (see below).
Usage: 
```bash
python FoEC.py -i [infolder] <options>
```
(N.B. make sure to provide the **absolute** path to this folder!)
Type `python FoEC.py -h` for a detailed help page including options.

## Dependencies
The pipeline relies a number of different 3rd party programs and libraries:
* *(optional)*[AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/)  
* [SignalP](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp)
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* R with the following libraries installed:
  * dendextend
  * gplots
  * ctc
  * extrafont
  * ade4
* Python with the following package installed:
  * [BioPython](http://biopython.org/wiki/Download)

## References
van Dam, P., Fokkens, L., Schmidt, S. M., Linmans, J. H. J., Kistler, H. C., Ma, L.-J., & Rep, M. (2016). Effector profiles distinguish <I>formae speciales </I>of <I>Fusarium oxysporum</I>. Environmental Microbiology. http://doi.org/10.1111/1462-2920.13445