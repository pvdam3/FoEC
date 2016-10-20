# FoEC
This script will take a folder with genome fasta files, find mimps and mimp terminal inverted repeats and try to identify candidate effectors. These will be clustered into families and then BLASTed against each of the genomes to identify presence (1) or absence (0). These binary patterns will be hierarchically clustered in R to produce a clustering figure using "heatmap.3.R".

## Concept
The following steps are executed in this pipeline (graphical overview available in [pipeline_overview.pdf](pipeline_overview.pdf)):

**1. Candidate effector identification in each of the provided genome fasta files:**
  * Miniature Impala (mimp) terminal inverted repeat (TIR) identification based on regular expression of the consensus sequence of this repeat.
  * Parsing a sequence (default 2500bp) downstream of this TIR.
  * Finding possible Open Reading Frames (ORFs) within this sequence using two methods; i) translating the sequence in three frames and finding the first Methionine (M) residue followed by a sequence of threshold length and a STOP codon or end of contig, or ii) using [AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/) gene prediction. 
  * These translated sequences are feeded to SignalP to identify potentially secreted proteins.
  * The records that pass this criterium are saved.
**2. Duplicate effector candidates are removed:**
  * All records are combined into a concatenated fastafile containing all the identified sequences from step 1. This file will probably contain many duplicates.
  * A BLAST database is created from this file and each of the fasta records inside the file are BLASTed against this database. This creates a network of 'gene families', thereby essentially marking redundancy.
  * The longest record from each gene family is extracted and saved to a new fastafile, which should contain far fewer records due to the duplicate removal step.
**3. Identifying presence-absence patterns in the genomes:**
  * The list of candidate effectors obtained in step 2 is used as a set of query sequences for BLASTing against a database of each of the genome fasta files. This will result (using some threshold values) in a binary presence or absence of each of the individual effector candidates in each of the genomes.
  * These binary values are stored in a table (.txt) with the effectors on one axis and the genomes on the other.
**4. Hierarchical clustering of binary effector presence patterns:**
  * To discover which genomes are most alike in terms of effector pallette, binary table is imported in an R script, which applies hierarchical clustering on the rows and columns.
  * The resulting matrix and tree are plotted using a script called heatmap.3.R (available on [GitHub](https://gist.github.com/nachocab/3853004))

## Usage
Please first make sure to have all the dependencies installed (see below).
Usage: 
```bash
python FoEC.py -i [infolder] <options>
```
N.B. make sure to provide the **absolute** path to this folder!

Type `python FoEC.py -h` for a detailed help page including options.

## Dependencies
The pipeline relies a number of different 3rd party programs and libraries:
* [AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/) *(optional)*
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

## Configuration
You can set up the paths and variables necessary for your analysis directly in the [FoEC.py](FoEC.py) file.

####paths
* *blastdatabasedir* (a central location where you store your blast databases)
* *contigprefix* (what string precedes the contig number in the fasta headers of your genome files? Suggestion: "contig_")
* *AUGUSTUS_path* (direct path to your AUGUSTUS binary)
* *AUGUSTUS_CONFIG_path* (eg ../augustus-3.1/config')
* *BLASTbindir* (where are your blast binaries stored?)
* *SignalPpath* (direct path to your SignalP binary)

####mimp-search variables
* *distance_MetStop* (sequence downstream of motif used for ORF prediction; default = 2500)
* *distance_Augustus* (sequence downstream of motif used for ORF prediction; default = 5000)
* *min_prot_len* (in aa; default = 30)
* *max_prot_len* (in aa; default = 600)
* *max_d2m* (max distance between mimp TIR and start-codon; default = 2000)
* *SignalP_threshold* (D-value; default = 0.550)

####blast variables
* *PERC_IDENTITY_THRESH* (threshold percentage used to define if a candidate effector is present or absent, which is calculated as follows: the number of identical and correctly aligned nucleotides divided by the query length; default = 30)
* *BLAST_task* ('blastn' or 'megablast')
* *buildblastdb* (should a new blast db be built for the genome files encountered? (recommended for first time this script is run on a set of genomes) )

####clustering variables:
* *distance_matrix_rows* (row distance matrix; default = 1)
* *clustering_method_rows* (row clustering method; default = average)
* *distance_matrix_cols* (col distance matrix; default = 1)
* *clustering_method_cols* (col clustering method; default = average)
  Please choose from the following distance matrices:
  * 1 = Jaccard index (1901) S3 coefficient of Gower & Legendre s1 = a / (a+b+c)
  * 2 = Simple matching coefficient of Sokal & Michener (1958) S4 coefficient of Gower & Legendre s2 = (a+d) / (a+b+c+d)
  * 3 = Sokal & Sneath(1963) S5 coefficient of Gower & Legendre s3 = a / (a + 2(b + c))
  * 4 = Rogers & Tanimoto (1960) S6 coefficient of Gower & Legendre s4 = (a + d) / (a + 2(b + c) +d)
  * 5 = Dice (1945) or Sorensen (1948) S7 coefficient of Gower & Legendre s5 = 2a / (2a + b + c)
  * 6 = Hamann coefficient S9 index of Gower & Legendre (1986) s6 = (a - (b + c) + d) / (a + b + c + d)
  * 7 = Ochiai (1957) S12 coefficient of Gower & Legendre s7 = a / sqrt((a + b)(a + c))
  * 8 = Sokal & Sneath (1963) S13 coefficient of Gower & Legendre s8 = ad / sqrt((a + b)(a + c)(d + b)(d + c))
  * 9 = Phi of Pearson S14 coefficient of Gower & Legendre s9 = (ad - bc) / sqrt((a + b)(a + c)(d + b)(d + c))
  * 10 = S2 coefficient of Gower & Legendre S10 = a / (a + b + c + d)

## References
van Dam, P., Fokkens, L., Schmidt, S. M., Linmans, J. H. J., Kistler, H. C., Ma, L.-J., & Rep, M. (2016). Effector profiles distinguish <I>formae speciales </I>of <I>Fusarium oxysporum</I>. Environmental Microbiology. http://doi.org/10.1111/1462-2920.13445