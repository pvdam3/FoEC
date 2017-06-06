# FoEC.py installation instructions for MacOS:


- to install biopython and bcbio-gff:
	* install anaconda (Python 2.7 version) via https://www.continuum.io/downloads
	* type in the terminal:
		conda install biopython
	* type in the terminal:
		pip install bcbio-gff


- install Homebrew (https://brew.sh/) by typing the below command in the terminal and following the instructions:
	/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" 


- install augustus using homebrew by typing the below command in the terminal:
	brew install homebrew/science/augustus


- fill in form at http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp to obtain SignalP academic download by email:
	* platform = Darwin (for MacOS)
	* uncompress and untar tar.gz package --> this will produce a directory called 'signalp-4.1'
		* open signalp file
		* edit paragraph general settings in the top of the signalp file:
				SIGNALP		full path to the signalp-4.1 directory on your system, e.g. '/Applications/signalp-4.1'
				
				Optional to change:
				outputDir				where to store temporary files (writable to all users)
				MAX_ALLOWED_ENTRIES		to limit the number of input sequences allowed per run


- install blast+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/):
	* download ncbi-blast-2.6.0+-x64-macosx.tar.gz


- to install R type the following commands in the terminal:
	* brew install homebrew/science/r
	or
	* brew tap homebrew/science
	* brew install Caskroom/cask/xquartz
	* brew install r


- to install the required R packages type the following commands in the terminal:
	r
	install.packages("dendextend")
		Question: select CRAN mirror --> 7
	install.packages("gplots")
	install.packages("extrafont")
	install.packages("ade4")
	source("https://bioconductor.org/biocLite.R")
	biocLite("ctc")
		 Question: update all/some/none? [a/s/n] --> a

- download effector pipeline "FoEC" from GitHub (pvdam3)
	* change paths in FoEC.py script:
		Examples of paths to be changed:
		blastdatabasedir			= '/Users/Mara/Documents/Genomes/Blastdbs'
		contigprefix 				= 'contig_' 	# default: contig_
		AUGUSTUS_path 				= '/usr/local/bin/augustus'
		AUGUSTUS_CONFIG_path 		= '/usr/local/opt/augustus/libexec/config'
		BLASTbindir 				= '/Applications/ncbi-blast-2.6.0+/bin'
		SignalPpath					= '/Applications/signalp-4.1/signalp'
	* If augustus is installed using homebrew, the config_PATH location is '/usr/local/opt/augustus/libexec/config',
	  you can check this by typing in the terminal:
	  	brew info augustus


- add blast to your path in .bash_profile:
	* type in the terminal:
		[name of your text editor] .bash_profile
			e.g. nano .bash_profile
	* add the following line to your bash_profile:
		export PATH="[absolute path to your blast bin folder]:$PATH"
			e.g. export PATH="/Applications/ncbi-blast-2.6.0+/bin:$PATH"
	* save bash_profile and open a new terminal window
