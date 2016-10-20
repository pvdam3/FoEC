from Bio.Blast.Applications import NcbiblastnCommandline, NcbitblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO

import csv, os, re, sys, glob, datetime


def BuildBlastDB(genome_folder, genome_fastafile, blastdatabasedir, BLASTbindir):
	database_in = genome_folder+'/'+genome_fastafile
	genome = genome_fastafile.split('.fa')[0]
	database_store = blastdatabasedir + '/' + genome
	cmnd = BLASTbindir+'/makeblastdb -dbtype nucl -in '+database_in+' -out '+database_store
 	if buildblastdb == 'yes':
		print "// ---BUILDING BLAST DB FROM FASTA FILE---"
		print cmnd, os.system(cmnd)	
	elif buildblastdb == 'no':
		print "// No BLAST db will be built\n\n"
		pass
	elif buildblastdb != 'yes':
		print "// No BLAST db will be built, use set buildblastdb to 'yes' to build a BLAST database from fastafiles."

	return database_store, genome

def BlastAndParse(blastdb, genome, queryfile, E_VALUE_THRESH, PERC_IDENTITY_THRESH, outputdir, BLAST_task):
	#blastoutput = outputdir+BLAST_type+"_presence_absence.blastout"
	blastoutput = outputdir+BLAST_type+queryfile.split('/')[-1]+"_vs"+genome+"_presence_absence.blastout"
	open(blastoutput, 'wb').close() #clearfile
	if BLAST_type == 'BLASTN':
		blastn_cline = NcbiblastnCommandline(query=queryfile, db=blastdb, evalue=E_VALUE_THRESH, outfmt=5, out=blastoutput, dust='no', task=BLAST_task)
		print blastn_cline
		#stdout, stderr = blastn_cline()
		blastn_cline()
	else:
		print "ERROR: NO BLAST TYPE HAS BEEN SPECIFIED (BLASTN OR TBLASTN)"
	
	result_handle = open(blastoutput)
	blast_records = NCBIXML.parse(result_handle) 

	effector2scoresForThisGenome = {}
	
	for blast_record in blast_records:
		for alignment in blast_record.alignments:
			for hsp in alignment.hsps:
				return_perc_id = int(100*hsp.identities/blast_record.query_length)
				#alignment_length = abs(int(hsp.query_end-hsp.query_start))
				#return_perc_id = int(100*hsp.identities/alignment_length)
				########
				return_contig_nr = "??"
				subject_name = alignment.title
				if "ontig_" in subject_name:
					return_contig_nr = subject_name.split("ontig_")[1]
					if "_" in return_contig_nr:
						return_contig_nr = return_contig_nr.split("_")[0]
				if not effector2scoresForThisGenome.has_key(blast_record.query): effector2scoresForThisGenome[blast_record.query]=[]
				if (hsp.expect < E_VALUE_THRESH) and (return_perc_id >= PERC_IDENTITY_THRESH):
					if output_contigs in no: 				#'', 'no' or 'n'
						effector2scoresForThisGenome[blast_record.query].append(return_perc_id)
						#print blast_record.query, '\t', effector2scoresForThisGenome[blast_record.query][0]
					elif output_contigs in yes:
						effector2scoresForThisGenome[blast_record.query].append(return_contig_nr+'('+str(hsp.sbjct_start)+'-'+str(hsp.sbjct_end)+')')
				#else: #requirements are not met, so key with value '[]' should be popped out: """this doesn't work and shows no output whatsoever 26-3-2015"""
					#effector2scoresForThisGenome.pop(blast_record.query, None)
	result_handle.close()
	return effector2scoresForThisGenome
	
	
if __name__ == "__main__":

	### arguments being passed from pipeline script: ###
	queryfile = sys.argv[1]
	queryfilename = queryfile.split('/')[-1]
	queryfolder = queryfile.split(queryfilename)[0]
	
	genome_folder = sys.argv[2]
	genome_folder_b = genome_folder.split('/')[-1]
	genome_folder_a = genome_folder.split(genome_folder_b)[0]
	
	blastdatabasedir = sys.argv[3]
	BLASTbindir = sys.argv[4]
	outputdir = sys.argv[5]+'/03.blastn_presence_absence/'
	if not os.path.exists(outputdir):
		os.makedirs(outputdir)
	###

	date_time_now 			= datetime.datetime.now().strftime("%d-%m-%Y %H:%M")
	queryfile 				= queryfolder+queryfilename
	outfilename 			= outputdir+genome_folder_b+'.vs.'+queryfilename.split('.fa')[0]+'.txt'
	E_VALUE_THRESH 			= 0.001
	PERC_IDENTITY_THRESH 	= int(sys.argv[6])#30 
	BLAST_type 				= 'BLASTN'
	BLAST_task				= sys.argv[7]#'blastn' or megablast
	buildblastdb 			= sys.argv[8]#'yes'
	transpose 				= True
	plot_copynumber 		= False #define whether presence is '1' or copy nr (# of blast hits)
	
	if buildblastdb == 'yes':
		print "// Will start with building BlastDBs of each genome encountered.\n"
	
	yes = set(['yes','y'])
	no = set(['no','n', ''])
	output_contigs = 'n'
	#hierarchicalclusteringquestion = raw_input("Would you like to get hierarchical clustering output in stead of %id? (y/n) \n>")
	hierarchicalclusteringquestion = 'y'
	if hierarchicalclusteringquestion in yes:
		hierarchicalclustering = True
		outfilename = outfilename.split('.txt')[0]+'_hierarch-clust.txt'
		
	else:
		hierarchicalclustering = False
		output_contigs = raw_input("Would you like to get contig numbers in stead of %id? (y/n) \n>")
		if output_contigs in yes:
			outfilename = outfilename.split('.txt')[0]+'_contignrs.txt'
			
	if transpose == True:
		outfilename = outfilename.split('.txt')[0]+'_transposed.txt'

  	query_list = []
   	for record in SeqIO.parse(queryfile, 'fasta'):
   		query_list.append(record.id)
   	
   	###############################################################
   	#cancels all above outputfile naming statements:
   	outfilename = outputdir+'blastn_presence_absence.txt'
  	###############################################################	
	effector2genome2scores = {}
	genomes = []
 	for genome_fastafile in os.listdir(genome_folder):
 		if genome_fastafile.endswith((".fa", ".fasta", ".fas")):
 			blastdb, genome = BuildBlastDB(genome_folder, genome_fastafile, blastdatabasedir, BLASTbindir)
 			genomes.append(genome)
 			#return a dictionary containing all query fastaheaders (keys) with BLAST% as values:
 			effector2scoresForThisGenome = BlastAndParse(blastdb, genome, queryfile, E_VALUE_THRESH, PERC_IDENTITY_THRESH, outputdir, BLAST_task)
 			
 			#order these back to the way they were in the original query fastafile:
 			n=1
 			for effector in query_list:
 				new_record_id = '{0:04}'.format(n)+'..'+effector 			#'..' to reduce mis-splitting
 				if not effector2genome2scores.has_key(new_record_id): effector2genome2scores[new_record_id] = {}
 				#in dictionary effector2scoresForThisGenome, add empty values (for effectors that were not found with BLAST:
 				if not effector2scoresForThisGenome.has_key(effector): 
 					effector2scoresForThisGenome[effector] = []
 				effector2genome2scores[new_record_id][genome] = effector2scoresForThisGenome[effector]
 				n+=1
 		else:
 			print '-'*20
 			print "// No more files with extension .fa, .fas, .fasta were found in directory '%s'" % (genome_folder)
			print '-'*20
			
	# Print to table:

	open(outfilename, 'w').close() 									#clear file
	outfile     = open(outfilename, 'a')	
	
	genomes.sort()
	header = ''
	if hierarchicalclustering == False:
		if output_contigs in no:
			header = BLAST_type+' PERCENT IDENTITY\nDate: '+date_time_now+'\nE-value: '+str(E_VALUE_THRESH)+'\n'+'%ID_threshold: '+str(PERC_IDENTITY_THRESH)+'%\n''query:'
		elif output_contigs in yes:
			header = BLAST_type+' CONTIG NR + SUBJECT LOCATION\nDate: '+date_time_now+'\nE-value: '+str(E_VALUE_THRESH)+'\n'+'%ID_threshold: '+str(PERC_IDENTITY_THRESH)+'%\n''query:'

	if transpose == True:
		for effector in sorted(effector2genome2scores.keys()):
			header += '\t'+effector
	else:
		for g in genomes:
			header += '\t'+g
	header+='\n'
	print header
	outfile.write(header)
	
	if transpose == True:
		for g in genomes:
			out = g
			for effector in sorted(effector2genome2scores.keys()):
				genome2scores = effector2genome2scores[effector]
				scores = genome2scores[g]
				if hierarchicalclustering == True:
					scorestr = 0
					for s in scores:
						if plot_copynumber == True:
							scorestr	+=1
						else:
							scorestr	=1
				else:
					scorestr = ''
					for s in scores:
						scorestr+=str(s)+';'
					scorestr = scorestr[:-1] # laatste ; eraf	
				out += '\t'+str(scorestr)
			print out
			outfile.write(out+'\n')
	else:
		for effector in sorted(effector2genome2scores.keys()):
			out = effector
			genome2scores = effector2genome2scores[effector]
			for g in genomes:
				#print g
				scores = genome2scores[g]
				if hierarchicalclustering == True:
					scorestr = 0
					for s in scores:
						scorestr=1
				else:
					scorestr = ''
					for s in scores:
						scorestr+=str(s)+';'
					scorestr = scorestr[:-1] # laatste ; eraf	
				out += '\t'+str(scorestr)
			print out
			outfile.write(out+'\n')
	outfile.close()
	print '-'*30
	print "// Written data to file: ", outfilename
	print '-'*30
