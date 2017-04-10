from datetime import datetime
startTime = datetime.now()

import sys, os, re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from BCBio import GFF

def MimpFinder(infile, sc_prefix, motiefje, motiefje_rc, datahandler, distance):
	datahandler_list = []
	print '-'*20
	for seq_record in SeqIO.parse(infile, 'fasta', IUPAC.ambiguous_dna):
		print '// Analyzing \'%s\'... (%i bp long)' % (seq_record.description, len(seq_record))
		if str(sc_prefix) in seq_record.description:
			new_id = seq_record.description.split(sc_prefix)[1].split('_')[0]  #this splitting of '_' is necessary for SignalP in order to not skip all queries in datahandler3 (because their headers might contain '_').
			if ' ' in new_id: #if 'contig_xx' is not the last phrase in the fasta header, resulting in ne_id to become 'contig_xx whole shotgun sequencing bladiebla'
				new_id = new_id.split(' ')[0]
		else:
			print  '-'*20 + '\nThis supercontig prefix was not found; add supercontig prefix as 1st argument..!\n' + '-'*20
			quit()
		match = re.finditer(motiefje, str(seq_record.seq))
		match_rc = re.finditer(motiefje_rc, str(seq_record.seq))
		ir_dict = {}
		i=0
		n=0 
		for m in match:
			if m.end()+distance >= len(seq_record): #ensure that no value bigger than the supercontig length will arise 
				end_pos_after_mimp = len(seq_record)
			else: 
				end_pos_after_mimp = m.end()+distance
			#ir_dict 0=start position, 1=end position, 2=motif sequence^location, 3=ORF-search-sequence
			ir_dict[i] = [m.start()+1, m.end(), seq_record.seq[m.start():m.end()], seq_record.seq[m.end():end_pos_after_mimp]]
			region_record = SeqRecord(seq=ir_dict[i][3], 
				id='contig' +new_id+ '_'+ 'mimp_downstreamregion' + str(ir_dict[i][1]+1) + '-' + str(end_pos_after_mimp) + '_strand1'+ '_' + str(ir_dict[i][2])+'^'+str(m.start()+1)+'-'+str(m.end()), 
				description='')
			datahandler_list.append(region_record)
			print '   >' + region_record.id
			
			###[comment these lines to not search _upstream_ of a forward motif:]###
			# if m.start() <= distance: #ensure that no negative value will arise 
			# 	start_pos_before_mimp = 0
			# else: 
			# 	start_pos_before_mimp = m.start()-distance
			# ir_dict[n] = [m.start()+1, m.end(), seq_record.seq[m.start():m.end()], seq_record.seq[start_pos_before_mimp:m.start()]]
			# region_record = SeqRecord(seq=ir_dict[n][3].reverse_complement(),  #make RC of selected area of 2000bp
			# 	id='contig' +new_id+ '_'+ 'mimp_upstreamregion' + str(start_pos_before_mimp) + '-' + str(ir_dict[n][0]-1) + '_strand-1'+ '_' + str(ir_dict[n][2].reverse_complement())+'^'+str(m.start()+1)+'-'+str(m.end()),
			# 	description='')
			# datahandler_list.append(region_record)
			# print '   >' + region_record.id
			#######################################################################
			i+=1

		for m in match_rc:
			if m.start() <= distance: #ensure that no negative value will arise 
				start_pos_before_mimp = 0
			else: 
				start_pos_before_mimp = m.start()-distance
			ir_dict[n] = [m.start()+1, m.end(), seq_record.seq[m.start():m.end()], seq_record.seq[start_pos_before_mimp:m.start()]]
			region_record = SeqRecord(seq=ir_dict[n][3].reverse_complement(),  #make RC of selected area of 2000bp
				id='contig' +new_id+ '_'+ 'mimp_upstreamregion' + str(start_pos_before_mimp) + '-' + str(ir_dict[n][0]-1) + '_strand-1'+ '_' + str(ir_dict[n][2].reverse_complement())+'^'+str(m.start()+1)+'-'+str(m.end()),
				description='')
			datahandler_list.append(region_record)
			print '   >' + region_record.id

			###[comment these lines to not search _upstream_ of a forward motif:]###
			# if m.end()+distance >= len(seq_record): #ensure that no value bigger than the supercontig length will arise 
			# 	end_pos_after_mimp = len(seq_record)
			# else: 
			# 	end_pos_after_mimp = m.end()+distance
			# #ir_dict 0=start position, 1=end position, 2=motif sequence^location, 3=ORF-search-sequence
			# ir_dict[i] = [m.start()+1, m.end(), seq_record.seq[m.start():m.end()], seq_record.seq[m.end():end_pos_after_mimp]]
			# region_record = SeqRecord(seq=ir_dict[i][3], 
			# 	id='contig' +new_id+ '_'+ 'mimp_downstreamregion' + str(ir_dict[i][1]+1) + '-' + str(end_pos_after_mimp) + '_strand1'+ '_' + str(ir_dict[i][2])+'^'+str(m.start()+1)+'-'+str(m.end()), 
			# 	description='')
			# datahandler_list.append(region_record)
			# print '   >' + region_record.id
			#######################################################################
			n+=1
	print '// Motif found: ' + str(len(datahandler_list)) + 'x.'
	SeqIO.write(datahandler_list, datahandler, 'fasta')
	print '\n// Wrote %s %sbp regions downstream of mimp IR motif to %s' % (len(datahandler_list),distance, datahandler)
	#return ir_dict

def PredictGenes(datahandler, datahandler2, datahandler3a, datahandler3b, datahandler3c, genome, max_d2m, AUGUSTUS_command, min_prot_len, max_prot_len):
	ORF_list = []
	CDS_list = []
	protein_list = []
	
	genecounter = 0
	
	print '// Running Augustus 3.1...'
	print AUGUSTUS_command
	cline = AUGUSTUS_command+' --species=fusarium --singlestrand=true %s > %s' % (datahandler, datahandler2) #--singlestrand=true
	print(cline),os.system(cline)

	in_seq_handle = open(datahandler)
	seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
	in_seq_handle.close()
 
	in_handle = open(datahandler2)

	for rec in GFF.parse(in_handle, base_dict=seq_dict):
		print '>'+rec.id
		feature_ids =[]
		for feature in rec.features:
			if feature.id != '': feature_ids.append(feature.id) # find out how many genes were predicted in this 2.5kb region
		
		for feature in rec.features: 
			#print feature.strand
			for feature_id in feature_ids:
				if feature.id == feature_id: # each gene is iterated in the gff file
					print "   "+"_"*20
					
					print '   gene found:', feature.id, feature.location, feature.strand
					#print feature.extract(rec).seq
					CDS = '' ##
					nrofexons = 0
					exonlocations = []
					for sub_feature in feature.sub_features:
						#print sub_feature.type, sub_feature.location
						if sub_feature.type == 'CDS':
							gene_orientation = sub_feature.strand # gene orientation compared to 2.5kb region > always with IR as 5' end
							
							nrofexons +=1
							if gene_orientation == 1:
								CDS+=sub_feature.extract(rec).seq #add all exons together to form full CDS
							else:
								CDS+=sub_feature.extract(rec).seq.reverse_complement()
							exonlocations.append([int(sub_feature.location.start), int(sub_feature.location.end)])
					#if str(CDS[:3]) == 'ATG':
					#print 'gene orientation:',gene_orientation
					#print 'exonlocations:', exonlocations
					nrofintrons = nrofexons-1 
					
					all_locations = []	
					if nrofexons > 1: #are there any introns?
						for exonlocation in exonlocations:
							all_locations.append(exonlocation[0])
							all_locations.append(exonlocation[1]) #add coordinates of exons to one list all_locations
					intronlocations = []
					for i in range(nrofintrons):
						intronlocations.append([all_locations[(i*2)+1], all_locations[(i+1)*2]]) #extract the coordinates of each intron
					
					ORF = '' ##
					for e in range(len(intronlocations)):
						ORF += rec.seq[exonlocations[e][0]:exonlocations[e][1]].upper()
						ORF += rec.seq[intronlocations[e][0]:intronlocations[e][1]].lower()
					ORF += rec.seq[exonlocations[len(exonlocations)-1][0]:exonlocations[len(exonlocations)-1][1]].upper()
					
					
					contig	= rec.id.split('contig')[1].split('_mimp')[0]
					region_start	= rec.id.split('region')[1].split('_')[0].split('-')[0]
					#print 'region_start:', region_start
					region_end		= rec.id.split('region')[1].split('_')[0].split('-')[1]
					#print 'region_end:', region_end
					strand 			= int(rec.id.split('_strand')[1].split('_')[0])
					#print 'strand', strand
					mimp_IR_seq 	= rec.id.split('_strand')[1].split('_')[1].split('^')[0]
					mimp_IR_loc 	= rec.id.split('_strand')[1].split('_')[1].split('^')[1]
					gene_start_relative		= feature.location.start
					gene_end_relative		= feature.location.end
					if gene_orientation == 1:
						if strand == 1:	
							real_orientation = 1  #gene is found on coding strand downstream of MIMP IR
							gene_start_absolute		= int(region_start)+int(gene_start_relative)
							gene_end_absolute		= int(region_start)+int(gene_end_relative)-1
						elif strand == -1:  
							real_orientation = -1 #gene is found on other strand upstream of MIMP IR
							gene_start_absolute		= int(region_end)-int(gene_start_relative)
							gene_end_absolute		= int(region_end)-int(gene_end_relative)+1
					elif gene_orientation == -1:
						ORF = ORF.reverse_complement()
						CDS = CDS.reverse_complement()
						if strand == 1:	
							real_orientation = -1
							gene_start_absolute		= int(region_start)+int(gene_start_relative)
							gene_end_absolute		= int(region_start)+int(gene_end_relative)-1
						elif strand == -1: 
							real_orientation = 1
							gene_start_absolute		= int(region_end)-int(gene_start_relative)
							gene_end_absolute		= int(region_end)-int(gene_end_relative)+1
					proteinseq = CDS.translate() ##
					proteinseq = proteinseq[:-1] #trim off the STOP codon ('*')
					
					#new_fastaheader = feature.id.split('.t')[0]+'_'+genome+'_contig'+contig+'_'+str(gene_start_absolute)+'-'+str(gene_end_absolute)+'_'+str(real_orientation)+'_d2m'+str(gene_start_relative)+'_'+str(mimp_IR_seq)+'_'+mimp_IR_loc+'_'+str(nrofexons)+'exons'
					#new_fastaheader_forsignalp = feature.id.split('.t')[0]+'_'+genome+'_contig'+contig+'_'+str(gene_start_absolute)+'-'+str(gene_end_absolute)
					
					new_fastaheader_id = 'sc'+contig+'|'+str(gene_start_absolute)+'-'+str(gene_end_absolute)+'|'+'d2m:'+str(gene_start_relative)+'bp|'+str(real_orientation)+'|f?'+'|len:'+str(len(proteinseq))
					new_fastaheader_description = str(mimp_IR_seq)+'^'+mimp_IR_loc+'_'+str(nrofexons)+'_exons'
					print new_fastaheader_id, new_fastaheader_description
					# SignalP only accepts short fastaheaders!!
					
					#print proteinseq[:1]
					if str(proteinseq[:1]) == 'M':
						if len(proteinseq) > min_prot_len:
							if len(proteinseq)<max_prot_len:
								if int(gene_start_relative) < max_d2m:
									ORF_list.append(SeqRecord(seq=ORF, id=new_fastaheader_id, description=new_fastaheader_description))
									CDS_list.append(SeqRecord(seq=CDS, id=new_fastaheader_id, description=new_fastaheader_description))
									protein_list.append(SeqRecord(seq=proteinseq, id=new_fastaheader_id, description=new_fastaheader_description))
									genecounter+=1

	print '\n// Putative genes near a mimp IR found in this genome: ' + str(genecounter)
	print '// Wrote ORF (EXONintronEXON), mRNA coding sequences and protein sequences to datahandler_3a/3b/3c.'
	in_handle.close()
	SeqIO.write(ORF_list, datahandler3a, 'fasta')
	SeqIO.write(CDS_list, datahandler3b, 'fasta')
	SeqIO.write(protein_list, datahandler3c, 'fasta')
	
def RunSignalP(datahandler3c, datahandler4, SignalPpath, SignalP_threshold): 
	print '// Running SignalP 4.1...'
	cline = SignalPpath+' -t euk -f summary -u %s %s > %s' % (SignalP_threshold, datahandler3c, datahandler4)
	os.system(cline)


def ParseSignalP(datahandler3a, datahandler3b, datahandler3c, datahandler4, datahandler5a, datahandler5b, datahandler5c, min_prot_len):
	SP_positives = []
	with open(datahandler4, 'r+') as sp:
		for line in sp:
			if line.startswith('Name='):
				sp_name = line.split('Name=')[1].split('\tSP=')[0]
				orfs = SeqIO.parse(datahandler3c, 'fasta')
				for orf in orfs:
					if orf.id == sp_name:
						sp_present = line.split('\tSP=\'')[1].split('\'')[0] # SP: YES OR NO
						orf.id += '|SP=' + sp_present
						orf.id += ';D='+line.split(' D=0')[1].split(' ')[0]  # SP D-value
						if sp_present == 'YES':
							cleavage_site_pos1 = line.split('Cleavage site between pos. ')[1].split(' and ')[0]
							cleavage_site_pos2 = line.split('Cleavage site between pos. ')[1].split(' and ')[1].split(': ')[0]
							if 'X' not in (orf.seq[:int(cleavage_site_pos2)-1]): #prevent 'NNN' translated to 'X' in signalpeptides to end up in the list
								signalpeptideseq = str(orf.seq[:int(cleavage_site_pos2)-1])
								#contig = orf.id.split('_contig')[1].split('_')[0]
								orfie = SeqRecord(seq=signalpeptideseq.lower()+orf.seq[int(cleavage_site_pos2)-1:].upper(), id=orf.id, description=orf.description.split(' ')[1]+'_signalpeptideseq='+signalpeptideseq)
								SP_positives.append(orfie)
								print '  >'+orfie.id+' contains a signal peptide.'
	sp.close()
	print '   Total # of sequences matching the criteria: %i' % len(SP_positives)
	print ''
	
	#Write all protein sequences that meet requirements (close to motif, longer than 30 aa and contains signal peptide) to datahandler5c:
	SeqIO.write(SP_positives, datahandler5c, 'fasta')

	#Write all Open Reading Frame sequences (DNA) to datahandler5a
	datahandler3a_handle = SeqIO.parse(datahandler3a, "fasta")
	SP_positives_ORFs=[]
	for rec in datahandler3a_handle:
		#print '    ???? ', rec.id
		for orfie in SP_positives:
			if rec.id == orfie.id.split('|SP=')[0]:
				#print '    !!!! ', orfie.id.split('|SP=')[0]
				new_rec = SeqRecord(seq=rec.seq, id=orfie.id, description=orfie.description)
				SP_positives_ORFs.append(new_rec)
	SeqIO.write(SP_positives_ORFs, datahandler5a, 'fasta')

	#Write all Coding Sequences (mRNA) to datahandler5b
	datahandler3b_handle = SeqIO.parse(datahandler3b, "fasta")
	SP_positives_CDS=[]
	for rec in datahandler3b_handle:
		#print '    ???? ', rec.id
		for orfie in SP_positives:
			if rec.id == orfie.id.split('|SP=')[0]:
				#print '    !!!! ', orfie.id.split('|SP=')[0]
				new_rec = SeqRecord(seq=rec.seq, id=orfie.id, description=orfie.description)
				SP_positives_CDS.append(new_rec)
	SeqIO.write(SP_positives_CDS, datahandler5b, 'fasta')

				
	

def ExtractOrfToFasta(proteinsfasta, uberinfile, puteff_dnaseqs, genome, puteff_logfile, combined_puteff_fasta, combined_puteff_logfile, filecounter, combined_puteff_logfile2, sc_prefix):
	#EXTRACT GENOMIC SEQUENSE:
	open(puteff_dnaseqs, 'wb').close() 						#clear genomic DNA sequence fastafile
	dnaoutfile = open(puteff_dnaseqs, 'a')
	open(puteff_logfile, 'wb').close()
	
	logheader="genome\tputeff_supercontig\tgenomic_start_pos\tgenomic_end_pos\tdist2mimp\torientation\tprotlength\tnrofexons\tD_value\tmimp_IR_seq\tmimp_IR_pos\tsignalpeptideseq\tproteinseq\tgenomicsequence\n"
	logheader2="genome\tnr_of_put_eff_Augustus\n"
	puteff_logfile_writer = open(puteff_logfile, 'a')
	puteff_logfile_writer.write(logheader)
	
	combined_putefffile = open(combined_puteff_fasta, 'a')
	combined_puteff_logfile_writer = open(combined_puteff_logfile, 'a')
	combined_puteff_logfile2_writer = open(combined_puteff_logfile2, 'a')
	if filecounter < 1: # During running of the first genome file, the log header should be added at the top.
		combined_puteff_logfile_writer.write(logheader)
		combined_puteff_logfile2_writer.write(logheader2)
	proteinsfastafile = SeqIO.parse(proteinsfasta, 'fasta')
	n=0
	print '// Putative effectors found in genome: ' + genome
	for seq_record in proteinsfastafile:
		n+=1
		puteff_supercontig = seq_record.description.split('sc')[1].split('|')[0]
		#MAKE SURE THE CONTIG IN THE ORIGINAL FASTA FILE IS LIKE THIS; contig_1 (space no comma behind it)
		genomic_start_pos = int(seq_record.description.split('|')[1].split('-')[0]) 
		genomic_end_pos = int(seq_record.description.split('|')[1].split('-')[1].split('|')[0])
		dist2mimp = seq_record.description.split('d2m:')[1].split('bp')[0]
		orientation = seq_record.description.split('bp|')[1].split('|')[0]  
		protlength = seq_record.description.split('|len:')[1].split('|')[0]
		D_value = seq_record.description.split('D=.')[1].split(' ')[0]
		mimp_IR_seq = seq_record.description.split(' ')[1].split('^')[0]
		mimp_IR_pos = seq_record.description.split('^')[1].split('_')[0]
		nrofexons = seq_record.description.split('^')[1].split('_')[1].split('_exons')[0]
		signalpeptideseq = seq_record.description.split('_signalpeptideseq=')[1]
		proteinseq = seq_record.seq

		all_contigs = list(SeqIO.parse(uberinfile, "fasta", IUPAC.unambiguous_dna))
 		for sc in all_contigs:
			sc_id = sc.description.split(sc_prefix)[1]
			if "_" in sc_id:
 				sc_id = sc_id.split("_")[0]
			if sc_id == puteff_supercontig:
				
				genomicsequence = sc.seq[genomic_start_pos-1:genomic_end_pos]
				
				print '   contig_'+str(puteff_supercontig)+'\tposition '+str(genomic_start_pos)+'-'+str(genomic_end_pos)+'\t'+signalpeptideseq
				if genomic_start_pos > genomic_end_pos:
					genomicsequence = sc.seq[genomic_end_pos-1:genomic_start_pos]
					
				
				
				putEff_fastaentry = ">"+str(n).zfill(3)+'.'+signalpeptideseq+"_"+genome+"_d2m"+str(dist2mimp)+"_len"+str(protlength)+"\n"+str(genomicsequence)+"\n\n"
				dnaoutfile.write(putEff_fastaentry)
				
				puteff_attributes = [genome, puteff_supercontig, genomic_start_pos, genomic_end_pos, dist2mimp, orientation, protlength, nrofexons, D_value, mimp_IR_seq, mimp_IR_pos, signalpeptideseq, proteinseq, genomicsequence]
				putEff_logentry = ('\t'.join(map(str,puteff_attributes)))+'\n'
				puteff_logfile_writer.write(putEff_logentry)
				
				#combined_putefffile will collect all the output from the mimpsearch; this means there will be many redundant put effectors.
				combined_putefffile.write(putEff_fastaentry)
				combined_puteff_logfile_writer.write(putEff_logentry)
	putEff_logentry2 = genome+'\t'+str(n)+'\n'
	combined_puteff_logfile2_writer.write(putEff_logentry2)			
	dnaoutfile.close() #collects inside genome out folder all DNA sequences of the putative effectors
	puteff_logfile_writer.close() #writes a log for all puteff found in the current genome (inside genome out folder)
	combined_putefffile.close() #collects inside the out folder all DNA sequences of the putative effectors of all genomes that are being processed by the script.
	combined_puteff_logfile_writer.close() #writes a general log file with more details of the putative effectors identified
	print '-'*20
	print "// Finished with genome of %s; wrote %i genomic DNA sequences of putEff ORFs to %s" % (genome, n, puteff_dnaseqs)
	

			


def MainDef(genome, genomefastafile, directory, folder, combined_puteff_fasta, combined_puteff_logfile, filecounter, combined_puteff_logfile2, combined_puteff_dir):
	# Add a string as an argument when running the Python script so to tell what is the supercontig prefix.
	# ex: "python mimpfinder.py Supercontig_"

	sc_prefix = sys.argv[3] #'contig_'
	AUGUSTUS_command = sys.argv[4]+' '+sys.argv[5]
	infilename, infileextension = os.path.splitext(genomefastafile)
	infile      	= directory+folder+'/'+genomefastafile
	outdirectory	= combined_puteff_dir+infilename+'_mimpfinder_out/'

	if not os.path.exists(outdirectory):
		os.makedirs(outdirectory)
	datahandler		= outdirectory+infilename+'_1_mimpfinder_downstreamregions.fasta'
	datahandler2	= outdirectory+infilename+'_2_mimpfinder_predictedgenes.gff'
	datahandler3a	= outdirectory+infilename+'_3a_mimpfinder_ORFs.fasta'
	datahandler3b	= outdirectory+infilename+'_3b_mimpfinder_CDS.fasta'
	datahandler3c	= outdirectory+infilename+'_3c_mimpfinder_proteins.fasta'
	datahandler4 	= outdirectory+infilename+'_4_mimpfinder_SignalP.summary_out'
	datahandler5a 	= outdirectory+infilename+'_5a_mimpfinder_ORFs_SP.fasta'
	datahandler5b 	= outdirectory+infilename+'_5b_mimpfinder_CDS_SP.fasta'
	datahandler5c 	= outdirectory+infilename+'_5c_mimpfinder_proteins_SP.fasta'
	datahandler6 	= outdirectory+infilename+'_6_mimpfinder_ORFs_SP.fasta'
	puteff_logfile 	= outdirectory+infilename+'_7_mimpfinder_puteff_logfile.txt'
	motiefje    	= 'TT[TA]TTGC..CCCACTG..'
	motiefje_rc 	= '..CAGTGGG..GCAA[TA]AA'
					   #..CAGT[GA]G[GA]..GCAA[TAG]AA (Mara Bergemann, 2008)
	orfs 			= []

	distance 		= int(sys.argv[6]) # sequence downstream of motif used for ORF prediction
	min_prot_len	= int(sys.argv[7]) #in aa
	max_prot_len	= int(sys.argv[8])
	max_d2m			= int(sys.argv[9]) # max distance between mimp TIR and start-codon
	SignalPpath		=str(sys.argv[10])
	SignalP_threshold = str(sys.argv[11])

	MimpFinder(infile, sc_prefix, motiefje, motiefje_rc, datahandler, distance) #ir_dict[i] = [m.start()+1, m.end(), seq_record.seq[m.start():m.end()], seq_record.seq[m.end():m.end()+distance]]
	PredictGenes(datahandler, datahandler2, datahandler3a, datahandler3b, datahandler3c, genome, max_d2m, AUGUSTUS_command, min_prot_len, max_prot_len)
	RunSignalP(datahandler3c, datahandler4, SignalPpath, SignalP_threshold)
	ParseSignalP(datahandler3a, datahandler3b, datahandler3c, datahandler4, datahandler5a, datahandler5b, datahandler5c, min_prot_len)
	ExtractOrfToFasta(datahandler5c, infile, datahandler6, infilename, puteff_logfile, combined_puteff_fasta, combined_puteff_logfile, filecounter, combined_puteff_logfile2, sc_prefix)
	


	
if __name__ == "__main__":

	######################
	directory_folder = sys.argv[1]
	folder = sys.argv[1].split('/')[-1]
	directory = directory_folder.split(folder)[0]
	output_dir = sys.argv[2]
	file_extensions = (".fa", ".fasta", ".fas", "fna") # Specify the suffix of the genome files (.fa, .fasta, etc)
	combined_puteff_dir = output_dir+'/01.mimpfinder/'+folder+'_AugustusOut/'


	if not os.path.exists(combined_puteff_dir):
		os.makedirs(combined_puteff_dir)
	combined_puteff_fasta	= combined_puteff_dir+'all_putative_effectors.fasta'
	open(combined_puteff_fasta, 'w').close() 
	combined_puteff_logfile = combined_puteff_dir+'all_putative_effectors_log.txt'
	open(combined_puteff_logfile, 'w').close()
	combined_puteff_logfile2 = combined_puteff_dir+'all_putative_effectors_log2.txt'
	open(combined_puteff_logfile2, 'w').close() 
	######################
	
	filecounter=0
 	for genomefastafile in os.listdir(directory_folder):
 		if genomefastafile.endswith(file_extensions): 
 			genome = genomefastafile.split('.fa')[0]
			print "\n// executing mimpfinder_AUGUSTUS script for "+genomefastafile
			MainDef(genome, genomefastafile, directory, folder, combined_puteff_fasta, combined_puteff_logfile, filecounter, combined_puteff_logfile2, combined_puteff_dir)
			filecounter+=1
	 
 	else:
		print '-'*20
		print '// THE END'
		print "No more files with extension '%s' were found in directory '%s'" % (file_extensions, (directory+folder))
		print "Executed the script for %i files." % (filecounter)
		print 'Total run time = %s' % (datetime.now()-startTime)
		print '-'*20	


