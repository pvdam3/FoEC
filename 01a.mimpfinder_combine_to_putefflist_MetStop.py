from datetime import datetime
startTime = datetime.now()

import sys, os, re, glob

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC

def MimpFinder(infile, sc_prefix, motiefje, motiefje_rc, datahandler, distance, mimpsequencesfile, infilename):
	datahandler_list = []
	nrofcompletemimps = 0
	complete_mimp_counter = 0
	nrofincompletemimps = 0
	open(mimpsequencesfile,'w').close()
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
		region_to_find_TIR_mate =''
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
				id='sc' +new_id+ '| '+ 'mimp_downstreamregion:' + str(ir_dict[i][1]+1) + '_' + str(end_pos_after_mimp) + ' (strand: +)'+ ' ' + str(ir_dict[i][2])+'^'+str(m.start()+1)+'-'+str(m.end()), 
				description='')
			datahandler_list.append(region_record)
			print '   >' + region_record.id
			i+=1

			if m.start()-400 > 0:
				region_to_find_TIR_mate_coord1 = m.start()-400
				region_to_find_TIR_mate_coord2 = m.start()
				region_to_find_TIR_mate = seq_record.seq[region_to_find_TIR_mate_coord1:region_to_find_TIR_mate_coord2]
			else:
				region_to_find_TIR_mate_coord1 = 0
				region_to_find_TIR_mate_coord2 = m.start()
				region_to_find_TIR_mate = seq_record.seq[region_to_find_TIR_mate_coord1:region_to_find_TIR_mate_coord2]
			match_TIR_mate = re.finditer(motiefje_rc, str(region_to_find_TIR_mate))
			mimp_completeness = False
			mimp_fasta_record = ''
			for TIR_mate in match_TIR_mate:
				mimp_completeness = True
				location_of_mimp_coord1=region_to_find_TIR_mate_coord1+TIR_mate.start()
				location_of_mimp_coord2=m.end()
				mimp_sequence=seq_record.seq[location_of_mimp_coord1:location_of_mimp_coord2]

			if mimp_completeness:
				nrofcompletemimps+=0.5 # 0.5 will also be added in the reverse complement part of the script
				complete_mimp_counter+=1
				mimp_fasta_record = '>'+infilename+'_mimp_'+str(complete_mimp_counter)+' _contig_'+new_id+'('+str(location_of_mimp_coord1)+':'+str(location_of_mimp_coord2)+')\n'+str(mimp_sequence)+'\n'
				with open(mimpsequencesfile,'a') as mimpsequencesfilewriter:
					mimpsequencesfilewriter.write(mimp_fasta_record)
			else:
				nrofincompletemimps+=1


		for m in match_rc:
			if m.start() <= distance: #ensure that no negative value will arise 
				start_pos_before_mimp = 0
			else: 
				start_pos_before_mimp = m.start()-distance
			ir_dict[n] = [m.start()+1, m.end(), seq_record.seq[m.start():m.end()], seq_record.seq[start_pos_before_mimp:m.start()]]
			region_record = SeqRecord(seq=ir_dict[n][3].reverse_complement(),  #make RC of selected area of 2000bp
				id='sc' +new_id+ '| '+ 'mimp_upstreamregion:' + str(start_pos_before_mimp) + '_' + str(ir_dict[n][0]-1) + ' (strand: -, rev. transc.)'+ ' ' + str(ir_dict[n][2].reverse_complement())+'^'+str(m.start()+1)+'-'+str(m.end()),
				description='')
			datahandler_list.append(region_record)
			print '   >' + region_record.id
			n+=1

			if m.end()+400 < len(seq_record.seq):
				region_to_find_TIR_mate = seq_record.seq[m.end():m.end()+400]
			else:
				region_to_find_TIR_mate = seq_record.seq[m.end():len(seq_record.seq)]
			match_TIR_mate = re.finditer(motiefje, str(region_to_find_TIR_mate))
			mimp_completeness = False
			for TIR_mate in match_TIR_mate:
				mimp_completeness = True

			if mimp_completeness:
				nrofcompletemimps+=0.5
			else:
				nrofincompletemimps+=1
	print 'nrofincompletemimps (1 TIR):', nrofincompletemimps, 'nrofcompletemimps (2 TIRs):', int(nrofcompletemimps)
	print '// Motif found: ' + str(len(datahandler_list)) + 'x.'
	SeqIO.write(datahandler_list, datahandler, 'fasta')
	nrofcompletemimps=int(nrofcompletemimps)
	print '\n// Wrote %s %sbp regions downstream of mimp IR motif to %s' % (nrofcompletemimps, distance, datahandler)
	return nrofcompletemimps, nrofincompletemimps


def Translator(infile, datahandler2, distance):
	datahandler_list2 = []
	for seq_record in SeqIO.parse(infile, 'fasta', IUPAC.ambiguous_dna):
		for i in range(3): 				#(i=0 - i=1 - i=2)
			s = seq_record.seq[i:]
			while len(s)%3 != 0: # add N to the end of the region if not devisable by 3
				s += 'N'
			translated_record = SeqRecord(seq=s.translate(table=1), id=seq_record.id, description=seq_record.description+' frame='+str(i))
			datahandler_list2.append(translated_record)
	SeqIO.write(datahandler_list2, datahandler2, 'fasta')
	print '-'*20
	print '// Wrote %s translated %iaa entries to %s' % (len(datahandler_list2), (distance/3), datahandler2)
	print '-'*20
			
def OrfFinder(datahandler2, min_prot_len, datahandler3, orfs, max_prot_len, max_d2m):
	for seq_record in SeqIO.parse(datahandler2, 'fasta', IUPAC.protein):
		genomic_region_up = int(seq_record.description.split('region:')[1].split('_')[0])
		genomic_region_down = int(seq_record.description.split('region:')[1].split('_')[1].split(' ')[0])
		strand = seq_record.description.split('strand:')[1][1]
		frame = int(seq_record.description.split('frame=')[1][0])
		mimpseq = seq_record.description.split(') ')[1].split(' ')[0]
		#print '\n'+genomic_region_up + ' - ' +genomic_region_down + ' ' +strand
		MetStop(seq_record.seq, seq_record.id, genomic_region_up, genomic_region_down, strand, frame, mimpseq, min_prot_len, datahandler3, orfs, max_prot_len, max_d2m)
	SeqIO.write(orfs, datahandler3, 'fasta')
	print '// Wrote %s protein sequences with a mimp IR motif in their promoter and an ORF >%i and <%i aa to %s' % (len(orfs), min_prot_len, max_prot_len, datahandler3)
	print '-'*20

def MetStop(sequence, ident, genomic_region_up, genomic_region_down, strand, frame, mimpseq123, min_prot_len, datahandler3, orfs, max_prot_len, max_d2m):
	met_location = sequence.find('M')
	stop_location = met_location+1
	while met_location >= 0 and stop_location>=0:
		stop_location = sequence.find('*', met_location+1)
		prot = sequence[met_location:stop_location+1] 	#+1 = add STOP
		dist_to_mimp = (met_location*3)+frame		 	#genomic distance
		dist_to_mimp_fromstop = (stop_location*3)+frame+2
		if len(prot) > min_prot_len:
			#print 'M:' + str(met_location) + '; *:' + str(stop_location) + ': '+ prot
			if strand == '+':
				startpos = genomic_region_up+dist_to_mimp
				endpos = genomic_region_up+dist_to_mimp_fromstop
				orf_record = SeqRecord(seq=prot.strip('*'), id=str(ident) +str(startpos)+'-'+str(endpos)+'|d2m:'+str(dist_to_mimp)+'bp|'+strand+'|f'+str(frame)+'|len:'+str(len(prot)-1), description=mimpseq123)
				if len(prot) < int(max_prot_len):
					if int(dist_to_mimp) < max_d2m:
						orfs.append(orf_record)
			elif strand == '-':
				endpos = genomic_region_down-dist_to_mimp
				startpos = genomic_region_down-dist_to_mimp_fromstop
				orf_record = SeqRecord(seq=prot.strip('*'), id=str(ident) +str(startpos)+'-'+str(endpos)+'|d2m:'+str(dist_to_mimp)+'bp|'+strand+'|f'+str(frame)+'|len:'+str(len(prot)-1), description=mimpseq123)
				if len(prot) < int(max_prot_len):
					if int(dist_to_mimp) < max_d2m:
						orfs.append(orf_record)
		met_location = sequence.find('M', met_location+1)
	#return met_location

def OrfWriter(datahandler3, signalpfile, min_prot_len, proteinoutfile, SignalPpath, SignalP_threshold):
	RunSignalP(datahandler3, signalpfile, 'euk', SignalPpath, SignalP_threshold)
	SPorfs = []
	#Parse SignalP output file:
	with open(signalpfile, 'r+') as sp:
		for line in sp:
			if line.startswith('Name='):
				sp_name = line.split('Name=')[1].split('\tSP=')[0]
				orfs = SeqIO.parse(datahandler3, 'fasta', IUPAC.protein)
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
								orfie = SeqRecord(seq=signalpeptideseq.lower()+orf.seq[int(cleavage_site_pos2)-1:].upper(), id=orf.id, description=orf.description.split(' ')[1]+'_signalpeptideseq='+signalpeptideseq)
								SPorfs.append(orfie)
	sp.close()
	print '   Total # of sequences matching the criteria: %i' % len(SPorfs)
	
	#Write all protein sequences that meet requirements (close to motif, longer than 30 aa and contains signal peptide) to proteinoutfile:
	SeqIO.write(SPorfs, proteinoutfile, 'fasta')
 
def RunSignalP(datahandler3, signalpfile, organism, SignalPpath, SignalP_threshold):
	print '// Running SignalP 4.1...'
	cline = SignalPpath+' -t %s -f summary -u %s %s > %s' % (organism, SignalP_threshold, datahandler3, signalpfile)
	os.system(cline)


def ExtractOrfToFasta(proteinsfasta, uberinfile, puteff_dnaseqs, genome, puteff_logfile, combined_puteff_fasta, combined_puteff_logfile, filecounter, nrofcompletemimps, nrofincompletemimps, combined_puteff_logfile2, sc_prefix):
	#EXTRACT GENOMIC SEQUENSE:
	open(puteff_dnaseqs, 'wb').close() 						#clear genomic DNA sequence fastafile
	dnaoutfile = open(puteff_dnaseqs, 'a')
	open(puteff_logfile, 'wb').close()
	
	logheader="genome\tputeff_supercontig\tgenomic_start_pos\tgenomic_end_pos\tdist2mimp\torientation\tprotlength\tD_value\tmimp_IR_seq\tmimp_IR_pos\tsignalpeptideseq\tproteinseq\tgenomicsequence\n"
	logheader2="genome\tnr_of_complete_mimps\tnr_of_incomplete_mimps\tnr_of_put_eff_MetStop\n"
	puteff_logfile_writer = open(puteff_logfile, 'a')
	puteff_logfile_writer.write(logheader) #for each genome, write a log with a header.
	
	combined_putefffile = open(combined_puteff_fasta, 'a')
	combined_puteff_logfile_writer = open(combined_puteff_logfile, 'a')
	combined_puteff_logfile2_writer = open(combined_puteff_logfile2, 'a')
	if filecounter < 1: # During running of the first genome file, the log header should be added at the top in the case of a combined log file.
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
		mimp_IR_pos = seq_record.description.split('^')[1].split('_signalpeptideseq=')[0]
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
				putEff_fastaentry = ">"+str(n).zfill(4)+'.'+signalpeptideseq+"_"+genome+"_d2m"+str(dist2mimp)+"_len"+str(protlength)+"\n"+str(genomicsequence)+"\n\n"
				dnaoutfile.write(putEff_fastaentry)
				
				puteff_attributes = [genome, puteff_supercontig, genomic_start_pos, genomic_end_pos, dist2mimp, orientation, protlength, D_value, mimp_IR_seq, mimp_IR_pos, signalpeptideseq, proteinseq, genomicsequence]
				putEff_logentry = ('\t'.join(map(str,puteff_attributes)))+'\n'
				puteff_logfile_writer.write(putEff_logentry)
				
				#combined_putefffile will collect all the output from the mimpsearch; this means there will be many redundant put effectors.
				combined_putefffile.write(putEff_fastaentry)
				combined_puteff_logfile_writer.write(putEff_logentry)
	putEff_logentry2 = genome+'\t'+str(nrofcompletemimps)+'\t'+str(nrofincompletemimps)+'\t'+str(n)+'\n'
	combined_puteff_logfile2_writer.write(putEff_logentry2)	
	dnaoutfile.close() #collects inside genome out folder all DNA sequences of the putative effectors
	puteff_logfile_writer.close() #writes a log for all puteff found in the current genome (inside genome out folder)
	combined_putefffile.close() #collects inside the out folder all DNA sequences of the putative effectors of all genomes that are being processed by the script.
	combined_puteff_logfile_writer.close() #writes a general log file with more details of the putative effectors identified
	combined_puteff_logfile2_writer.close()
	print '-'*20
	print "// Finished with genome of %s; wrote %i genomic DNA sequences of putEff ORFs to %s" % (genome, n, puteff_dnaseqs)

def MainDef(genomefastafile, directory, folder, combined_puteff_fasta, combined_puteff_logfile, filecounter, combined_puteff_logfile2, combined_puteff_dir):

	 #(default contig prefix)
	sc_prefix = sys.argv[3] #'contig_'
	infilename, infileextension = os.path.splitext(genomefastafile)
	infile      	= directory+folder+'/'+genomefastafile
	outdirectory	= combined_puteff_dir+infilename+'_mimpfinder_out/'
	if not os.path.exists(outdirectory):
		os.makedirs(outdirectory)
	mimpsequencesfile = outdirectory+infilename+'_0_mimpfinder_completemimpsequences.fasta'
	datahandler		= outdirectory+infilename+'_1_mimpfinder_downstreamregions.fasta'
	datahandler2	= outdirectory+infilename+'_2_mimpfinder_translateddownstreamregions.fasta'
	datahandler3	= outdirectory+infilename+'_3_mimpfinder_putativeORFs.fasta'
	signalpfile 	= outdirectory+infilename+'_4_mimpfinder_SignalP.summary_out'
	proteinoutfile 	= outdirectory+infilename+'_5_mimpfinder_proteinseq_out.fasta'
	blastoutfile 	= outdirectory+infilename+'_6_mimpfinder_blast_out.csv'
	puteff_dnaseqs 	= outdirectory+infilename+'_7_mimpfinder_puteff_genomicseq_out.fasta'
	puteff_logfile 	= outdirectory+infilename+'_8_mimpfinder_puteff_logfile.txt'
	motiefje    	= 'TT[TA]TTGC..CCCACTG..'
	motiefje_rc 	= '..CAGTGGG..GCAA[TA]AA'
					   #..CAGT[GA]G[GA]..GCAA[TAG]AA (Mara Bergemann, 2008)
	orfs 			= []


	distance 		= int(sys.argv[4]) # sequence downstream of motif used for ORF prediction
	min_prot_len	= int(sys.argv[5]) #in aa
	max_prot_len	= int(sys.argv[6])
	max_d2m			= int(sys.argv[7]) # max distance between mimp TIR and start-codon
	SignalPpath 	= str(sys.argv[8])
	SignalP_threshold = str(sys.argv[9])

	nrofcompletemimps, nrofincompletemimps = MimpFinder(infile, sc_prefix, motiefje, motiefje_rc, datahandler, distance, mimpsequencesfile, infilename) #ir_dict[i] = [m.start()+1, m.end(), seq_record.seq[m.start():m.end()], seq_record.seq[m.end():m.end()+distance]]
	
	Translator(datahandler, datahandler2, distance)
	OrfFinder(datahandler2, min_prot_len, datahandler3, orfs, max_prot_len, max_d2m)
	OrfWriter(datahandler3, signalpfile, min_prot_len, proteinoutfile, SignalPpath, SignalP_threshold)
	ExtractOrfToFasta(proteinoutfile, infile, puteff_dnaseqs, infilename, puteff_logfile, combined_puteff_fasta, combined_puteff_logfile, filecounter, nrofcompletemimps, nrofincompletemimps, combined_puteff_logfile2, sc_prefix)
	
	return mimpsequencesfile

	
if __name__ == "__main__":
	######################
	directory_folder = sys.argv[1]

	folder = sys.argv[1].split('/')[-1]
	directory = directory_folder.split(folder)[0]
	output_dir = sys.argv[2]
	

	file_extensions = (".fa", ".fasta", ".fas", "fna") # Specify the suffix of the genome files (.fa, .fasta, etc)
	combined_puteff_dir = output_dir+'/01.mimpfinder/'+folder+'_MetStopOut/'
	
	if not os.path.exists(combined_puteff_dir):
		os.makedirs(combined_puteff_dir)
	combined_puteff_fasta	= combined_puteff_dir+'all_putative_effectors.fasta'
	open(combined_puteff_fasta, 'w').close() 
	combined_puteff_logfile = combined_puteff_dir+'all_putative_effectors_log.txt'
	open(combined_puteff_logfile, 'w').close() 
	combined_puteff_logfile2 = combined_puteff_dir+'all_putative_effectors_log2.txt'
	open(combined_puteff_logfile2, 'w').close() 
	######################
	mimpsequencesfile_list = []
	filecounter=0
 	for genomefastafile in os.listdir(directory_folder):
 		if genomefastafile.endswith(file_extensions): 
			print "\n// executing mimpfinder_MetStop script for "+genomefastafile
			mimpsequencesfile = MainDef(genomefastafile, directory, folder, combined_puteff_fasta, combined_puteff_logfile, filecounter, combined_puteff_logfile2, combined_puteff_dir)
			mimpsequencesfile_list.append(mimpsequencesfile)
			filecounter+=1
	 	else:
			print '-'*20
			print '// THE END'
			print "No more files with extension '%s' were found in directory '%s'" % (file_extensions, (directory+folder))
			print "Executed the script for %i files." % (filecounter)
			print 'Total run time = %s' % (datetime.now()-startTime)
			print '-'*20	

	ConcatenateMimpSequences_cmnd = 'cat '+' '.join(mimpsequencesfile_list)+' > '+combined_puteff_dir+'all_mimp_sequences.fasta'
	print ConcatenateMimpSequences_cmnd, os.system(ConcatenateMimpSequences_cmnd)
