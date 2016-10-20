#this script takes in the output of mimp-identification script "all_putative_effectors.fasta"
#It will take each put effector in this fastafile and blast it against itself. Then it takes these groups 
#together and clusters them with single_linkage. Each of the clusters is assessed for the longest
#sequence, and these longest sequences are written to the outfile (_clustered.fasta)


from Bio import SeqIO
import os, sys
debug=0

def cluster_homologous_effectors(infile, E_VALUE_THRESH, PERC_IDENTITY_THRESH, LENGTH_THRESH, blastdatabasedir, BLASTbindir, infiledir):

	database_store = blastdatabasedir+'/'+infile.split('/')[-1].split('.fa')[0]
	cmnd = BLASTbindir+'/makeblastdb -dbtype nucl -in '+infile+' -out '+database_store
	print "---BLASTDB---\n", cmnd, os.system(cmnd), "---\n" #uncomment if you want to rebuild a blastdb #untag # if you want to build the BLAST database
	
	blastoutfilename = infiledir+infile.split('/')[-1].split('.fa')[0]+'.vs.'+infile.split('/')[-1].split('.fa')[0]+'.blastout'
	cmnd = BLASTbindir+"/blastn -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -query "+infile+' -db '+database_store+' -out '+blastoutfilename+' -evalue '+str(E_VALUE_THRESH)
	#									0	  1	  2	  3	  4		5	      6	 	 7		8	  9	10	 11	   12
	effector2homologs = {}
	
	if os.system(cmnd) == 0:
		
		lines = open(blastoutfilename).readlines()
		#'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen'
		for line in lines:
			tabs = line.strip().split('\t')
			if float(tabs[2]) > PERC_IDENTITY_THRESH and (int(tabs[7]) - int(tabs[6]))/float(tabs[12]) > LENGTH_THRESH: #with last check: from 358 to 354 put effectors with same setup
				if effector2homologs.has_key(tabs[0]): 
					effector2homologs[tabs[0]].add(tabs[1]) #add all BLAST-associated hits to this entry
				else: effector2homologs[tabs[0]] = set([tabs[1]])
				
				if effector2homologs.has_key(tabs[1]): 
					effector2homologs[tabs[1]].add(tabs[0]) #the other way around; also add the BLAST query as an association to each BLAST hit.
				else: effector2homologs[tabs[1]] = set([tabs[0]])
				
	return single_linkage(effector2homologs)

def single_linkage(node_partners):
	
	clusters=[]
	
	nodes=node_partners.keys()
	while len(nodes)>0:
		if debug:
			print 'new cluster, seed node:', nodes[0]
			print '# keys in hash:', len(nodes)
		cluster=update_cluster(node_partners, [nodes[0]], set([]))
		clusters.append(cluster)
		nodes=node_partners.keys()
		if debug:
			print 'finished cluster:', cluster
			print '# keys in hash:', len(nodes)
		
	return clusters

def update_cluster(node_partners, todo, cluster=set([])):
	if debug: print todo, cluster 
	new_todo=set([])
	for t in todo:
		partners=node_partners[t]
		new_todo=new_todo.union(set(partners))
		cluster.add(t)
		del node_partners[t]

	new_todo=new_todo.difference(cluster)
	if len(new_todo)>0: return update_cluster(node_partners, new_todo, cluster)
	else: return cluster				
						
if __name__ == "__main__":

	### arguments being passed from pipeline script: ###
	infile = sys.argv[1]
	infile_filename = infile.split('/')[-1]
	infiledir = infile.split(infile_filename)[0]
	outfile = infile.split('.fa')[0]+'_clustered.fasta'
	blastdatabasedir = sys.argv[2]
	leave_put_eff_identifiers_during_clustering = sys.argv[3]
	BLASTbindir = sys.argv[4]
	###

	outfilewriter = open(outfile, 'w')
	E_VALUE_THRESH = 0.001
	PERC_IDENTITY_THRESH = 60
	LENGTH_THRESH = 0.3

	
	print '\n'
	print "// Running clustering script on file %s" % infile
	clusters = cluster_homologous_effectors(infile, E_VALUE_THRESH, PERC_IDENTITY_THRESH,LENGTH_THRESH, blastdatabasedir, BLASTbindir, infiledir)
	longest_elements = []
	handle = open(infile, "rU")
	record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
	handle.close()
	for c in clusters:
		print '// # of nodes in this cluster = ', str(len(c))

		longest_elem = [0, '']
		
		for elem in c:
			print elem
			#elem_len = int(elem.split('len')[1])
			elem_len = len(record_dict[elem].seq)

			######################################################################
			### in case these pieces of information are going to be necessary: ###
			#elem_d2m = elem.split('_len')[0].split('d2m')[1]
			#elem_isolate = elem.split('_', 1)[1].split('_d2m')[0]
			#elem_identifier = elem.split('.', 1)[0]
			#elem_idmethod = ''
			#if len(str(elem_identifier)) == 3: 
				#elem_idmethod = 'augustus'
			#elif len(str(elem_identifier)) == 4: 
				#elem_idmethod = 'metstop'
			#print elem_d2m, elem_isolate, elem_identifier, elem_idmethod
			######################################################################
			
			if elem_len > longest_elem[0]:
				longest_elem[0] = elem_len
				longest_elem[1] = elem
			print '\t', elem_len, '\t', elem
		print '\t-------'
		print '\tLongest element in cluster is: ', longest_elem[1], '\n'
		longest_elements.append(record_dict[longest_elem[1]])

	#SeqIO.write(longest_elements, outfile, "fasta")
	for longest_element in longest_elements:
		if leave_put_eff_identifiers_during_clustering == 'TRUE':
			longest_element.id = longest_element.id #in case you're clustering a clustered file (that does not contain any description, but only a .id)
		else:
			longest_element.id = longest_element.description.split('.', 1)[1].split('_', 1)[0] #split 1x after identifier (001.) and parse the rest to be the new name

		longest_element.description = ''
		if "\x00" in longest_element.seq:
			print longest_element.id
		
	SeqIO.write(longest_elements, outfile, "fasta")
	outfilewriter.close()
	print '-'*30
	print '// Found %i putative effectors in these genomes: ' % len(clusters)
	for longest_elem in longest_elements:
		print '\t'+longest_elem.id
	print '-'*30
	
	
	
	