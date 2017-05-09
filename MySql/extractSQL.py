from Bio import SeqIO
import os
import re
import gzip

#Identifes which genomes to parse
genomeList = ['../E_Coli', '../Agrobacterium']
fileList = []

#Creates the variables used in each table
genome_id, tax_id, short_name, long_name, size_bp, domain, genbank_accession, genbank_date = '', '', '', '', '', '', '', ''
gene_id, replicon_id, locus_tag, gene_name, strand, exon_number, geneSize, product = '', '', '', '', '', '', '', ''
replicon_id, rep_name, num_genes, replicon_type, replicon_structure =  '', '', '', '', ''
exon, left_pos, right_pos, size_bp = '', '', '', ''

#Defines variables used in the functions, synonym and references tables
if False:
	function = ''
	synonym = ''
	external_db, external_id = '',''

# Gets the genbank files for each genome
for genome in genomeList:
	for filename in os.listdir(genome):
		if  re.findall('genomic.gbff.gz$', filename):
			fileList.append(genome+'/'+filename)

#goes through each genbank file to extract data
print fileList
writeFile = open('genomeTables.txt', 'w')

#Opens all files that will  be written to
genomesFile = open('genomesTable.txt','w')
genesFile = open('genesTable.txt', 'w')
repliconsFile = open('repliconsTable.txt','w')
referencesFile = open('referencesTable.txt','w')
synonymsFile = open('synonymsTable.txt','w')
exonsFile = open('exonsTable.txt','w')
functionsFile = open('functionsTable.txt','w')

genome_id = 1
gene_id = 1
num_genes = 0
for filename in fileList:
	#Unzips and gets the file ready to parse
	f = gzip.open(filename, 'rt')
	entryInfo = []
	replicon_id = 1
	for record in SeqIO.parse(f, 'genbank'):
		#record.annotations
		replicon_type = ''
		replicon_structure = ''
		#Finds the shortest name for the genome
		beg = genomeList[genome_id -1].find('../') + 3
		short_name = genomeList[genome_id -1]
		short_name = short_name[beg:len(short_name)]
		rep_name = short_name
		#finds the long name for the genome
		long_name = record.description
		entryInfo.append(long_name)
		#Finds the  size of basepairs
		lenBP = record.annotations['contig']
		beg  = lenBP.find('..')
		size_bp =  lenBP[beg+2: len(lenBP)-1]
		#Finds the domain, accession number and genbank dates
		domain = record.annotations['taxonomy'][0]
		genbank_accession = record.annotations['accessions'][0]
		genbank_date = record.annotations['date']
		#Loops through the entries in the file
		for feat in record.features:
			#Extracts the tax ID
			if feat.type == 'source':
				tax_id = feat.qualifiers.get('db_xref')[0]
				tax_id = tax_id[len('taxon:'):len(tax_id)]
				entryInfo.append(tax_id)
				#adds the replicon structure
				replicon_structure = feat.qualifiers.get('chromosome')
				if replicon_structure:
					replicon_structure = ', '.join(replicon_structure)
				else:
					replicon_structure = 'circular'
				#Extracts the replicon type
				replicon_type = feat.qualifiers.get('plasmid')
				print type(replicon_type)
				if replicon_type:
					replicon_type = 'plasmid'
				else:
					replicon_type = 'chromosome'	
			if feat.type == 'CDS':
				q = feat.qualifiers
				#extracts the strand location
				strand = str(feat.location)
				beg = strand.find('](')+2
				strand = strand[beg:beg+1]
				gene_name = q.get('gene')
				if gene_name:
					gene_name = ', '.join(gene_name)
				synonym = q.get('gene_synonym')
				#Extracts the references
				references = ''
				if q.get('db_xref'):	
					references = q.get('db_xref')
					references =  ''.join(references)
					references = references[3:len(references)]
				#Gets the locus tag
				locus_tag = ', '.join(q.get('locus_tag'))
				exon = locus_tag
				#Gets the gene product, if it exists
				product= ''	
				if q.get('product'):
					product = q.get('product')
					product = ', '.join(q.get('product'))
				#Gets the gene size
				geneSize = str(feat.location)
				beg = geneSize.find('[')
				mid = geneSize.find(':')
				end = geneSize.find(']')
				first = geneSize[beg+1:mid]
				last = geneSize[mid+1:end]
				#Excludes sizes for pseudogenes
				if first[0] != '<' and last[0] != '>':
					geneSize = int(last) - int(first) + 1
				else:
					geneSize = ''
				#sets the default for genes with 1 exon
				exon_split = [feat.location]
				exon = 1
				#finds all exons if they exist
				if 'join' in str(feat.location):
					exon_split = str(feat.location).strip('join{').split(',')
					exon_num = len(exon_split)
					print exon_split
				#gets information for each exon and writes to file
				for exn in exon_split:
					exn = str(exn)
					beg = exn.find('[')
					mid = exn.find(':')
					end = exn.find(']')
					first = exn[beg+1:mid]
					last = exn[mid+1:end]
					#finds exon length for non-pseudogenes
					if first[0] != '<' and last[0] != '>':
						size_bp = int(last) - int(first) + 1
						left_pos = first
						right_pos = last
					#Stores zero location values for pseudogenes
					else:
						size_pb = 0
						left_pos = 0 
						right_pos = 0
					#Writes to  the exons table	
					exonCols = [gene_id, exon, left_pos, right_pos, size_bp]
					exonCols = [str(item) for item in exonCols]
					exonsFile.write('\t'.join(exonCols) + '\n')
					exon += 1
				#Writes to the Gene containing file\
				exon_number = exon-1
				geneCols = [gene_id, genome_id, replicon_id, locus_tag, gene_name, strand, exon_number, geneSize, product]
				geneCols = [str(item) for item  in geneCols]
				genesFile.write('\t'.join(geneCols) + '\n')
				gene_id += 1
				num_genes  += 1
			#Creates the infrastructure for other tables
			if False:
				referenceCols = [gene_id, external_db, external_id]
				referencesFile.write('\t'.join(referenceCols) + '\n')
				
				functionCols = [gene_id, function]
				functionsFile.write('\t'.join(functionCols) + '\n')
				
				synonymCols = [gene_id, synonym]
				synonymsFile.write('\t'.join(synonymCols) + '\n')	
		#Writes to the replicons table
		repliconCols = [replicon_id, genome_id, rep_name, num_genes, replicon_type, replicon_structure]
		repliconCols = [str(item) for item in repliconCols]
		repliconsFile.write('\t'.join(repliconCols) + '\n')
		replicon_id += 1
		num_genes = 0
	#Writes the genomes table
	genomeCols = [genome_id, tax_id, short_name, long_name, size_bp, domain, genbank_accession, genbank_date]
	genomeCols = [str(item) for item in genomeCols]
	genomesFile.write('\t'.join(genomeCols)+'\n')
	genome_id += 1
