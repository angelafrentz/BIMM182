from Bio import SeqIO
import gzip

f =  gzip.open("GCF_000005845.2_ASM584v2_genomic.gbff.gz", 'rt')
genome = SeqIO.parse(f, 'genbank')
featList = {'db_xref for taxon','protein_id','gene','locus_tag','gene_synonym', 'product',
'EC_number', 'db_xref' }
for record in genome:
	print record.id, len(record.seq)
	feature = record.features
	for feat in feature:
		ln = ''
		sourceInfo = ''
		if feat.type == 'source':
			sourceInfo+= '; TaxID: '+ ''.join(feat.qualifiers.get('db_xref'))
		if feat.type == "CDS":
			q = feat.qualifiers
			#Checks for pseudogenes
			ln += 'ProtID: '
			if q.get('protein_id'):
				ln += ''.join(q.get('protein_id'))
			else:
				ln += 'Pseudogene'
			#Appends other fields to the println
			ln += '; Location/Strand: '+ str(feat.location) 
			ln += '; GeneName: '+ ''.join(q.get('gene'))
			ln += '; GeneTag: '+ ''.join(q.get('locus_tag'))
			ln += '; GeneSynonym: ' + ''.join(q.get('gene_synonym'))
			ln += '; ExternalReferences: ' + ''.join(q.get('db_xref'))
			#Checks for pseudogenes
			ln += '; ProteinName: '
			if q.get('product'):
				ln += ''.join(q.get('product'))
			else:
				ln += 'Pseudogene'
			#Checks For EC Number
			ln += '; ECNumber: '
			if q.get('EC_number'):
				ln += ''.join(q.get('EC_number'))
			else:
				ln += 'None'
			#Gets Externam References
			ln += '; ExternalReferences: ' + ''.join(q.get('db_xref'))
			#prints all information
			print ln + sourceInfo
