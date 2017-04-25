from Bio import SeqIO
import gzip
f =  gzip.open("GCF_000005845.2_ASM584v2_protein.faa.gz", 'rt')
for record in SeqIO.parse(f, 'fasta'):
	print record.id
	print record.seq

