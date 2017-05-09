#opens and reads the entire genome file
fastaFile = open('genome.fna','r')
seq = ''
fastaFile.readline()
for line in fastaFile:
	seq+= line.strip()
#opens the table file and a file to write to
tblFile = open('protein_tbl.txt', 'r')
compDict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
tblFile.readline()
writeFile = open('tag_seq.txt','w')
for line in tblFile:
	#accesses data from the table and genome files
	tblInfo = line.strip().split('\t')
	ln = tblInfo[7] + '\t'
	protSeq = list(seq[int(tblInfo[2])-1:int(tblInfo[3])])
	#Reverses reverse strand sequences
	if tblInfo[4] == '-':
		i=0
		while i < len(protSeq):
			letter = protSeq[i]
			protSeq[i] = compDict[letter]
			i += 1
		protSeq = reversed(protSeq)
	#Formats data and writes to file
	protSeq = ''.join(protSeq)
	ln += protSeq + '\n'
	writeFile.write(ln)
#Creates a list of all possible codons
codonList = []
nucList = ['A','C','G','T']
for nuc1 in nucList:
	for nuc2 in nucList:
		for nuc3 in nucList:
			codonList.append(nuc1+nuc2+nuc3)
stopCodons = ['TAG', 'TAA', 'TGA']
#reads in the newly created file with an ID and sequence
readFile = open('tag_seq.txt','r')
#creates a counts array
idArr = []
countsArr = []
#error array catches genes with stop codons before the end of the array
#and genes that are not multipules of 3 nucleotides
errorArray = []
print '\n'
print "Warnings:"
for line in readFile:
	lnInfo = line.strip().split('\t')
	oneln = [0]*64
	seq = lnInfo[1]
	i = 0
	while i < len(seq):
		#catches genes that  do not end in a codon of length three
		if (i == len(seq)-1) or (i == len(seq)-2):
			print lnInfo[0] + " does not end in a codon of length three"
			if lnInfo[0] not in errorArray:
				errorArray.append(lnInfo[0])
			break
		codon = seq[i:i+3]
		#catches stop codons that occur before the end of the sequence
		if (codon in stopCodons) and (i != (len(seq)-3)):
			print lnInfo[0] + " has a stop codon before the end."
			if lnInfo[0] not in errorArray:
				errorArray.append(lnInfo[0])
		iID = codonList.index(codon)
		oneln[iID] += 1
		i+= 3
	#adds a new line to the counts array, and adds the gene to the total list of genes
	countsArr.append(oneln)
	idArr.append(lnInfo[0])
print "\n"
print "The following genes should be excluded in further analysis: " + ', '.join(errorArray)
print "\n"

#gets the total number of counts per gene, and for the entire array
countsPerGene = []
for gene in countsArr:
	geneCounts = 0
	for count in gene:
		geneCounts += count
	countsPerGene.append(geneCounts)
totalCodons = sum(countsPerGene)
print "The total number of codons analyzed is", totalCodons, '\n'

#Gets the total counts per codon
codonCounts = []
i = 0
while i < 64:
	codonCount = 0
	for gene in countsArr:
		codonCount += gene[i]
	codonCounts.append(codonCount)
	i += 1
#calculates the frequencies for the entire genome
genomeFreq = []
for count in codonCounts:
	genomeFreq.append(count/float(totalCodons))
#calculates the frquencies  for each gene
freqArr = []
for gene in countsArr:
	geneCount = sum(gene)
	oneGeneFreq = []
	for count in gene:
		oneGeneFreq.append(count/float(geneCount))
	freqArr.append(oneGeneFreq)
#Calculates the CUI for all genes
cuiArr = []
#loops through all genes
i = 0
while i < len(freqArr):
	#stores the result of each codon in a gene multiplied by the corresponding codon freq in thw whole genome
	multArrys = []
	#loops through all codons
	j = 0
	while j < 64:
		#print freqArr[i][j], genomeFreq[j]
		multArrys.append((freqArr[i][j])*(genomeFreq[j]))
		j+= 1
	cuiArr.append(sum(multArrys))
	i+= 1

#Writes a table of CUI's to a file that is sorted using bash commands
print "Codons and their correponding CUI's are printed in a new file:\n"
writeFile = open('cui_freq.tsv','w')
i = 0 
printCodons = []
for gene in cuiArr:
	printCodons.append(idArr[i]+':\t'+str(cuiArr[i]))
	i += 1
writeFile.write('\n'.join(printCodons))
