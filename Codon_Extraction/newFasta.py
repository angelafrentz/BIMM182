#Reads the genome Sequence
fastaFile = open('genome.fna', 'r')
seq = ''
fastaFile.readline()
for line in fastaFile:
	seq += line.strip()
#Gets individual protein data from the table
tblFile = open('protein_tbl.txt','r')
compDict = {'A': 'T', 'C' : 'G', 'G':'C', 'T':'A'}
tblFile.readline()
writeFile = open('newFasta.txt','w')
for line in tblFile:
	#Writes the header line
	tblInfo = line.strip().split('\t')
	header = '>' + tblInfo[8] +'|' +  tblInfo[6] + '|' +  tblInfo[7] + '\n'
	writeFile.write(header)
	protSeq = list(seq[int(tblInfo[2])-1:int(tblInfo[3])])
	#finds reverse strand sequences
	if tblInfo[4] == '-':
		i = 0
		while i < len(protSeq):
			letter = protSeq[i]
			protSeq[i] = compDict[letter]
			i += 1
		protSeq = reversed(protSeq)
	protSeq = ''.join(protSeq)
	#Prints the sequence	
	i = 0
	while i < len(protSeq):
		if (i + 71) <= len(protSeq):
			writeFile.write(protSeq[i:i+71] + '\n')
		else:
			writeFile.write(protSeq[i:len(protSeq)] + '\n')
		i +=  71
