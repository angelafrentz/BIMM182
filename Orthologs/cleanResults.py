fileList = ['blast_AvsE.out', 'blast_EvsA.out']
for filename in fileList:
	read = open(filename, 'r')
	out = 'clean_' + filename
	outfile = open(out, 'w')
	for line in read:
		ln = line.strip().split()
		original = ln[0]
		end = original.find('.')
		ln[0] = original[0:end]
		original = ln[1]
		beg = original.find('|') +1
		end = original.find('.')
		ln[1] = original[beg:end]
		outfile.write('\t'.join(ln) + '\n')
