import MySQLdb as sql
import os

db = sql.connect( host = 'bm185s-mysql.ucsd.edu', user = 'afrentze', passwd = '', db = 'afrentze_db')

#Extracts from the operons table
curs_operons = db.cursor()
curs_operons.execute('Select * from Operons')
posCtrlDist = []
negCtrlDist = []
geneID = {}

#goes through every operon
for line in curs_operons:
	operon = line[0]
	#print operon
	operonGenes =  line[5].strip().split(',')
	bOperonGenes = []
	#Gets the bNum for every operon in the gene
	for gene in operonGenes:
		if '<' not in gene:
			#print gene
			curs_bNum = db.cursor()
			curs_bNum.execute('Select * from genes where name = \'' + gene + '\' limit 1;')
			ln = curs_bNum.fetchall()
			#only acesses genes with a valid bNumber
			if ln:	
				bNum = ln[0][3]
				bOperonGenes.append(bNum)
				#Creates dictionary of operons with the associated bNumbers
				geneID.setdefault(operon, []).append(bNum)
	#Gets positive control distances for all genes in an operon
	curs_pos = db.cursor()
	command = 'select g.gene_id, e.left_pos, e.right_pos, g.strand from genes g join exons e using(gene_id) where g.locus_tag in (\'' + '\',\''.join(bOperonGenes) + '\') order by e.left_pos asc;'
	curs_pos.execute(command)
	operonGeneInfo =  curs_pos.fetchall()
	#Accesses the distances
	length = len(bOperonGenes)
	for i in range(1, length):
		firstRight = operonGeneInfo[i-1][2]
		secondLeft = operonGeneInfo[i][1]
		dist = secondLeft - firstRight +1
		posCtrlDist.append(dist)

#Creates a table of gene id's, bNumbers and left/right positions
command = 'select g.gene_id, g.strand, g.locus_tag, e.left_pos, e.right_pos from genes g join exons e using(gene_id) order by e.left_pos asc;'
curs_join = db.cursor()
curs_join.execute(command)
joinTable = curs_join.fetchall()
for line in curs_join:
	bNum = line[2]
#Creates a table that accesses all the operons
command = 'select name, strand, left_pos, right_pos from Operons order by left_pos asc;'
curs_operon = db.cursor()
curs_operon.execute(command)
operonsTable = curs_operon.fetchall()
#Loops through operons to get the negative control distances
i = 0
count = 0
while i < len(operonsTable):
	#Accesses the operon names
	operon = operonsTable[i][0]
	#operon = operonsTable[i+1][0]
	#Only acesses opeorns that contain genes associated with a bNum
	if (operon in geneID):
		#Gets the list of exon coordinates for genes in the first operon
		geneList = geneID[operon]
		command = 'select g.locus_tag, e.left_pos, e.right_pos, g.strand from genes g join exons e using(gene_id) where g.locus_tag in (\''+ '\',\''.join(geneList) + '\') order by e.right_pos desc limit 1;'
		curs_prev = db.cursor()
		curs_prev.execute(command)
		prev = curs_prev.fetchall()[0]
		strandPrev = prev[3]
		maxRight = prev[2]
		#Finds the closest command
		curs_closestGene = db.cursor()
		command = 'select g.strand, e.left_pos from genes g join exons e using(gene_id) where e.left_pos > ' + str(maxRight) + ' order by left_pos asc limit 1;'
		curs_closestGene.execute(command)
		curr = curs_closestGene.fetchall()[0]
		strandCurr = curr[0]
		minLeft = curr[1]
		#Stores Negative Control Distances
		if strandPrev == strandCurr:
			dist = minLeft - maxRight + 1
			negCtrlDist.append(dist)
	i += 1
print "\nThe number of positive controls analyzed is:", len(posCtrlDist)
print "\nThe number of negatibe controls analyzed is:", len(negCtrlDist), '\n'
#Probability determination
posProb = [str(float(num)/len(posCtrlDist)) for num in posCtrlDist]
negProb = [str(float(num)/len(negCtrlDist)) for num in posCtrlDist]
negProbF = open( 'negativeCtrlProb.txt', 'w')
posProbF = open('positiveCtrlProb.txt', 'w')
negProbF.write('\n'.join(negProb) )
posProbF.write('\n'.join(posProb) )
#Output diostance values to another file because matplotlib is not supported in ieng6
posCtrl = open('posCtrlDist.txt', 'w')
negCtrl = open('negCtrlDist.txt','w')
posCtrlDist = [ str(num) for num in posCtrlDist ]
negCtrlDist = [str(num) for num in negCtrlDist ]
posCtrl.write('\n'.join(posCtrlDist))
negCtrl.write('\n'.join(negCtrlDist))
#outputs a table of bNumbers
operonDict = open('operonBNums.txt', 'w')
for operon in geneID:
	for gene in geneID[operon]:
		operonDict.write( operon + '\t' + gene + '\n' )

