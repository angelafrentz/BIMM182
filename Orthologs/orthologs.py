import MySQLdb as sql
db = sql.connect( host = 'bm185s-mysql.ucsd.edu', user = 'afrentze', passwd = 'Sql$pwd', db = 'afrentze_db')
tableList = ['Blast_EvsA', 'Blast_AvsE']
outfile = open('orthologsTablei.txt','w')
for i in range(0, len(tableList)):	
	tablename = tableList[i]
	otherTable = tableList[0]
	if i == 0:
		otherTable = tableList[1]
	cursor = db.cursor()
	cursor.execute('SELECT DISTINCT(qseqid) FROM ' + tablename + ' WHERE qcovs>=60 OR scov >=60;');
	distinct = cursor.fetchall()
	for line in distinct:
		uniqID = line[0]
		curs = db.cursor()
		curs.execute('select qseqid, sseqid, bitscore from ' + tablename + ' where qseqid = \'' + uniqID + '\' order by bitscore desc limit 1;' )
		sorted = curs.fetchall()
		highestMatch1 = sorted[0][1]
		curs = db.cursor()
		curs.execute('select qseqid, sseqid, bitscore from ' + otherTable + ' where qseqid = \'' + highestMatch1 + '\' order by bitscore desc limit 1;' )
		sorted = curs.fetchall();
		if sorted:
			highestMatch2 = sorted[0][1]
		if uniqID == highestMatch2:
			#command =  'insert into orthologs(id1, id2, type) values (\'' + uniqID + '\', \'' + highestMatch1 + '\', \'ortholog\');'
			#print command
			outfile.write(uniqID + '\t' + highestMatch1 + '\t' + 'ortholog' + '\n')
			
db.close()


