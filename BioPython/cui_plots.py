import matplotlib.pyplot as plt
sortedCUI = open('../Week2/cui_freq.txt','r')
sortedChr = open('../Week2/cui_freq.tsv','r')
#Maps Hist of genes by CUI
names = []
cui = []
for line in sortedCUI:
	ln = line.strip().split()
	names.append(ln[0])
	cui.append(ln[1])
plt.scatter(range(0, len(cui)), cui, marker = 'o', color = 'c' )
plt.title("Genes Mapped by Descending CUI")
plt.xlabel("Genes")
plt.ylabel("CUI")
plt.show()
#Maps hist of genes by chromosome
names = []
cui = []
for line in sortedChr:
	ln = line.strip().split()
	names.append(ln[0])
	cui.append(ln[1])
plt.scatter(range(0, len(cui)), cui, marker = 'o', color = 'c')
plt.title("Genes CUI Mapped by Chromosomal Order")
plt.xlabel("Genes")
plt.ylabel("CUI")
plt.show()
