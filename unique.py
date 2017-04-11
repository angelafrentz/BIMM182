import sys
import bz2

filename = sys.argv[1]
protInfo = bz2.BZ2File(filename)
protDict = {}
for line in protInfo:
	protLN = line.strip().split()
	if (len(protDict.keys()) == 2000) and (protLN[0] not in protDict):
		break
	protDict.setdefault(protLN[0], []).append(protLN[3])
#Protein,Num Interactions, Max Interaction
protInfo = []
maxI = -float('inf')
maxProt = None
for prot in protDict:
	numI = len(protDict[prot])
	maxVal = max(protDict[prot])
	protInfo.append([prot, numI, maxVal])
	if numI > maxI:
		maxI = numI
		maxProt = prot
print maxProt, maxI
