# Reading genome reference fasta files:
def read_fasta(fastaFile):
	#initialize empty dictionary
	fastaDict = {}
	#go over every line in the file
	chromNum = 0
	for line in fastaFile:
		#if it's got a > at the start, we know it's the name of a sequence
		if line[0] == ">":
			chromNum += 1
			chromName = chromNum #line.strip()[1:]
			#initialize newly seen sequence names with empty lists 
			#this way, we can append every line as it comes in
			fastaDict[chromName] = []
		else:
			#append all the lines that belong to that chromosome
			fastaDict[chromName].append(line.strip())
	#now we need to merge all the lines together
	#loop over every chromosome and merge all of its lines
	for chrom in fastaDict:
		#''.join(list_of_strings) will join all the strings in the list into a new string with no separating characters
		fastaDict[chrom] = ''.join(fastaDict[chrom])
	return fastaDict

x = open("/mnt/c/Users/SMG/Desktop/Final/S_cerevisiae.fa")
ScerGenomeDict = read_fasta(x)
output = open("Scer_transcriptome.fa","w")

bedFileScer = open("/mnt/c/Users/SMG/Desktop/Final/S_cerevisiae_genes.bed")
for line in bedFileScer:
	line = line.strip().split()
	if line[0] == 'chrI':
		line.insert(0,1)
	if line[0] == 'chrII':
                line.insert(0,2)
	if line[0] == 'chrIII':
                line.insert(0,3)
	if line[0] == 'chrIV':
                line.insert(0,4)
	if line[0] == 'chrV':
                line.insert(0,5)
	if line[0] == 'chrVI':
                line.insert(0,6)
	if line[0] == 'chrVII':
                line.insert(0,7)
	if line[0] == 'chrVIII':
                line.insert(0,8)
	if line[0] == 'chrIX':
                line.insert(0,9)
	if line[0] == 'chrX':
                line.insert(0,10)
	if line[0] == 'chrXI':
                line.insert(0,11)
	if line[0] == 'chrXII':
                line.insert(0,12)
	if line[0] == 'chrXIII':
                line.insert(0,13)
	if line[0] == 'chrXIV':
                line.insert(0,14)
	if line[0] == 'chrXV':
                line.insert(0,15)
	if line[0] == 'chrXVI':
                line.insert(0,16)
	geneName = line[4]
	geneSeq = ScerGenomeDict[line[0]][int(line[2]):int(line[3])+1]
	output.write(">%s\n%s\n" % (geneName, geneSeq))
	
# Define function to give complementary nucleotides
def complement(nuc):
	if nuc == "A":
		return "T"
	elif nuc == "C":
		return "G"
	elif nuc == "G":
		return "C"
	elif nuc == "T":
		return "A"

def reverse_complement(seq):
	#first reverse the sequence using a weird python trick
	revSeq = seq[::-1]
	#now go through and complement everything
	revCom = []
	for pos in revSeq:
		curComplement = complement(pos)
		revCom.append(curComplement)
	#now join the reverse complement back together
	revComSeq = ''.join(revCom)
	return revComSeq

	
def fasta_read(x):
	names_line_num = []
	names = []
	genome = []
	chromosome_number = 0
	for i, line in enumerate(x):
		line = line.strip()
		genome.append(line)
		if line.startswith(">"):
			chromosome_number += 1
			names_line_num.append(i)
			chromosome = line[1:]
			chromosome = chromosome.replace("_"," ")
			names.append(chromosome)

	names_line_num.append(len(genome))
	sequences = []
	for i in range(len(names_line_num)-1):
		seq = [''.join(genome[names_line_num[i]+1:names_line_num[i+1]])]
		sequences.extend(seq)

	dictionary = dict(zip(names, sequences))
	return dictionary

