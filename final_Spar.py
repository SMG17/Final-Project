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

x = open("/mnt/c/Users/SMG/Desktop/Sequencing_class/Final/S_paradoxus.fa")
SparGenomeDict = read_fasta(x)
output = open("Spar_transcriptome.fa","w")

bedFile_Spar = open("/mnt/c/Users/SMG/Desktop/Sequencing_class/Final/S_paradoxus_genes.bed")
for line in bedFile_Spar:
	line = line.strip().split()
	if line[0] == 'Spar_1':
		line.insert(0,1)
	if line[0] == 'Spar_2':
                line.insert(0,2)
	if line[0] == 'Spar_3':
                line.insert(0,3)
	if line[0] == 'Spar_4':
                line.insert(0,4)
	if line[0] == 'Spar_5':
                line.insert(0,5)
	if line[0] == 'Spar_6':
                line.insert(0,6)
	if line[0] == 'Spar_7':
                line.insert(0,7)
	if line[0] == 'Spar_8':
                line.insert(0,8)
	if line[0] == 'Spar_9':
                line.insert(0,9)
	if line[0] == 'Spar_10':
                line.insert(0,10)
	if line[0] == 'Spar_11':
                line.insert(0,11)
	if line[0] == 'Spar_12':
                line.insert(0,12)
	if line[0] == 'Spar_13':
                line.insert(0,13)
	if line[0] == 'Spar_14':
                line.insert(0,14)
	if line[0] == 'Spar_15':
                line.insert(0,15)
	if line[0] == 'Spar_16':
                line.insert(0,16)
	stop = line[4].find('.')
	if stop == -1:
		geneName = line[4]
	else:
		geneName = line[4][:stop]
	geneSeq = SparGenomeDict[line[0]][int(line[2]):(int(line[3])+1)]
	output.write(">%s\n%s\n" % (geneName,geneSeq))

x.close()
output.close()
