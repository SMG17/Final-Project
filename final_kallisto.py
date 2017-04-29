import os

ScerRiboSamples = ["Scer_ribo_seq_1.fastq.gz","Scer_ribo_seq_2.fastq.gz"]

#for file in os.listdir("/mnt/c/Users/SMG/Desktop/Sequencing_class/Final/Scer_ribo_seq/"):
#    if file.endswith(".fastq.gz"):
#        ScerRiboSamples.append(file)

for sample in ScerRiboSamples:
	command = "kallisto quant -i Scer_transcriptome.idx -o Scer_ribo_output%s --single -l 180 -s 20 %s" % (ScerRiboSamples.index(sample), sample)
	print "Currently running: %s" % command
	os.system(command) 

ScerRNASamples = ["Scer_RNA_seq_1.fastq.gz","Scer_RNA_seq_2.fastq.gz"]

#for file in os.listdir("/mnt/c/Users/SMG/Desktop/Sequencing_class/Final/Scer_RNA_seq/"):
#    if file.endswith(".fastq.gz"):
#        ScerRNASamples.append(file)

for sample in ScerRNASamples:
        command = "kallisto quant -i Scer_transcriptome.idx -o Scer_RNA_output%s --single -l 180 -s 20 %s" % (ScerRNASamples.index(sample), sample)
        print "Currently running: %s" % command
        os.system(command)

SparRiboSamples = ["Spar_ribo_seq_1.fastq.gz","Spar_ribo_seq_2.fastq.gz"]

#for file in os.listdir("/mnt/c/Users/SMG/Desktop/Sequencing_class/Final/Spar_ribo_seq/"):
#    if file.endswith(".fastq.gz"):
#        SparRiboSamples.append(file)

for sample in SparRiboSamples:
        command = "kallisto quant -i Spar_transcriptome.idx -o Spar_ribo_output%s --single -l 180 -s 20 %s" % (SparRiboSamples.index(sample), sample)
        print "Currently running: %s" % command
        os.system(command)

SparRNASamples = ["Spar_RNA_seq_1.fastq.gz","Spar_RNA_seq_2.fastq.gz"]

#for file in os.listdir("/mnt/c/Users/SMG/Desktop/Sequencing_class/Final/Spar_RNA_seq/"):
#    if file.endswith(".fastq.gz"):
#        SparRNASamples.append(file)

for sample in SparRNASamples:
        command = "kallisto quant -i Spar_transcriptome.idx -o Spar_RNA_output%s --single -l 180 -s 20 %s" % (SparRNASamples.index(sample), sample)
        print "Currently running: %s" % command
        os.system(command)
