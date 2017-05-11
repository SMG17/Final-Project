# Load libraries:
library(DESeq2)
library(GenomeInfoDb)
library(tximport)
library(rhdf5)
library(data.table)

# Use DESeq2 for differential mRNA expression:
# Create a named vector of files to import:
files_RNA <- c("C:/Users/SMG/Desktop/Sequencing_class/Final/Scer_RNA_output0/abundance.tsv",
               "C:/Users/SMG/Desktop/Sequencing_class/Final/Scer_RNA_output1/abundance.tsv",
               "C:/Users/SMG/Desktop/Sequencing_class/Final/Spar_RNA_output0/abundance.tsv",
               "C:/Users/SMG/Desktop/Sequencing_class/Final/Spar_RNA_output1/abundance.tsv")
names(files_RNA) <- c("Scer1","Scer2","Spar1","Spar2")

# Read in the Kallisto files using tximport:
txdat_RNA <- tximport(files_RNA, type = "kallisto", txOut = TRUE)

# Generate the condition matrix:
coldata_RNA <- data.frame(condition = c("Scer","Scer","Spar","Spar"))
rownames(coldata_RNA) = names(files_RNA)
coldata_RNA

# Turn it into a DESeq2 object:
dds_RNA <- DESeqDataSetFromTximport(txdat_RNA, colData = coldata_RNA, design =~ condition) #Tells it to group samples by the "condition" column of colData

# Run differential expression:
dds_RNA <- DESeq(dds_RNA)

# Summarize results:
res_RNA <- results(dds_RNA)
res_RNA_FDR10 <- subset(res_RNA, padj < 0.1)
write.csv(as.data.frame(res_RNA_FDR10[order(res_RNA_FDR10$log2FoldChange),]), file = "C:/Users/SMG/Desktop/Sequencing_class/Final/mRNA_diff_exp_sig.csv")

res_RNA_table <- as.data.table(as.data.frame(res_RNA), keep.rownames="id")
res_RNA_raw <- res_RNA_table[,.(id,log2FoldChange)]
write.table(res_RNA_raw, file = "C:/Users/SMG/Desktop/Sequencing_class/Final/mRNA_raw.txt", col.name=FALSE, row.name=FALSE)

res_RNA_sig <- res_RNA_table[padj< 0.1,.(id,log2FoldChange)]
write.table(res_RNA_sig, file = "C:/Users/SMG/Desktop/Sequencing_class/Final/mRNA_diff_exp_sig.txt", col.name=FALSE, row.name=FALSE)

# A plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts. Points will be colored red if the adjusted p-value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.
# An "MA-plot" provides a useful overview for an experiment with a two-group comparison (here the two conditions: WT and SNF2). The log2 fold change (using the mean of the counts for all replicates within WT, and within SNF2) for a particular comparison is plotted on the y-axis and the average of the counts normalized by size factor (by group/condition) is shown on the x-axis ("M" for minus, because a log ratio is equal to log minus log, and "A" for average).
# Each gene is represented with a dot. Genes with an adjusted p-value below a threshold (here 0.1, the default) are shown in red. The DESeq2 package incorporates a prior on log2 fold changes, resulting in moderated log2 fold changes from genes with low counts and highly variable counts, as can be seen by the narrowing of spread of points on the left side of the plot. This plot demonstrates that only genes with a large average normalized count contain sufficient information to yield a significant call.
plotMA(res_RNA) #, main="DESeq2", ylim=c(-2,2)

# Whether a gene is called significant depends not only on its LFC but also on its within-group variability, which DESeq2 quantifies as the dispersion. For strongly expressed genes, the dispersion can be understood as a squared coefficient of variation: a dispersion value of 0.01 means that the gene's expression tends to differ by typically ???0.01 = 10% between samples of the same treatment group. For weak genes, the Poisson noise is an additional source of noise.
# The function plotDispEsts visualizes DESeq2's dispersion estimates: each dot is a gene, and on the y-axis shows how much variation there is in how the gene is expressed in the group/condition.
# The black points are the dispersion estimates for each gene as obtained by considering the information from each gene separately. Unless one has many samples, these values fluctuate strongly around their true values. Therefore, we fit the red trend line, which shows the dispersions' dependence on the mean, and then shrink each gene's estimate towards the red line to obtain the final estimates (blue points) that are then used in the hypothesis test. The blue circles above the main "cloud" of points are genes which have high gene-wise dispersion estimates which are labelled as dispersion outliers. These estimates are therefore not shrunk toward the fitted trend line.
plotDispEsts(dds_RNA)

# Calculate sum of p-value below 0.05:
nom_pvalue_nb_RNA <- sum(res_RNA$pvalue < 0.1, na.rm = TRUE) # Should missing values (including NaN) be removed?
nom_pvalue_nb_RNA
# 2883; Kristy had 2884

# Calculate sum of adjusted p-values:
false_discovery_nb_RNA <- sum(res_RNA$padj < 0.1, na.rm = TRUE)
false_discovery_nb_RNA
# 2472; Kristy had 2469


# Use DESeq2 for differential ribosome occupancy (basically number of ribosomes working on mRNA, so measure of degree of protein synthesis):
# Create a named vector of files to import:
files_Ribo <- c("C:/Users/SMG/Desktop/Sequencing_class/Final/Scer_ribo_output0/abundance.tsv",
                "C:/Users/SMG/Desktop/Sequencing_class/Final/Scer_ribo_output1/abundance.tsv",
                "C:/Users/SMG/Desktop/Sequencing_class/Final/Spar_ribo_output0/abundance.tsv",
                "C:/Users/SMG/Desktop/Sequencing_class/Final/Spar_ribo_output1/abundance.tsv")
names(files_Ribo) <- c("Scer1","Scer2","Spar1","Spar2")

# Read in the Kallisto files using tximport:
txdat_Ribo <- tximport(files_Ribo, type = "kallisto", txOut = TRUE)

# Generate the condition matrix:
coldata_Ribo <- data.frame(condition = c("Scer","Scer","Spar","Spar"))
rownames(coldata_Ribo) = names(files_Ribo)
coldata_Ribo

# Turn it into a DESeq2 object:
dds_Ribo <- DESeqDataSetFromTximport(txdat_Ribo, colData = coldata_Ribo, design =~ condition) #Tells it to group samples by the "condition" column of colData

# Run differential expression:
dds_Ribo <- DESeq(dds_Ribo)

# Summarize results:
res_Ribo <- results(dds_Ribo)
res_Ribo_FDR10 <- subset(res_Ribo, padj < 0.1)
write.csv(as.data.frame(res_Ribo_FDR10[order(res_Ribo_FDR10$log2FoldChange),]), file = "C:/Users/SMG/Desktop/Sequencing_class/Final/Ribo_diff_occup_sig.csv")

# MA-plot:
plotMA(res_Ribo) #, main="DESeq2", ylim=c(-2,2)

# Dispersion plot:
plotDispEsts(dds_Ribo)

# Calculate sum of p-value below 0.05:
nom_pvalue_nb_Ribo <- sum(res_Ribo$pvalue < 0.1, na.rm = TRUE) # Should missing values (including NaN) be removed?
nom_pvalue_nb_Ribo
# 1179; Kristy has 1178

# Calculate sum of adjusted p-values:
false_discovery_nb_Ribo <- sum(res_Ribo$padj < 0.1, na.rm = TRUE)
false_discovery_nb_Ribo
# 556; Kristy has 563

# Use DESeq2 to calculate translational efficiency is the ribosome occupancy (counts from ribosome profiling experiment) normalized by the mRNA abundance (counts from mRNA experiment) (to disentangle whether degree of protein synthesis results from a few ribosomes working on many mRNAs, or a lot of ribosomes working on a few mRNAs):
files_TE <- c("C:/Users/SMG/Desktop/Sequencing_class/Final/Scer_RNA_output0/abundance.tsv",
                 "C:/Users/SMG/Desktop/Sequencing_class/Final/Scer_RNA_output1/abundance.tsv",
                 "C:/Users/SMG/Desktop/Sequencing_class/Final/Spar_RNA_output0/abundance.tsv",
                 "C:/Users/SMG/Desktop/Sequencing_class/Final/Spar_RNA_output1/abundance.tsv",
                 "C:/Users/SMG/Desktop/Sequencing_class/Final/Scer_ribo_output0/abundance.tsv",
                 "C:/Users/SMG/Desktop/Sequencing_class/Final/Scer_ribo_output1/abundance.tsv",
                 "C:/Users/SMG/Desktop/Sequencing_class/Final/Spar_ribo_output0/abundance.tsv",
                 "C:/Users/SMG/Desktop/Sequencing_class/Final/Spar_ribo_output1/abundance.tsv")
names(files_TE) <- c("ScerRNA1","ScerRNA2","SparRNA1","SparRNA2","ScerRibo1","ScerRibo2","SparRibo1","SparRibo2")

# Read in the Kallisto files using tximport:
txdat_TE <- tximport(files_TE, type = "kallisto", txOut = TRUE)

# Generate the condition matrix:
coldata_TE<- data.frame(condition = c("Scer","Scer","Spar","Spar","Scer","Scer","Spar","Spar"), assay = c("RNA","RNA","RNA","RNA","Ribo","Ribo","Ribo","Ribo"))
rownames(coldata_TE) = names(files_TE)
coldata_TE

# Turn it into a DESeq2 object:
dds_TE <- DESeqDataSetFromTximport(txdat_TE, colData = coldata_TE, design= ~ assay + condition + assay:condition)

# Run differential expression:
dds_TE <- DESeq(dds_TE, test="LRT", reduced= ~ assay + condition) # likelihood ratio test = comparing the normal hypothesis to the alternative hypothesis, Ha states that the two expression (or dispersion?) values are different; H0: they are the same, so for to test it against the H0, we need to remove the interaction term.

# Summarize results:
res_TE <- results(dds_TE)
res_TE_FDR10 <- subset(res_TE, padj < 0.1)
write.csv(as.data.frame(res_TE_FDR10[order(res_TE_FDR10$log2FoldChange),]), file = "C:/Users/SMG/Desktop/Sequencing_class/Final/TE_diff_sig.csv")

res_TE_table <- as.data.table(as.data.frame(res_TE), keep.rownames="id")
res_TE_raw <- res_TE_table[,.(id,log2FoldChange)]
write.table(res_TE_raw, file = "C:/Users/SMG/Desktop/Sequencing_class/Final/TE_raw.txt", col.name=FALSE, row.name=FALSE)

res_TE_sig <- res_TE_table[padj< 0.1,.(id,log2FoldChange)]
write.table(res_TE_sig, file = "C:/Users/SMG/Desktop/Sequencing_class/Final/TE_sig.txt", col.name=FALSE, row.name=FALSE)

# MA-plot:
plotMA(res_TE) #, main="DESeq2", ylim=c(-2,2)

# Dispersion plot:
plotDispEsts(dds_TE)

# Calculate sum of p-value below 0.05:
nom_pvalue_nb_TE<- sum(res_TE$pvalue < 0.1, na.rm = TRUE) # Should missing values (including NaN) be removed?
nom_pvalue_nb_TE
# 1188; same as Kristy

# Calculate sum of adjusted p-values:
false_discovery_nb_TE <- sum(res_TE$padj < 0.1, na.rm = TRUE)
false_discovery_nb_TE
# 332; same as Kristy

# Scatterplot showing log fold change in mRNA abundance between the species vs. log fold change in translational efficiency.
plot(res_RNA$log2FoldChange, res_TE$log2FoldChange, main = "Log fold change in mRNA abundance vs. log fold change in translation efficiency", xlab= "Log fold change in mRNA abundance", ylab = "Log fold change in translation efficiency")
#plot(RNA_raw[,2],TE_raw[,2],main = "Log fold change in mRNA abundance vs. log fold change in translation efficiency", xlab= "Log fold change in mRNA abundance", ylab = "Log fold change in translation efficiency")

# Determining compensatory evolution (mRNA up, TE down, or mRNA down, TE up) and coordinated evolution (mRNA up, TE up, or mRNA down, TE down)
RNA_raw <- read.table("C:/Users/SMG/Desktop/Sequencing_class/Final/mRNA_raw.txt")
TE_raw <- read.table("C:/Users/SMG/Desktop/Sequencing_class/Final/TE_raw.txt")
table <- as.data.table(merge(RNA_raw,TE_raw,by="V1"))

RNAup_TEdown <- length(which(table$V2.x>0&table$V2.y<0))
RNAup_TEdown
# 892
RNAdown_TEup <- length(which(table$V2.x<0&table$V2.y>0))
RNAdown_TEup
# 1024
RNAup_TEup <- length(which(table$V2.x>0&table$V2.y>0))
RNAup_TEup
# 1815
RNAdown_TEdown <- length(which(table$V2.x<0&table$V2.y<0))
RNAdown_TEdown
# 1707

RNA_sig <- read.table("C:/Users/SMG/Desktop/Sequencing_class/Final/mRNA_diff_exp_sig.txt")
TE_sig <- read.table("C:/Users/SMG/Desktop/Sequencing_class/Final/TE_sig.txt")
table_sig <- as.data.table(merge(RNA_sig,TE_sig,by="V1", all=FALSE))
table_sig_num <- transform(table_sig, V2.x = as.numeric(V2.x), V2.y = as.numeric(V2.y))
sapply(table_sig_num, mode)

RNAup_TEdown_sig=length(which(table_sig_num$V2.x>0&table_sig_num$V2.y<0))
RNAup_TEdown_sig
# 15
RNAdown_TEup_sig=length(which(table_sig_num$V2.x<0&table_sig_num$V2.y>0))
RNAdown_TEup_sig
# 27
RNAup_TEup_sig=length(which(table_sig_num$V2.x>0&table_sig_num$V2.y>0))
RNAup_TEup_sig
# 169
RNAdown_TEdown_sig=length(which(table_sig_num$V2.x<0&table_sig_num$V2.y<0))
RNAdown_TEdown_sig
# 47