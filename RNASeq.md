## Install packages (This is more than you need)
```{BASH}
conda create -n analysis -c bioconda -c conda-forge fastqc \
             trimmomatic spades megahit quast \
             bowtie2 java-jdk samtools muscle hmmer cufflinks-py tophat2 salmon --yes
 ```
 
 ## Activate the environment
 ```{BASH}
 conda activate analysis
 ```
 
## Quality Trim the reads
 
### Check the quality with fastqc.  

Select a subset of the reads
 
 ```{BASH}
 fastqc Reads.fastqc.gz
 ```

### Quality trim the reads with trimmomatic

```{BASH}
trimmomatic PE -phred33 G11_S2_L001_R1_001.fastq.gz G11_S2_L001_R2_001.fastq.gz G11_R1_paired.fq.gz G11_R1_unpaired.fq.gz G11_R2_paired.fq.gz G11_R2_unpaired.fq.gz ILLUMINACLIP:~/miniconda3/pkgs/trimmomatic-0.36-6/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
## Mapping of reads with `salmon`

### Making a quantification directory
```{BASH}
mkdir quants
```

### Index the genome sequences for mapping 
```{BASH}
salmon index -t StrainM_transcripts.fasta -i StrainM_index
```

### Map and count reads to the index using salmon.
```{BASH}
salmon quant -i StrainM_index -l A -1 M1_A_Paired_R1.fq.gz -2 M1_A_Paired_R2.fq.gz -p 8 --validateMappings -o  quants/M1_A_quant
```
## Import into R THIS IS NOW ALL IN R ON THE CONSOLE
based on https://www.hadriengourle.com/tutorials/rna/

### Set working directory.

1. Session
2. Set working directory
3. Choose directory
4. Navigate to Biocide_RNASeq

### Install packages

```{R}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("GenomicFeatures")
BiocManager::install("readr")
```

### Initialize packages
```{R}
library(tximport)
library(GenomicFeatures)
library(readr)
```

### Upload quant data

Make a sample.csv file (comma separated file with samples as one column)

```{R}
samples <- read.csv("sample.csv", header = TRUE)
files <- file.path("quant_M", samples$sample, "quant.sf")
names(files) <- paste0(samples$sample)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
```
N.B. If you don't import annotations rmove `, tx2gene = tx2gene`

### Preparing for DESeq
```{BASH}
dds <- DESeqDataSetFromTximport(txi.salmon, samples, ~condition)
```
Subset the DESeq object
```{BASH}
dds_C_G<- dds[,c(columns_to_keep] ## Columns to keep is a comma separated list of the columns in the comparisions you're interested in.
```

Run DESeq2
```{BASH}
dds_C_G <- DESeq(dds_C_G, test="Wald", fitType="parametric")
```
Export the results
```{BASH}
res = results(dds_C_G, cooksCutoff = FALSE)
```
Set Significance cutoff
```{BASH}
alpha = 0.05
```
Export significant results
```{BASH}
sigtab = res[which(res$padj < alpha), ]
```

# Alternative Path for Quant

### Upload quant data

Make a sample.csv file (comma separated file with samples as one column)

**Check to make sure I labeled these correctly**
```{R}
samples <- read.csv("sample.csv", header = TRUE)
samples_CG<-filter(samples,condition=="control"|condition="GA")
samples_CD<-filter(samples,condition=="control"|condition="DBNPA")
samples_CB<-filter(samples,condition=="control"|condition="BAC")
```
**I can't remember if the header for the sample names is sample or names. So change to what you have in your R script**
```{R}
files_CG <- file.path("quant_M", samples_CG$sample, "quant.sf")
files_CD <- file.path("quant_M", samples_CD$sample, "quant.sf")
files_CB <- file.path("quant_M", samples_CB$sample, "quant.sf")
names(files_CG) <- paste0(samples_CG$sample)
names(files_CD) <- paste0(samples_CD$sample)
names(files_CB) <- paste0(samples_CB$sample)
```
Make 3 import files
```{R} 
txi.salmon_CG <- tximport(files_CG, type = "salmon", tx2gene = tx2gene)
txi.salmon_CD <- tximport(files_CD, type = "salmon", tx2gene = tx2gene)
txi.salmon_CB <- tximport(files_CB, type = "salmon", tx2gene = tx2gene)
```
N.B. If you don't import annotations rmove `, tx2gene = tx2gene`

Then import these three files separatly into DESeq



# Plotting the DESeq Results

## Add the gene annotation data to the `res` object

## Make Volcano plot

### Add a label of if the genes are significant or not.
```{R}
All_res<-dplyr::mutate(All_res,treatment = ifelse(log2FoldChange<=-2 & padj<=0.01,"Biocide",
                                        ifelse (log2FoldChange>=2 & padj<=0.01, "Control", "Non-significant")))
```
### Make a volcano plot
```{R}
ggplot(All_res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=significance)) +
  ggtitle("Control versus Biocide") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  scale_fill_manual(values=c("blue", "orange", "red"))+
  geom_hline(yintercept = 2, colour="#990000", linetype="dashed") + geom_vline(xintercept = 2, colour="#990000", linetype="dashed") + geom_vline(xintercept = -2, colour="#990000", linetype="dashed")+
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  +
  theme_bw()
```

## Plot PCA
You can do this with each of the pairs or you could this with the full dataset
### Log transform the data
```{R}
rld=rlog(dds_C_G)
```
### Make PCA
```{R}
plotPCA(rld, intgroup = "condition")
```

### Extract information from the PCA file

```{R}
pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
```

### Make a Pretty PCA
```{R}
ggplot(pcaData, aes(PC1, PC2, color=dex, shape=dex)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme_bw()
```

## Plot expression of individual genes

## Plotting Heatmaps.
**I've not run this before. So it may not work**

```{R}
library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", 
     trace="none", dendrogram="column", 
     col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
     ColSideColors = c( Control="gray", GA="darkgreen", BAC="orange", DBNPA="darkblue" )[
        colData(rld)$condition ] )
```

# Old Don't Use

## Mapping of reads with `bowtie2`
### Index reference database
- Use bowtie2 to index your reference genome to make it possible for you to search it
```{BASH}
bowtie2-build C.hydrogenoformans_Z2901.fasta C.hydro
```

### Map reads your indexed reference database
- Use bowtie 2 to map your reads to your indexed reference genome
```{BASH}
bowtie2 -x C.hydro -1 C.hydro_6008_R1.fastq -2 C.hydro_6008_R2.fastq -U C.hydro_6008_unpaired.fastq -S 6008.sam
```
### Add annotations with Salmon
Upload the .gff file to the Biocide_RNASeq directory

N.B. Try this, but if it gives errors move on.
```{R}
txdb_M <- makeTxDbFromGFF("Strain-M.gff")
k <- keys(txdb_M, keytype = "ID")
tx2gene <- select(txdb_M, keys = k, keytype = "ID", columns = "TXNAME")
head(tx2gene)
``` 


