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


