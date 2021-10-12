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

Make a sample.csv file (comma separated)

```{R}
samples <- read.csv("sample.csv", header = TRUE)
files <- file.path("quant_M", samples$sample, "quant.sf")
names(files) <- paste0(samples$sample)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
```
N.B. If you don't import annotations rmove `, tx2gene = tx2gene`




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
