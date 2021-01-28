# Code for de novo assembly

## Install Packages

```{BASH}
conda create -n de_novo -c bioconda -c conda-forge fastqc=0.11.5 \
             trimmomatic=0.36 spades=3.11.1 quast=5.0.2 \
             bowtie2=2.2.5 java-jdk=8.0.112 samtools --yes
```

## Activate the environment
```{BASH}
conda activate de_novo
```


## Fastqc
```{BASH}
fastqc 
```

## Trimmomatic
```{BASH} 
trimmomatic-0.36.jar PE -phred33 G11_S2_L001_R1_001.fastq.gz G11_S2_L001_R2_001.fastq.gz G11_R1_paired.fq.gz G11_R1_unpaired.fq.gz G11_R2_paired.fq.gz G11_R2_unpaired.fq.gz ILLUMINACLIP:/softwares/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

## SPAdes
```{BASH}
spades.py -k 21,51,71,91,111,127 --careful --pe1-1 T1_R1.cutadapt.fastq --pe1-2 T1_R2.cutadapt.fastq -o T1_spades_output
```

## Prokka

#### Format the SPAdes assembly to be compatible with prokka
Prokka does not like long names as is the natural output of SPAdes.  To prepare our data for processing with prokka, we must change the names of the contigs.  We will use a unix filter known as `awk` to perform this change.
```{BASH}
awk '/^>/{print ">T1_" ++i; next}{print}' < contigs.fasta > contigs_names.fasta
```
This command is performing the following function.
-  search for lines that start with `>`
-  replace the `>` with `>T1_`
-  after the `T1_` add a number with the first instance being 1 and each instance increasing by 1
-  read in `contigs.fasta`
-  output a new file `contigs_names.fasta

#### Run the prokka pipeline on the `contigs_names.fasta` file.
```{BASH}
prokka --outdir T1 --prefix T1 contigs_names.fasta
```
