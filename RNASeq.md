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

 
 
