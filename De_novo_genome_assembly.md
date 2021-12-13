# Code for de novo assembly
## Prepare the package manager

Log on to the smtechtm-pc.sabu.mtu.edu using putty

```{BASH}
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
- Make sure to look through the miniconda EULA and agree by typing `yes`  
- Install in your home directory by hitting return  
- Select yes to initialize the installation

Then activate miniconda
```{BASH}
source ~/.bashrc
```
You should see the command prompt have a `(base)` before your `username@smtechtm-pc`

## Install Packages

```{BASH}
conda create -n de_novo -c bioconda -c conda-forge fastqc=0.11.5 \
             trimmomatic spades megahit quast \
             bowtie2=2.2.5 java-jdk=8.0.112 prokka=1.14.6 samtools --yes
```

## Activate the environment
To access all of those packages you must activate your environment.

```{BASH}
conda activate de_novo
```

# Trimming

To trim off low quality reads you will need to run trimmomatic.

## Fastqc
You can first check the quality of your reads. You don't need to do this with everyone, but it might be good to check a few.

```{BASH}
fastqc G11_S2_L001_R1_001.fastq.gz
```

## Trimmomatic
Now you'll need to trim your reads.

### Make a trimming directory for a bit of housekeeping
```{BASH}
mkdir Trimming
```
### Perform the trimming

```{BASH} 
trimmomatic PE -phred33 G11_S2_L001_R1_001.fastq.gz G11_S2_L001_R2_001.fastq.gz Trimming/G11_R1_paired.fq.gz Trimming/G11_R1_unpaired.fq.gz Trimming/G11_R2_paired.fq.gz Trimming/G11_R2_unpaired.fq.gz ILLUMINACLIP:~/miniconda3/pkgs/trimmomatic-0.36-6/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
Here are details about the script above.
`G11_S2_L001_R1_001.fastq.gz` - This is your forward read for your raw reads
`G11_S2_L001_R2_001.fastq.gz` - This is your reverse read for your raw reads 
`G11_R1_paired.fq.gz` - This is the new forward paired read
`G11_R1_unpaired.fq.gz` - This is the new forward unpaired read
`G11_R2_paired.fq.gz` - This is the new reverse paired read
`G11_R2_unpaired.fq.gz` - This is the new reverse unpaired read
`ILLUMINACLIP:/home/pkokate/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10`  - This is the path to the location of a file that contains the Illumina adapters
`LEADING:3`- This trims off the first 3 bases
`TRAILING:3` - This trims off the last 3 bases
`SLIDINGWINDOW:4:15` - This trims your reads with a sliding window of four bases and an average quality within that window of 15.  If the quality within the four base window drops below 15, it will trim off everything downstream.
`MINLEN:36`- If after trimming the reads are less than 36 bp, they will be discarded.

## SPAdes

No it's time to perform the assembly with your trimmed reads


```{BASH}
spades.py -k 21,51,71,91,111,127 --careful --pe1-1 G11_R1_paired.fq.gz --pe1-2 G11_R2_paired.fq.gz --pe1-s G11_R1_unpaired.fq.gz-o G11_spades_output
```
Here are details on the command
`spades.py` - the command
`-k 21,51,71,91,111,127` - the kmers to use
`--careful` - use careful error correction
`--pe1-1 G11_R1_paired.fq.gz` - the forward paired reads
`--pe1-2 G11_R2_paired.fq.gz` - the reverse paired reads
`--pe1-s G11_R1_unpaired.fq.gz` - the unpaired reads
`-o G11_spades_output` - the output directory to make

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
