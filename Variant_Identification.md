# Code for mapping and variant identification

## Install packages
```{bash}
conda create -n assembly -c bioconda -c conda-forge sra-tools fastqc=0.11.5 \
             trimmomatic=0.36 spades=3.11.1 quast=5.0.2 \
             bowtie2=2.2.5 prokka java-jdk=8.0.112 samtools bcftools --yes
```
## Activate the environment with all of the programs

```{bash}
conda activate assembly
```

## Quality Filter with Trimmomatic
You need to make sure that you start in the directory with the reads.
You need to make sure that you change the following 
- `Forward_reads.fastq.gz` - should be the R1 file
- `Reverse_reads.fastq.gz` - should be the R2 file
- `Forward_paired.fq.gz` - should be changed to be something like the following `genome_name_forward_paired.fq.gz`
- `Forward_unpaired.fq.gz` - should be changed to be something like the following `genome_name_forward_unpaired.fq.gz`
- `Reverse_paired.fq.gz` - should be changed to be something like the following `genome_name_reverse_paired.fq.gz`
- `Reverse_unpaired.fq.gz`- should be changed to be something like the following `genome_name_reverse_unpaired.fq.gz`

```{BASH} 
java -jar /softwares/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 Forward_reads.fastq.gz Reverse_Reads.fastq.gz Forward_paired.fq.gz Forward_unpaired.fq.gz Reverse_paired.fq.gz Reverse_unpaired.fq.gz ILLUMINACLIP:/softwares/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

## Map the reads to the reference with bowtie2

### Build the reference
```{bash}
bowtie2-build /home/ralhajja/old/Adapted_Genomes/E.coli_K-12_substr_MG1655.fasta /home/ralhajja/old/Adapted_Genomes/E.coli_ref
```

### Map your genome to the reference
You need to change the following
`Forward_paired.fq.gz` should be the forward reads for the strain that you're working with
`Reverse_paired.fq.gz` should be the reverse reads for the strain that you're working with 
`Genome_name` is the name of the genome

```{bash}
bowtie2 -x /home/ralhajja/old/Adapted_Genomes/E.coli_ref -1 genome_name_forward_paired.fq.gz -2 genome_name_reverse_paired.fq.gz -S Genome_name.sam
```
## Variant identification

### Convert `.sam` file to a `.bam` file
```{bash}
samtools view -b Genome_name.sam | samtools sort - > Genome_name.bam
```
### Index the `.bam` file
```{bash}
samtools index Genome_name.bam
```

## Variant calling

### make a new directory for variants
```{BASH}
mkdir variant
cd variant
```
### Copy the reference genome to the `variant` directory
```{BASH}
cp /home/ralhajja/old/Adapted_Genomes/E.coli_K-12_substr_MG1655.fasta .
```
- index the reference for variant calling
```{BASH}
samtools faidx E.coli_K-12_substr_MG1655.fasta
```
- Copy the indexed .bam file to the variant calling dir
```{BASH}
cp ../Genome_name.bam .
```
### Perform variant calling on the `.bam` file
```{BASH}
bcftools mpileup -f E.coli_K-12_substr_MG1655.fasta Genome_name.bam | bcftools call -mv -Ob --ploidy 1 -o genome_name_raw.calls.bcf
```
- Filter the variant calls
```{BASH}
bcftools filter --exclud 'QUAL < 30' genome_name_raw.calls.bcf > genome_name_calls.vcf
```
- Examine all of the variant calls
```{BASH}
less genome_name_calls.vcf
```
- Extract only SNPs
```{BASH}
bcftools view -v snps genome_name_calls.vcf > genome_name_snp.vcf
```
