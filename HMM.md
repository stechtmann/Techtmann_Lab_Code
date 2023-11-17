# Install hmmer
Create a conda environment for hmmer

```{BASH}
conda create -n HMM -c bioconda hmmer
```

# Activate the HMM environment.

```{BASH}
conda activate HMM
```
# Build an HMM
#### Background on HMMs and Multiple Sequence Alignments.
Hidden Markov Models (HMMs) must start with a multiple sequence alignment (MSA).  An MSA is similar to a pairwise alignment, however it attempts to align multiple sequences to find conserved regions across as set of sequences.  In the case of the use HMMs for protein family identification an MSA is the first step.  The HMM will use this MSA to learn the important sequence regions that can be used for placing a protein into a particular family.  This class we're going to build an HMM to search for PETase genes.  PETases are proteins involved in the depolymerization of PET (plastic #1).  We're also going to build one for alkB (an important gene in oil biodegradation).

#### Set things up

```{BASH}
mkdir HMMs
cd HMMs
```

#### Download an the sequence file for PETase
We're going to use a set of sequences from [Write et al 2021](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01054-5)

```{BASH}
wget https://raw.githubusercontent.com/R-Wright-1/PET-Plastisphere/master/5_RW_thesis/Degradation_Isolates_HMMer/Alignment.fasta
```
#### Download the sequence file for alkB

```{BASH}
wget https://raw.githubusercontent.com/stechtmann/BL2700/master/data/alkB.txt
```

#### Construct a multiple sequence alignment

```{BASH}
muscle -align Alignment.fasta -output PETase.aln
muscle -align alkB.txt -output alkB.aln
```
#### Build the HMM
```{BASH}
hmmbuild PETase.hmm PETase.aln
hmmbuild alkB.hmm alkB.aln
```
# Search set of amino acid sequences

If you have a premade HMM, you can download that and run this through the `hmmsearch` command as well as use the HMMs that we just built.

## Download the amino acid sequence file for one of the bins that we want to annotate.

```{BASH}
wget https://raw.githubusercontent.com/stechtmann/BL4300-5300/master/data/In_class/bin3.faa
```


## Search the `.faa` file with the HMM using `hmmsearch`
```{BASH}
hmmsearch PETase.hmm bin3.faa > PETase.out
hmmsearch alkB.hmm bin3.faa > alkB.out
```

# Examine output.

```{BASH}
less PETase.out
```
There are several hits, but they all have relatively low e-values.

```{BASH}
less alkB.out
```
This bin has two very good hits (e-values of less that 10^-100
