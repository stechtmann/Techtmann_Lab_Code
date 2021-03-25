# Normalize data  
## Calculate minimum number of reads in a sample
```{R}
min(sample_sums(merged_data))
```
## Create an object with all of the sample sums
```{R}
sum=as.data.frame(sample_sums(merged_data))
```
## Delete samples with a mean of less than 1000 (this will depend on your sequencind depth, same number you would use to rarefy in qiime)
```{R}
samplesover1000 = subset_samples(merged_data, sample_sums(merged_data) > 1000)
```
## Check if there are OTUs with no counts
```{R}
any(taxa_sums(samplesover1000) == 0)
```
## If True, how many?
```{R}
sum(taxa_sums(samplesover1000) == 0)
```
## Prune OTUs with no counts 
```{R}
prune_samplesover1000 = prune_taxa(taxa_sums(samplesover1000) > 0, samplesover1000)
```
## Did it work?
```{R}
any(taxa_sums(prune_samplesover1000) == 0)
min(sample_sums(prune_samplesover1000))
```
