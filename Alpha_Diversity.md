## Calculate alpha diversity

### Determine the minimum number of samples to rarify the table.  You shoudl use the full non-rarified table for this.  It's good to trim out samples with less than 1000 sequences first
```{R}
min_lib <- min(sample_sums(Site3_C_D_T_S))
```
### Determine number of samples in your dataset
```{R}
nsamp = nsamples(Site3_C_D_T_S)
```

### Set the number of tables to construct
```{R}
trials = 100
```

### Create matrices for the diferent alpha diversity metrics with sample names to be filled (Here we calaculate richness=Observed, Shannon, and eveness=Inverse simpson
```{R}
richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(Site3_C_D_T_S)

shannon <- matrix(nrow = nsamp, ncol = trials)
row.names(shannon) <- sample_names(Site3_C_D_T_S)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(Site3_C_D_T_S)
```

### It is always important to set a seed when you subsample so your result is replicable 
```{R}
set.seed(81)
```

### Calculate alpha diversity for each metric for each of the 100 tables.
```{R}
for (i in 1:100) {
  # Subsample
  r <- rarefy_even_depth(Site3_C_D_T_S, sample.size = min_lib, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even
  
  # Calculate shannon
  shan <- as.numeric(as.matrix(estimate_richness(r, measures = "Shannon")))
  shannon[ ,i] <- shan
}
```
### Calculate stats for each of the metrics
```{R}
SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

SampleID <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)

SampleID <- row.names(shannon)
mean <- apply(shannon, 1, mean)
sd <- apply(shannon, 1, sd)
measure <- rep("Shannon", nsamp)
Shannon_stats <- data.frame(SampleID, mean, sd, measure)
```

### Combine alpha diversity mean dataframes together
```{R}
alpha <- rbind(rich_stats, even_stats,Shannon_stats)

s <- data.frame(sample_data(Site3_C_D_T_S))
names(s)[names(s) == "X.SampleID"] <- "SampleID"
alphadiv <- cbind(alpha,s)
```

### Separate summarized dataframes for each metrix
```{R}
S3_rarefy_alpha_rich=subset(alphadiv,measure=="Richness")
S3_rarefy_alpha_even=subset(alphadiv,measure=="Inverse Simpson")
S3_rarefy_alpha_shannon=subset(alphadiv,measure=="Shannon")
```

## Plot alpha diversity
### Set color pallete
```{R}
Color_22_Palette <- c("darkgray", "blue", "chocolate", "green","gold","purple","black","tan1","orangered","orange","blue","orangered","magenta","orange","olivedrab3","navy","seagreen","royalblue","red","black","white")
```

### Plot
```{R}
ggplot(S3_rarefy_alpha_shannon, aes(x=Days,y=mean, group=Treatments))+
  geom_line()+
  geom_point(aes(fill=Treatments),shape=21,size=3)+
  scale_fill_manual(values = Color_22_Palette)+
  scale_x_discrete(limits=c("Day0","Day1","Day7","Day14"))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.0001))
```
