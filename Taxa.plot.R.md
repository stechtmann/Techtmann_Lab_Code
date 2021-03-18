## Read in phyloseq object
If you already of have your data in phyloseq format, you don't need to do this 
```{R}
ps<-phyloseq(otu_table(otus,taxa_are_rows=T),
             sample_data(sampldata),
             tax_table(taxa))
```

## Convert phyloseq object into relative abundance in a melted format.

The purpose this step is to convert your phyloseq object into relative abundance by aglomerating your classes into relative abundance.
To change the rank that you want to use, change the `taxrank` to be another taxonomic level.  This also assumes that you taxa in your taxa table are ranks that correspond to the ranks that are labeled as "Phylum". "Class", "Order", etc. Some outputs have the Ranks labeled as "Rank1", "Rank2", etc.  Make sure to check if you get an error.

This will also melt your table so that it is in an appropriate format for later steps.

```{R}
classabundance <- ps %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()%>%                                          # Melt to long format
  dplyr::arrange(Class) 
```
Check to make sure that this worked by looking at your file.
```{R}
head(classabundance)
```
## Filter out classes with abundance of less than 1%

- Here we want to remove all of the classes that are low abundance (< 1%).    
- You can change the percent cutoff by changing the filte rcommand
- The select function also will pull other metadata categories from your metadata here.  So you should adjust the `select` to be the appropriate metadata categories from your data

```{R}
class <- classabundance %>%
  select(Class, Sample, Biome, Abundance, Oil,Depth) %>%
  group_by(Class, Sample, Biome) %>%
  filter(Abundance > 0.01)
```
Check to make sure that this worked by looking at your file.
```{R}
head(class)
```

## Create an "other"category by extracting all of the classes with abundance of less than 1% to a new object.

Here the all of the claases with low abundance are added together.

```{R}
other<-classabundance %>%
  select(Class, Sample, Biome, Abundance, Oil) %>%
  group_by(Class, Sample, Biome) %>%
  filter(Abundance < 0.01 & Abundance>0)
```

## Add the low abundance classes together and make a new category called "Other"

```{R}
other2<-other%>%
  group_by(SampleID)%>%
  summarise(Abundance = sum(Abundance))%>%
  mutate(Class="Other")%>%
  inner_join(metadata,by="SampleID")%>%
  select(Class, SampleID, Type, Abundance)
```
## Add the "Other" category to the trimmed taxa table

```{R}
full<-bind_rows(class, other2)
```

# Make a plot

## Load a manual color scheme
```{R}
Color_22_Palette <- c("darkgray", "blue", "aquamarine", "chocolate1", "deepskyblue", "gray88", "green","gold","darkred","orangered","magenta","orange","black","cyan","seagreen3","royalblue","red","purple","white","sienna4","darkolivegreen","darkseagreen","mediumvioletred")
```

## Make the plot

```{R}
ggplot(class, aes(x=Sample, y=Abundance, fill= Class))+
  geom_bar(stat="identity",position="stack",color="black")+
  scale_x_discrete(limits=c("St2T5Bakrep3","St2T5Bitrep2","Sup1NSrep1","St2T5Bakrep1","St2T5Ctrlrep2","St2T3Bitrep3","St2T1Ctrlrep2","St2T1Bitrep1","St2T5Bitrep1","St2T1Bakrep2","St2T5Bakrep2","St2T5Ctrlrep1","St1NSrep3","St2T3Bitrep1","St2T3Bitrep2","St2T1Bitrep2","St2T1Bakrep1","St2T1Ctrlrep3","St2T5Ctrlrep3","St2T1Bakrep3","St2T3Ctrlrep1","St2T3Bakrep3","St2T3Bakrep1","St2T1Bitrep3","St2T3Ctrlrep3","CST01401","CST07401","CST03401","CST0540","CST0412001","CST0712001","CST061200","CST0212001","CST0512001","CST081200","CST0632001","CST0432001","CST0532001"))+
  scale_fill_manual(values = Color_22_Palette)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  theme(legend.position = "bottom")+
  theme(legend.text=element_text(size=15))
  ```
