## Install packages for dunn.test
```{R}
install.packages("dunn.test")
```
## Activate dunn.test package
```{R}
library(dunn.test)
```
## Subset data to only include the data that you need
```{R}
qPCR_all_LG<-subset(qPCR3, Site =="LG")
```
## Use the shapiro wilk test for normality to determine if the data is normally distributed.
This test looks at the data that you are using in this case gene copies.
```{R}
shapiro.test(qPCR_all_LG$Genecopies)
```
If the p-value is significant it means that the data is NOT normally distrubuted and you cannot use an ANOVA
## Perform the Kruskal Wallis Test
This is test compares a single category at a time.

### Compare by nanoparticle and biocides
```{R}
kruskl.test(Genecopies~Nano_bio,data=qPCR_all_LG)
```
### Compare between treatments
```{R}
kruskal.test(Genecopies~Treatments,data=qPCR_all_LG)
````

### Compare between time
```{R}
kruskal.test(Genecopies~Time,data=qPCR_all_LG)
{R}
```
None of these were significant so we want to see if there is a significant impact of treatment by time.

### Create a new column in your file that is nanoparticles and time (e.g. Nanoparticles time 0)
```{R}
qPCR_all_LG_time<-qPCR_all_LG%>%
  unite(Nano_bio_time,Time, Nano_bio,sep="_")
```
### Perform a kruskal wallis test with Nano/Bio and time
```{R}
kruskal.test(Genecopies~Nano_bio_time,data=qPCR_all_LG_time)
```

This test was significant.  We will then use the dunn test as a post hoc test to see between which categories is the test significant.

### Perform the Dunn test
```{R}
dunn.test(x=qPCR_all_LG_time$Genecopies,g=qPCR_all_LG_time$Nano_bio_time,method="bh")
```
Here we corrected for multiple comparisons using the Benjamini-Hochberg correction for multiple comparisons (`bh`)
