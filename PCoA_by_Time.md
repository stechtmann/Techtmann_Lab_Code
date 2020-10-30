## Create an ordination with your phyloseq object 
```{R}
MD.ord1 <- ordinate(straitsr,"PCoA", "bray")
```
`straitsr' is the phyloseq object

## Plot PCoA of all data
```{R}
plotstrr= plot_ordination(straitsr, MD.ord1)#, type = "SITE", color = "TIME", shape= "OIL_TYPE")
```
`SITE`, `TIME`, and `OIL_TYPE` are metadata fields in the metadata file

## Plot the PCoA

### Set color pallet
```{R}
Color_22_Palette <- c("darkgray", "blue", "chocolate", "green","gold","purple","black","tan1","orangered","orange","blue","orangered","magenta","orange","olivedrab3","navy","seagreen","royalblue","red","black","white")
```
### Make plot
```{R}
plotstrr+
  geom_point(aes(fill = OIL_TYPE, shape= SITE),size=5, color="black") + 
  scale_shape_manual(values=c(21,22,23,24,25,26,19))+
  scale_fill_manual(values=Color_22_Palette)+
  theme_set(theme_bw())+
  labs(title = "PCoA Straits ")+
  theme(plot.title = element_text(size=20),
        axis.text.y.left=element_text(size=20),
        axis.text.x=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20))+
  theme(legend.text = element_text(size = 20))+
  guides(fill=guide_legend(override.aes=list(shape=21)))
```
## Make Plot with Axis1 as Y-axis and Time as X-axis

### Extract the vectors from the PCoA ordination
vectors=as.data.frame(MD.ord1$vectors)

### Extract Axis 1
vectors$Axis.1

### Extract the time field from metadata file
```{R}
Time=as.data.frame(straitsr@sam_data$TIME)
names(Time)=c("Time")
```

### Extract other important metadata file that will be needed for plotting 
```{R}
Oil_type=as.data.frame(straitsr@sam_data$OIL_TYPE)
Site=as.data.frame(straitsr@sam_data$SITE)
names(Oil_type)=c("OIL_TYPE")
names(Site)=c("SITE")
```

### Create new metadata file for plotting
```{R}
New=cbind(Time,Oil_type,Site,vectors)
```

### Format new object for plotting
In our data file, time was written as numeric value. I wanted to convert that to a character value that corresponded to timepoint
```{R}
New_Day<-dplyr::mutate(New,Day = ifelse(Time=="T0", 0,
                                    ifelse (Time=="T1", 7,
                                            ifelse(Time=="T2", 14,
                                                   ifelse(Time=="T3", 21,
                                                          ifelse(Time=="T4",28,35))))))
```
### Draw plot
You wil have to adjust the breaks to correspond to what you're calling time points.
```{R}
ggplot(New, aes(x = Time, y = Axis.1)) +
  geom_point(aes(fill = OIL_TYPE, shape = OIL_TYPE), size = 10,color="black") + 
  scale_shape_manual(values=c(21,22,23,24,25,26,19))+
  scale_fill_manual(values=Color_22_Palette)+
  scale_color_manual(values = c("#E96446", "#302F3D", "#87CEFA")) +
  scale_x_discrete(
    breaks = c("T0","T1", "T2", "T3", "T4", "T5"),
    drop = FALSE
  ) +
  facet_grid(~SITE)+
  ylab("PCoA Axis 1 [26%]")+
  theme(plot.title = element_text(size=20),
        axis.text.y.left=element_text(size=20),
        axis.text.x=element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        legend.text=element_text(size=20),legend.title=element_text(size=20))
````
