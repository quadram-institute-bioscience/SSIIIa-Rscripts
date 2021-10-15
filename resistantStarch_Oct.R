# rs_oct_2021

library(readxl)
library(data.table)
library(ggplot2)
library(tidyr)

#### QUestions
## Estimate difference between ABD and WT for each gene.
## Comparing the different triple mutants to each other.
###### Ranking the effects of the genes.

###### Do the WT sibs vary across the sample?

datRS <- read_excel("data/2021-09-08RSonstarch.xlsx",
                         sheet="Calc",
                         range="A9:J165")

names(datRS) <- c("crossID", "SIIIagenotype", 
                  "flourweight", "RSdigestvolume", "volumeassayed",
                  "repA", "repB", "repAvg", "RSper100g", "RSaverage")



## How is amylose linked to resistant starch.
## attributable fraction amylose/RS

## Independent mutants for each gene - can we estimate the differences between them?

setDT(datRS)
datRS <- separate(datRS, crossID, c("Genotype", "Cross1", "Cross2"),sep = " ")
datRS <- separate(datRS, Genotype, c("strain", "mutant"),sep = "-")
setDT(datRS)
datRS[, CrossID := as.numeric(factor(Cross1)), by=.(strain, mutant)]
View(datRS)

gRS <- ggplot(datRS,aes(mutant, RSper100g)) + 
  geom_point(aes(shape=factor(CrossID))) + 
  coord_flip() + scale_y_log10() + 
  facet_grid(rows=vars(strain)) + 
  theme_bw()
gRS


