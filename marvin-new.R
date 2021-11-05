# Marvin Oct 2021.

library(tidyr)
library(readxl)
library(ggplot2)
library(lmerTest)

marvin1 <- read_excel(path="data/marvin.xlsx",
           sheet=2,
           range="A2:N158")

marvin1 <- separate(marvin1, Plant, c("Genotype", "Cross", "Plant"), sep=" ")
marvin1 <- separate(marvin1, Genotype, c("Gene", "Genotype"), sep="-")
setDT(marvin1)
marvin1[ , Genotype := factor(Genotype, levels=c("wt", "abd"))]

names(marvin1)[11] <- "Area"
names(marvin1)[12] <- "Width"
names(marvin1)[15] <- "Length"

marvin2 <- marvin1[`Total grain Weight (g)`>1]

## Remove the outlier?


results <- lapply(c("Main Seeds (total seed number)", 
         "Total grain Weight (g)",
         "Area", "Width", "Length"), function(v){

g1 <- ggplot(data=marvin1, aes(gsub("x","\n",Cross), get(v) )) + 
  geom_point(aes(shape=Genotype, col=Genotype)) + 
             ylab(v) + 
  scale_color_manual(values=c("wt"="black", "abd"="red"))+
  facet_grid(~Gene, scale="free", space="free") + labs(x="Cross")

lmer_all <- lmer(data= marvin2, get(v) ~ Genotype*Gene + (Genotype | Cross))
lm_all   <- lm(  data= marvin2, get(v) ~ Genotype*Cross)
lm_all_noint   <- lm(  data= marvin2, get(v) ~ Genotype+Cross)
lm_all_gene_int   <- lm(  data= marvin2, get(v) ~ Genotype+Gene)
lm_all_nogenotype   <- lm(  data= marvin2, get(v) ~ Cross)

anova(lm_all, lm_all_noint, lm_all_nogenotype)
anova(lm_all, lm_all_gene_int)

c1 <- as.data.frame(emmeans(lmer_all , trt.vs.ctrl ~ Genotype | Gene )$contrast)
c2 <- as.data.frame(emmeans(lm_all , trt.vs.ctrl ~   Genotype | Cross )$contrast)

list(g1,v, c1, c2)
})


results
