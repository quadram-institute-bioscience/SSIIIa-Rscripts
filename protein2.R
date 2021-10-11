
library(data.table)
library(ggplot2)
library(readxl)

proteinDat <- read_excel(path="data/ProteinContent201119BHazard.xlsx", 
                         sheet="Data Summary",
                         range="A2:E81",
                         col_names = c("ID", "rep", "genotype", "protein","average"))


proteinDat[, genotype := factor(genotype, 
                                levels=c("WT", "A", "B", "D",
                                         "AB", "BD", "AD", "ABD"))  ]

setDT(proteinDat)

ggplot(proteinDat, aes(genotype, average)) + geom_point()

lmprotein <- lm(data=proteinDat, average ~ genotype)
anova(lmprotein)

## No effect of protein 