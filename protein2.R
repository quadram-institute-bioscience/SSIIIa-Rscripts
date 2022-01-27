
library(data.table)
library(ggplot2)
library(readxl)

proteinDat <- read_excel(path="data/protein_jan2022.xlsx", 
                         sheet="Data Summary",
                         range="A2:F81",
                         col_names = c("ID", "rep", "genotype", "protein","blankCol","average"))

# proteinDat <- read_excel(path="data/ProteinContent201119BHazard.xlsx", 
#                          sheet="Data Summary",
#                          range="A2:E81",
#                          col_names = c("ID", "rep", "genotype", "protein","average"))

setDT(proteinDat)
proteinDat[, genotype := factor(genotype, 
                                levels=c("WT", "A", "B", "D",
                                         "AB", "BD", "AD", "ABD"))  ]
proteinDat[, A := grepl(pattern = "A",genotype)]
proteinDat[, B := grepl(pattern = "B",genotype)]
proteinDat[, D := grepl(pattern = "D",genotype)]

ggplot(proteinDat, aes(genotype, average)) + geom_point()

lmprotein <- lm(data=proteinDat, average ~ genotype)
lmprotein2 <- lm(data=proteinDat, average ~ (A+B+D)^3)

modelsummary(list(lmprotein))
modelsummary(list(lmprotein), statistic = NULL, estimate="{estimate} ({std.error}) p={p.value} {stars}",stars=TRUE, fmt=2)
modelsummary(list(lmprotein2), statistic = NULL, estimate="{estimate} ({std.error}) p={p.value} {stars}",stars=TRUE, fmt=2)

summary(lmprotein)

anova(lmprotein)

emmeans::emmeans(lmprotein, ~genotype)

proteinDat[!is.na(average), 
           .(mean=mean(average), 
             sd=sd(protein)),
           by=genotype]
## No effect of protein 



