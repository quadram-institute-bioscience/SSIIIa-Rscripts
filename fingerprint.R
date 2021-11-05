
library(readxl)
library(data.table)
library(ggplot2)
library(patchwork)
library(tidyr)

fingerprint <- read_excel("data/ssIIIaFingerprintData.xlsx",
                          sheet="Data summary", 
                          range="A1:F41")

fingerprint2 <- read_excel("data/ssIIIaFingerprintData.xlsx",
                          sheet="Raw data and calculations", 
                          range="B3:AR203")


## From Brittany:
# Monosubstituted XA3XX + XA3A3XX + XA3XA3XX + XA3A2+3XX + XA3XA2+3XX

# Disubstituted XA2+3XX + XA3A2+3XX + XA3XA2+3XX

# Unsubstituted Xylose + Xyl2 + Xyl3 + Xyl5

names(fingerprint2)
setDT(fingerprint2)
fingerprint2 <- fingerprint2[!is.na(Run)]
fingerprint2 <- fill(fingerprint2,`Sample No:`, `Sample Name`)
fingerprint2 <- separate(fingerprint2, `Sample Name`, into=c("Sample", "genotype"), sep=" ")
setDT(fingerprint2)

fingerprint2[, Monosubstituted := XA3XX + XA3A3XX + XA3XA3XX + `XA3A2+3XX` + `XA3XA2+3XX`]
fingerprint2[, Disubstituted := `XA2+3XX` + `XA3A2+3XX` + `XA3XA2+3XX`]
fingerprint2[, Unsubstituted := Xylose + Xyl2 + Xyl3 + Xyl5]
fingerprint2[, SumAX := XA3XX + XA3A3XX + XA3XA3XX  + 
                `XA2+3XX` + `XA3A2+3XX` + `XA3XA2+3XX`+
                Xylose + Xyl2 + Xyl3 + Xyl5]
fingerprint2[, SumBG := G3 + G4]


fingerprint2[, Monosubstituted := Monosubstituted / `Melibiose IS`]
fingerprint2[, Disubstituted := Disubstituted / `Melibiose IS`]
fingerprint2[, Unsubstituted := Unsubstituted / `Melibiose IS`]
fingerprint2[, SumAX := SumAX / `Melibiose IS`]
fingerprint2[, SumBG := SumBG / `Melibiose IS`]


table(fingerprint2$Genotype)


fingerprint2[,A:= genotype %in% c("A", "AB", "AD", "ABD")]
fingerprint2[,B:= genotype %in% c("B", "AB", "BD", "ABD")]
fingerprint2[,D:= genotype %in% c("D", "BD", "AD", "ABD")]
fingerprint2[, genotype := factor(genotype, levels=c("WT", "A", "B", "D", "AB", "AD" ,"BD", "ABD"))]
fingerprint2[, mutations := A+B+D]
fingerprint2[, run :=factor(Run)]

library(lmerTest)

ggplot(fingerprint2, aes(genotype, Monosubstituted)) + geom_point()
summary(lmer(data=fingerprint2, Monosubstituted ~ A*B*D + (1|Sample)))
summary(lmer(data=fingerprint2, Monosubstituted ~ genotype + (1|Sample)))

ggplot(fingerprint2, aes(genotype, Disubstituted)) + geom_point()
summary(lmer(data=fingerprint2, Disubstituted ~ A*B*D + (1|Sample)))
summary(lmer(data=fingerprint2, Disubstituted ~ genotype + (1|Sample)))

ggplot(fingerprint2, aes(genotype, Unsubstituted)) + geom_point()
summary(lmer(data=fingerprint2, Unsubstituted ~ A*B*D + (1|Sample)))
summary(lmer(data=fingerprint2, Unsubstituted ~ genotype + (1|Sample)))

ggplot(fingerprint2, aes(genotype, SumBG)) + geom_point()
summary(lmer(data=fingerprint2, SumBG ~ A*B*D + (1|Sample)))
summary(lmer(data=fingerprint2, SumBG ~ genotype + (1|Sample)))

ggplot(fingerprint2, aes(genotype, SumAX)) + geom_point()
summary(lmer(data=fingerprint2, SumAX ~ A*B*D + (1|Sample)))
summary(lmer(data=fingerprint2, SumAX ~ genotype + (1|Sample)))


## Check the runs.
ggplot(fingerprint2, aes(x = run, y=`Melibiose IS`)) + 
  geom_boxplot() + geom_point()


## OLD ANALYSIS
setDT(fingerprint)
names(fingerprint) <- c("sampleNo", "sample", "genotype","run",
                        "sumAX", "sumBG")

fingerprint[,A:= genotype %in% c("A", "AB", "AD", "ABD")]
fingerprint[,B:= genotype %in% c("B", "AB", "BD", "ABD")]
fingerprint[,D:= genotype %in% c("D", "BD", "AD", "ABD")]
fingerprint[, genotype := factor(genotype, levels=c("WT", "A", "B", "D", "AB", "AD" ,"BD", "ABD"))]
fingerprint[, mutations := A+B+D]
fingerprint[, run :=factor(run)]
gax <- ggplot(fingerprint, aes(genotype, sumAX)) + geom_point(aes(shape=run)) + facet_grid(~mutations, space="free", scales="free")+ 
  labs(x="Genotype", y="sum ax")
gbg <- ggplot(fingerprint, aes(genotype, sumBG)) + geom_point(aes(shape=run)) + facet_grid(~mutations, space="free", scales="free")+ 
  labs(x="Genotype", y="sum bg")

ggplot(fingerprint, aes(sumAX, sumBG)) + geom_point(aes(shape=run)) + facet_grid(~mutations, space="free", scales="free")+ 
  labs(x="Genotype", y="sum bg")

gax + gbg + plot_layout(guides="collect")

mutations <- sum()
lmax1 <- lm(data=fingerprint , sumAX ~ A*B*D )
lmbg1 <- lm(data=fingerprint , sumBG ~ A*B*D )

lmax_2nd_order <- lm(data=fingerprint , sumAX ~ (A+B+D)^2 )
lmbg_2nd_order <- lm(data=fingerprint , sumBG ~ (A+B+D)^2 )

lmax_1st_order <- lm(data=fingerprint , sumAX ~ (A+B+D) )
lmbg_1st_order <- lm(data=fingerprint , sumBG ~ (A+B+D) )

lmax_null <- lm(data=fingerprint , sumAX ~ 1 )
lmbg_null <- lm(data=fingerprint , sumBG ~ 1 )

lmax2 <- lm(data=fingerprint , sumAX ~ genotype )
lmbg2 <- lm(data=fingerprint , sumBG ~ genotype )


summary(lmax1)
summary(lmbg1)

summary(lmax_2nd_order)
summary(lmbg_2nd_order)

summary(lmax2)
summary(lmbg2)

### Type 1 ANOVA reveals the order in which the main effects and interactions matter
anova(lmax1)
anova(lmbg1)

summary(lmax_1st_order)
summary(lmbg_1st_order)

summary(lmax_2nd_order)
summary(lmbg_2nd_order)
## For both outcomes it is clear that each of the mutations is important.
## B and D are consistently more important than A.

lmax_muts <- lm(data=fingerprint , sumAX ~ factor(mutations) )
lmbg_muts <- lm(data=fingerprint , sumBG ~ factor(mutations) )

summary(lmax_muts)
summary(lmbg_muts)

emmeans(lmax_muts, consec ~ mutations , adjust="none")
emmeans(lmbg_muts, consec ~ mutations , adjust="none")

anova(lmax_null, lmax_1st_order, lmax_2nd_order, lmax1)
anova(lmbg_null, lmbg_1st_order, lmbg_2nd_order, lmbg1)

## There isn't much correlation between amylose and RS beyond this, when taken at the RS level.




