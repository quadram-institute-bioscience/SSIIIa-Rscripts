
library(readxl)
library(data.table)
library(ggplot2)
library(patchwork)
library(tidyr)

fingerprint2 <- read_excel("data/ssIIIaFingerprintData.xlsx",
                           sheet="Raw data and calculations", 
                           range="B3:AR203")


## From Brittany:
# Monosubstituted XA3XX + XA3A3XX + XA3XA3XX + XA3A2+3XX + XA3XA2+3XX

# Disubstituted XA2+3XX + XA3A2+3XX + XA3XA2+3XX

# Unsubstituted Xylose + Xyl2 + Xyl3 + Xyl5

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


fingerprint2[,A:= genotype %in% c("A", "AB", "AD", "ABD")]
fingerprint2[,B:= genotype %in% c("B", "AB", "BD", "ABD")]
fingerprint2[,D:= genotype %in% c("D", "BD", "AD", "ABD")]
fingerprint2[, genotype := factor(genotype, levels=c("WT", "A", "B", "D", "AB", "AD" ,"BD", "ABD"))]
fingerprint2[, mutations := A+B+D]
fingerprint2[, run :=factor(Run)]


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

