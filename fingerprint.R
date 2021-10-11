
library(readxl)
library(data.table)
library(ggplot2)
library(patchwork)

fingerprint <- read_excel("data/ssIIIaFingerprintData.xlsx",
                          sheet="Data summary", 
                          range="A1:F41")


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


## For BG the 'significant' effects are the A, B and D main effects, the A/B interactions and possibly the A/B/D interaction
## For AX the 'significant' effects are the A, B and D main effects, the B/D interaction possibly the A/B and A/D interactions


