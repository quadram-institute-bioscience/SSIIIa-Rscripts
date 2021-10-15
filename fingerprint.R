
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
