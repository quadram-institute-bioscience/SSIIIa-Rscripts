#' ---
#' title: "Analysis of resistant starch"
#' author: "George Savva"
#' date: "November 2021"
#' ---
#' 

knitr::opts_chunk$set(fig.width=12, fig.height=8, warning=FALSE, message=FALSE)


library(readxl)
library(lmerTest)
library(data.table)
library(modelsummary)
library(patchwork)
library(ggplot2)

dat_rs <- read_excel(path="data/RS_SSIIIa_all_samples2021-11-11.xlsx", sheet=5,
                     range="A1:G81")
names(dat_rs) <- c("Sample", "CrossID", "Genotype", "RSflour","FlourAvg","RSstarch","StarchAvg" )

setDT(dat_rs)
dat_rs[, Genotype:=factor(Genotype, levels=c("WT.", "A", "B", "D", "AB", "AD", "BD", "ABD"))]
dat_rs[,A:= Genotype %in% c("A", "AB", "AD", "ABD")]
dat_rs[,B:= Genotype %in% c("B", "AB", "BD", "ABD")]
dat_rs[,D:= Genotype %in% c("D", "BD", "AD", "ABD")]

dat_rs[,nA:= !Genotype %in% c("A", "AB", "AD", "ABD")]
dat_rs[,nB:= !Genotype %in% c("B", "AB", "BD", "ABD")]
dat_rs[,nD:= !Genotype %in% c("D", "BD", "AD", "ABD")]

dat_rs[, mutations := A+B+D]

g_flour_geno <- ggplot(dat_rs, aes(Genotype, FlourAvg)) + 
  geom_point() + scale_y_log10() +
  theme_bw() + labs(y="RS content of flour (g/100g)")
g_flour_muts <- ggplot(dat_rs, aes(mutations, FlourAvg)) + 
  geom_point() + scale_y_log10() +
  theme_bw() + labs(y="RS content of flour (g/100g)",x="Mutations")
g_flour_geno + g_flour_muts


g_starch_geno <- ggplot(dat_rs, aes(Genotype, StarchAvg)) + 
  geom_point() + scale_y_log10() +
  theme_bw() + labs(y="RS content of starch (g/100g)")
g_starch_muts <- ggplot(dat_rs, aes(mutations, StarchAvg)) + 
  geom_point() + scale_y_log10() +
  theme_bw() + labs(y="RS content of starch (g/100g)",x="Mutations")
g_starch_geno + g_starch_muts



mod1_flour <- lmer(data=dat_rs, RSflour ~ Genotype + (1|CrossID))
mod2_flour <- lmer(data=dat_rs, RSflour ~ factor(mutations) + (1|CrossID))
mod3_flour <-lmer( data=dat_rs, RSflour ~ nA*nB*nD + (1|CrossID))
mod4_flour <-lmer(data=dat_rs, RSflour ~ A*B*D + (1|CrossID))

summary(mod1_flour)


mod1_starch <-lmer(data=dat_rs, RSstarch ~ Genotype + (1|CrossID))
mod2_starch <-lmer(data=dat_rs, RSstarch ~ factor(mutations) + (1|CrossID))
mod3_starch <-lmer(data=dat_rs, RSstarch ~ nA*nB*nD + (1|CrossID))
mod4_starch <-lmer(data=dat_rs, RSstarch ~ A*B*D + (1|CrossID))

summary(mod1_starch)


#' Test each genotype and each number of mutations against the WT.
modelsummary(list("RS in flour"=mod1_flour, "RS in starch"=mod1_starch), statistic = NULL, estimate="{estimate} ({conf.low}-{conf.high}) p={p.value} {stars}",stars=TRUE, fmt=2)
modelsummary(list("RS in flour"=mod2_flour, "RS in starch"=mod2_starch), statistic = NULL, estimate="{estimate} ({conf.low}-{conf.high}) p={p.value} {stars}",stars=TRUE, fmt=2)


modelsummary(list("RS in flour"=mod1_flour, "RS in starch"=mod1_starch), statistic = NULL, estimate="{estimate} ({std.error}) p={p.value} {stars}",stars=TRUE, fmt=2)
modelsummary(list("RS in flour"=mod4_flour, "RS in starch"=mod4_starch), statistic = NULL, estimate="{estimate} ({std.error}) p={p.value} {stars}",stars=TRUE, fmt=2)


#' Now look at the effect of each gene, each 2- and 3-way interaction.  This is difficult to understand.
modelsummary(list("RS in flour"=mod4_flour, "RS in starch"=mod4_starch), statistic = NULL, estimate="{estimate} ({conf.low}-{conf.high}) p={p.value} {stars}",stars=TRUE, fmt=2)

#' Instead we could focus on the triple mutant as the reference, and look at the impact of having WT versions of A, B and D.
#' Here, for example, nATRUE refers to gene A *not* being a mutant.
#' This is probably easier to interpret.  You can see the effect of each gene being intact, and that there is some syngery in that having more than one intact isn't quite the sum of each individually.
modelsummary(list("RS in flour"=mod3_flour, "RS in starch"=mod3_starch), statistic = NULL, estimate="{estimate} ({conf.low}-{conf.high}) p={p.value} {stars}",stars=TRUE, fmt=2)



rs_for_graph <- dat_rs[ !is.na(StarchAvg) , .(Genotype, StarchAvg,FlourAvg) ]


amylose <- melt(read_excel(path="data/data.xlsx", sheet="Amylose (%)", range="A1:H6"))
starch2 <- melt(read_excel(path="data/data.xlsx", sheet="Starch Content (% of flour)", range="A1:H6"))

setDT(amylose)
setDT(starch2)

amylose[ , mutations := 0]
amylose[ variable %in% c("A","B", "D"), mutations := 1]
amylose[ variable %in% c("AB","BD", "AD"), mutations := 2]
amylose[ variable %in% c("ABD"), mutations := 3]

starch2[ , mutations := 0]
starch2[ variable %in% c("A","B", "D"), mutations := 1]
starch2[ variable %in% c("AB","BD", "AD"), mutations := 2]
starch2[ variable %in% c("ABD"), mutations := 3]

rs_for_graph[ , mutations := 0]
rs_for_graph[ Genotype %in% c("A","B", "D"), mutations := 1]
rs_for_graph[ Genotype %in% c("AB","BD", "AD"), mutations := 2]
rs_for_graph[ Genotype %in% c("ABD"), mutations := 3]

  gam <-ggplot(data=amylose , aes(y=value , x=variable)) + 
    stat_summary(geom="bar",aes(fill=factor(mutations)), col="black") + 
    geom_point(color="black")+
    scale_fill_manual(values=c("white","#ffbbbb","#ff7777","#ff2222"))+
    stat_summary(geom="errorbar", width=0.5) + 
    ylab("Amylose (%)") + xlab("ssIIIa") + theme_bw()  + 
    scale_x_discrete(labels=c(`WT`="Sib.Ctrl"))+
    theme(axis.title.x = element_text(face="italic")) + 
    theme(legend.position = "none", axis.text.x = element_text(angle=30, hjust=1))

  gst <- ggplot(data=starch2 , aes(y=value , x=variable)) + 
    stat_summary(geom="bar",aes(fill=factor(mutations)), col="black") + 
    geom_point(color="black")+
    scale_fill_manual(values=c("white","#ffbbbb","#ff7777","#ff2222"))+
    stat_summary(geom="errorbar", width=0.5) + 
    ylab("Starch (%)") + xlab("ssIIIa") + theme_bw()  + 
    scale_x_discrete(labels=c(`WT`="Sib.Ctrl"))+
    theme(axis.title.x = element_text(face="italic")) + 
    theme(legend.position = "none", axis.text.x = element_text(angle=30, hjust=1))
  
  gsf <- ggplot(data=rs_for_graph , aes(y=FlourAvg , x=Genotype)) + 
    stat_summary(geom="bar",aes(fill=factor(mutations)), col="black") + 
    geom_point(color="black")+
    scale_fill_manual(values=c("white","#ffbbbb","#ff7777","#ff2222"))+
    stat_summary(geom="errorbar", width=0.5) + 
    ylab("Resistant Starch (g/100g flour)") + xlab("ssIIIa") + theme_bw()  + 
    scale_x_discrete(labels=c(`WT.`="Sib.Ctrl"))+
    theme(axis.title.x = element_text(face="italic")) + 
    theme(legend.position = "none", axis.text.x = element_text(angle=30, hjust=1))

  library(patchwork)  
  
  gst + gam + gsf + patchwork::plot_annotation(tag_levels = "a")
  
  ggsave("triple.png", width=8, height=5, dpi="retina")
  