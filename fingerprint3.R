#' ---
#' title: "Analysis of fibre composition"
#' author: "George Savva"
#' date: "November 2021"
#' ---
knitr::opts_chunk$set(fig.width=12, fig.height=8, warning=FALSE, message=FALSE)

library(readxl)
library(data.table)
library(ggplot2)
library(patchwork)
library(tidyr)
library(lmerTest)
library(modelsummary)
library(emmeans)
library(knitr)
library(kableExtra)


#' Loading, cleaning and coding.
#' 
fingerprint2 <- read_excel("data/ssIIIaFingerprintData.xlsx",
                           sheet="Raw data and calculations", 
                           range="B3:AR203")
setDT(fingerprint2)

fingerprint2 <- fingerprint2[!is.na(Run)]
fingerprint2 <- fill(fingerprint2,`Sample No:`, `Sample Name`)
fingerprint2 <- separate(fingerprint2, `Sample Name`, into=c("Sample", "genotype"), sep=" ")
setDT(fingerprint2)
fingerprint2[, Monosubstituted := XA3XX + XA3A3XX + XA3XA3XX + `XA3A2+3XX` + `XA3XA2+3XX`]
fingerprint2[, Disubstituted := `XA2+3XX` + `XA3A2+3XX` + `XA3XA2+3XX`]
fingerprint2[, Substituted := XA3XX + XA3A3XX + XA3XA3XX +`XA2+3XX` + `XA3A2+3XX` + `XA3XA2+3XX`]
fingerprint2[, Unsubstituted := Xylose + Xyl2 + Xyl3 + Xyl5]
fingerprint2[, SumAX := XA3XX + XA3A3XX + XA3XA3XX  + `XA2+3XX` + `XA3A2+3XX` + `XA3XA2+3XX`+                Xylose + Xyl2 + Xyl3 + Xyl5]
fingerprint2[, LogSumAX := log(SumAX)]
fingerprint2[, SumBG := G3 + G4]
fingerprint2[, RatioBG := G3 / G4]
fingerprint2[, LogRatioBG := log(RatioBG)]
fingerprint2[, Monosubstituted := Monosubstituted / `Melibiose IS`]
fingerprint2[, Disubstituted := Disubstituted / `Melibiose IS`]
fingerprint2[, Unsubstituted := Unsubstituted / `Melibiose IS`]
fingerprint2[, Substituted := Substituted / `Melibiose IS`]
fingerprint2[, RatioSubsUnsubs := Substituted / Unsubstituted]
fingerprint2[, LogRatioSubsUnsubs := log(RatioSubsUnsubs)]
fingerprint2[, SumAX := SumAX / `Melibiose IS`]
fingerprint2[, SumBG := SumBG / `Melibiose IS`]


fingerprint2[,A:= genotype %in% c("A", "AB", "AD", "ABD")]
fingerprint2[,B:= genotype %in% c("B", "AB", "BD", "ABD")]
fingerprint2[,D:= genotype %in% c("D", "BD", "AD", "ABD")]
fingerprint2[, genotype := factor(genotype, levels=c("WT", "A", "B", "D", "AB", "AD" ,"BD", "ABD"))]
fingerprint2[, mutations := A+B+D]
fingerprint2[, run :=factor(Run)]


#'
#'  Make a function that estimates models and makes a graph for each outcome.
#'  
doanalysis <- function(class){
  mod2 <- (lmer(data=fingerprint2, get(class) ~ genotype + (1|Sample)))
  mod3way <- (lmer(data=fingerprint2, get(class) ~ (A+B+D)^3 + (1|Sample)))
  mod2way <- (lmer(data=fingerprint2, get(class) ~ (A+B+D)^2 + (1|Sample)))
  modmain <- (lmer(data=fingerprint2, get(class) ~ (A+B+D) + (1|Sample)))
  modnull <- (lmer(data=fingerprint2, get(class) ~ 1 + (1|Sample)))
  modmutants <- (lmer(data=fingerprint2, get(class) ~ factor(mutations) + (1|Sample)))
  ests <- as.data.frame(emmeans(mod2, ~genotype))
  ests$mutations = c(0,1,1,1,2,2,2,3)
  g_sumbg <- ggplot(fingerprint2, aes(genotype, get(class))) + 
    stat_summary(aes(group=Sample), geom="point", pch=1) + 
    geom_errorbar(data=ests, aes(ymin=lower.CL, ymax=upper.CL, y=emmean), width=0.2) + 
    geom_point(data=ests, aes(y=emmean), size=2) + labs(y=class, x="Genotype") + 
    #facet_grid(~mutations, scale="free_x", space="free_x") + 
    theme_bw()
  list(mod2,mod3way, mod2way, modmain,modnull, g_sumbg,ests, modmutants)
  }

#' Now specify the outcomes we want, and generate the models.
#'
listOfOutcomes <- c("SumAX","LogSumAX","Substituted", "Unsubstituted", "SumBG", "LogRatioBG","LogRatioSubsUnsubs" )
mods <- lapply(listOfOutcomes, doanalysis)
names(mods) <- listOfOutcomes
listofOutcomes2 <- c("SumAX", "Substituted", "Unsubstituted", "SumBG")

#' Make a table of ANOVAs for each outcome
#' This tests whether there is any effect of genotype on outcome.  
#' There is an effect of genotype on all outcomes except the BG ratio and the substituted / unsubstituted ratio.
#' 
do.call(rbind,lapply(mods, function(m) anova(m[[1]])))

#' Now for each outcome, check whether the  2-way, 3-way interactions are significant.
#' This shows that the main effects and 2-way interactions are important, but in 
#' general the ABD effect can be explained by adding the lower-order interactions together
#' without any synergy.
#' 
anovas <- lapply(mods, function(m) do.call(anova,m[2:5]))
do.call(rbind, anovas)

#' Check the 3-way vs 1-way for completeness.  
#' This confirms that there are synergies between pairs of mutations.
anovas2 <- lapply(mods, function(m) do.call(anova,m[c(2,4)]))
do.call(rbind, anovas2)


#' Table of model coefficients (compare each genotype vs WT)
#' This is useful if you want to predict the actual effect of specific combinations of 
#' mutations compared to WT.
#' 
#' First all the data (estimates, confidence intervals, p-values, stars) to 3.d.p
#' 
modelsummary(lapply(mods, `[[`,1), 
             estimate  = "{estimate} [{conf.low}, {conf.high}]",
             fmt=3, statistic="p={p.value} {stars}", stars=TRUE)

modelsummary(lapply(mods[listofOutcomes2], `[[`,1), 
             estimate  = "{estimate} ({std.error}) p={p.value} {stars}",
             fmt=2, statistic=NULL, stars=TRUE)
modelsummary(lapply(mods[listofOutcomes2], `[[`,2), 
             estimate  = "{estimate} ({std.error}) p={p.value} {stars}",
             fmt=2, statistic=NULL, stars=TRUE)




#' Now a simpler table 
modelsummary(lapply(mods, `[[`,1), 
             estimate  = "{estimate} [{conf.low}, {conf.high}] {stars}",
             fmt=1, stars=TRUE)

#' Models per mutation number.
modelsummary(lapply(mods, `[[`,8), 
             estimate  = "{estimate} [{conf.low}, {conf.high}] {stars}",
             fmt=1, stars=TRUE)


#' Graph of each outcome.  Error bars correspond to confidence intervals.
wrap_plots( lapply(mods, `[[`, 6))

#' Finally the estimates and confidence intervals for each outcome as plotted on the graph.
ests <- as.data.frame(lapply(lapply(mods, `[[`, 7), function(d) sprintf("%0.2f (%0.2f,%0.2f)", d$emmean, d$lower.CL,d$upper.CL)))
rownames(ests) <- c("WT", "A", "B", "D", "AB", "AD" ,"BD", "ABD")
knitr::kable(ests) |> column_spec (1:6,border_left = T, border_right = T) %>%
  kable_styling()

#' And a table of the ratios (as opposed to the log-ratios above)

ests2 <- as.data.frame(lapply(lapply(mods[5:6], `[[`, 7), function(d) sprintf("%0.2f (%0.2f,%0.2f)", exp(d$emmean), exp(d$lower.CL),exp(d$upper.CL))))
rownames(ests2) <- c("WT", "A", "B", "D", "AB", "AD" ,"BD", "ABD")
knitr::kable(ests2, col.names = c("Ratio BG", "Ratio Sub/Unsubs")) |> column_spec (1:3,border_left = T, border_right = T) %>%
  kable_styling()


