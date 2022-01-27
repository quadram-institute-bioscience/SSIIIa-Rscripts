
##
## ALWAYS RESTART R BEFORE RUNNING THIS FILE!!!!!!!!!!!!!!!!!!!

## Intention of this file is to be sourced from the report file.
## No reporting done here, just generation of tests, graphs,
## tables for printing later.

## Load libraries
library(readxl)
library(dplyr)
library(tidyr)
library(lmerTest)
library(ggplot2)
library(emmeans)
library(car)
library(sjPlot)
library(reshape2)
library(patchwork)
library(car)
library(mixdist)
library(broom)
library(broom.mixed)

####
# Each individual dataset is loaded here or another file is sourced.
#### DSC
dsc <- read_excel("copy of DSC data.xlsx",
                  sheet = 1,
                  range = "B4:M34")
dsc <- dsc[, c(1,2,3,9:12)]
dsc<-remove_missing(dsc)
dsc$ID <- dsc$`Plant/Line`
dsc$`Plant/Line` <- NULL

### Create 'starches' table with this information
amylose <- melt(read_excel(path="data/data.xlsx", sheet="Amylose (%)", range="A1:H6"))
starch2 <- melt(read_excel(path="data/data.xlsx", sheet="Starch Content (% of flour)", range="A1:H6"))
starchpergrain <- melt(read_excel(path="data/data.xlsx", sheet="Starch Content (mg per grain)", range="A1:H6"))
RS <- melt(read_excel(path="data.xlsx", sheet="RS", range="A1:H6"))
RSgPer100 <- melt(read_excel(path="data.xlsx", sheet="RS g per 100 g starch", range="A1:H6"))
starches <- cbind(amylose, 
                  starch2, 
                  starchpergrain,
                  RS,
                  RSgPer100)[,c(1,2,4,6,8,10)]

names(starches) <- c("genotype","amylose",
                     "starch","starchPerGrain",
                     "RS","RSgPer100")

starch <- read_excel(path="2018-12-04 SSIIIa Starch Content Mature Grain.xlsx", 
                     range = "a22:am69", sheet="Starch Content")
starch <- select(starch, c(1,2,3,35,37,38))
names(starch) <- c("cross","ID","genotype","flour3","flour","pergrain")
starches$ID <- subset(starch,cross=="2074x0905x0291")$ID
# starches


marvin2 <- read_excel("marvin.xlsx",
                     sheet=2,
                     range="A3:N97")
names(marvin2)[8] <- "Area"
names(marvin2)[9] <- "Width"
names(marvin2)[12] <- "Length"
marvin2 <- subset(marvin2, cross="2074x0905x0291")
## This has all the outcomes in rows, and the genotypes in cols.
## We need it the other way around

germination <-  tidyr::separate(remove_missing(read_excel("germination.xlsx",
                          sheet=1,
                          range="A10:G59")[,c(1,7)]), 1,c("Cross","ID","Genotype"), " ")
germination <- subset(germination, Cross=="2074x0905x0291")
germination$Genotype<-NULL


### FILES BELOW THIS POINT ARE NO LONGER USED, BUT SUPPORT A PREVIOUS VERSION OF THE PAPER.

source("SEC.R")
source("protein.R")
source("mixtures2.R")
source("sugars.R")

allData <- Reduce(  f = function(x,y) merge(x,y,by="ID", all=TRUE),
                    list(marvin2, 
                         germination, 
                         starches, 
                         dsc, 
                         SECsummary, 
                         coulterResults,
                         protein))

allData <- subset(allData, Cross.x=="2074x0905x0291")

# Get block data
blocks <- read_excel("map.xlsx",sheet = "map")
blocks <- reshape2::melt(blocks, id.vars=c("block"))
names(blocks) <- c("block", "position" , "id")
blocks$id <- substr(blocks$id, 4,100)

allData <- merge( allData , blocks , all.x = TRUE , by.x="ID", by.y="id")

# hardness data is hardcoded
allData$hardness <- NA
allData$hardness <- dplyr::recode(allData$ID,
                                  `66-2-126`	=	71.8545472,
                                  `66-7-103`	=	66.25995491,
                                  `66-7-102`	=	51.51449047,
                                  `66-8-39` =	48.78159639,
                                  `66-7-43`	=	59.63320403,
                                  `66-2-95`	=	84.1222995,
                                  `66-11-120`	=	73.549719,
                                  `66-3-37`	=	82.14648296,
                                  `66-4-150`	=	77.81663616,
                                  `66-7-22`	=	77.25538834
                                  )
allData$genotype <- NULL
allData$genotype <- allData$SSIIIa

allData$A <- grepl("A",allData$genotype)>0
allData$B <- grepl("B",allData$genotype)>0
allData$D <- grepl("D",allData$genotype)>0
allData$mutations <- allData$A + allData$B + allData$D
allData$genotype.x<- NULL
allData$genotype.y<- NULL
allData$genotype[allData$genotype=="WT."] <- "WT"
allData$genotype[allData$genotype=="WT"] <- "Sib.Ctrl"
allData$genotype <- factor(allData$genotype, levels=c("Sib.Ctrl","A","B","D","AB","AD","BD","ABD"))
names(allData) <- make.names(names(allData))


# statistical analysis.
# list of outcomes
starchoutcomes <- c("RS","RSgPer100","amylose","starch","starchPerGrain")
dscoutcomes <- c("Enthalpy.J.g...9", "Onset.Temp..C...10","Peak.Temp..C...11","Conclusion.Temp..C...12")
marvinoutcomes <- c("Area", "Width","TGW.g.", "Length")
coulteroutcomes <- c("pi1","mu1","mu2")
secoutcomes <- c("peakratio","arearatio1","arearatio2")
germinationoutcomes <- c("Germination.Index","hardness")
proteinoutcomes <- c("Protein.Dry.basis....NIR")

outcomes <- c(starchoutcomes, dscoutcomes, marvinoutcomes, coulteroutcomes,
              secoutcomes, germinationoutcomes, proteinoutcomes)
outcomeslist <- list(starchoutcomes, dscoutcomes, marvinoutcomes, coulteroutcomes,
                     secoutcomes, germinationoutcomes, proteinoutcomes)


### These are all used for exploring relationships with different outcomes.

modelsFULL <- lapply(outcomes,function(o){
  lmer(data=allData , get(o) ~ genotype + (1|block))
})

modelsLOG <- lapply(outcomes,function(o){
  lmer(data=allData , log(get(o)) ~ genotype + (1|block))
})

modelsMUTATION <- lapply(outcomes,function(o){
  lmer(data=allData , get(o) ~ A*B*D + (1|block))
})

modelsNO_INTERACTION <- lapply(outcomes,function(o){
  lmer(data=allData , get(o) ~ A+B+D + (1|block))
})

modelsLM <- lapply(outcomes,function(o){
  lm(data=allData , get(o) ~ A*B*D)
})

modelsLM2 <- lapply(outcomes,function(o){
  lm(data=allData , get(o) ~ genotype)
})

modelsLM3 <- lapply(outcomes,function(o){
  lm(data=allData , log(get(o)) ~ genotype)
})


#lapply(modelsFULL, summary)
names(modelsFULL) <- outcomes
names(modelsMUTATION) <- outcomes
names(modelsLOG) <- outcomes
names(modelsLM) <- outcomes
names(modelsLM2) <- outcomes
names(modelsLM3) <- outcomes
names(modelsNO_INTERACTION) <- outcomes




outlabs <- outcomes
names(outlabs) <- outcomes
outlabs["RS"] <- "Resistant Starch (g/100g flour)"
outlabs["RSgPer100"] <- "Resistant Starch (g/100g starch)"
outlabs["amylose"] <- "Amylose (%)"
outlabs["starch"] <- "Starch (%)"
outlabs["starchPerGrain"] <- "Starch per grain (mg)"
outlabs["Enthalpy.J.g...9"] <- "Enthalpy (J/g)"
outlabs["Onset.Temp..C...10"] <- "Onset temp (\u00B0C)"
outlabs["Peak.Temp..C...11"] <- "Peak temp (\u00B0C)"
outlabs["Conclusion.Temp..C...12"] <- "Conclusion temp (\u00B0C)"
outlabs["Area"] <- "Grain area"
outlabs["Width"] <- "Grain width"
outlabs["TGW.g."] <- "Thousand grain weight (g)"
outlabs["Length"] <- "Grain length"
outlabs["pi1"] <- "Proportion of volume\nas 'B'-type particles"
outlabs["mu1"] <- "Mean diam. of 'B'-type particles (nm)"
outlabs["mu2"] <- "Mean diam. of 'A'-type particles (nm)"
outlabs["peakratio"] <- "Ratio of SEC peaks"
outlabs["arearatio1"] <- "% amy"
outlabs["arearatio2"] <- "Ratio of SEC areas 2"
outlabs["Germination.Index"] <- "Germination Index"
outlabs["Protein"] <- "Protein (%)"
outlabs["Hardness"] <- "Hardness"

graphs <- lapply(outcomes,function(o){
  ggplot(data=subset(allData,!is.na(get(o))) , aes(y=get(o) , x=genotype)) + 
    stat_summary(geom="bar",fill="red") + 
    geom_point()+
    ylab(o) + theme_bw()
})

graphs2 <- lapply(outcomes,function(o){
  ggplot(data=subset(allData,!is.na(get(o))) , aes(y=get(o) , x=genotype)) + 
    geom_point(color="grey")+
    stat_summary(geom="pointrange") + 
        ylab(o) + theme_bw() + ggtitle(paste0("P-value (lmer)=", round(anova(modelsFULL[[o]])[[6]],5),
                                              ",\nP-value (lm)=", round(anova(modelsLM2[[o]])[1,5],5) ))
})

graphs2b <- lapply(outcomes,function(o){
  ggplot(data=subset(allData,!is.na(get(o))) , aes(y=get(o) , x=genotype)) + 
    geom_point(color="grey")+
    stat_summary(geom="pointrange") + 
    theme_bw() +
    ylab(outlabs[o]) + xlab("SIIIa Genotype Class") 
})

graphs3 <- lapply(outcomes,function(o){
  ggplot(data=subset(allData,!is.na(get(o))) , aes(y=get(o) , x=genotype)) + 
    stat_summary(geom="bar",aes(fill=factor(mutations)), col="black") + 
    geom_point(color="black")+scale_fill_manual(values=c("white","#ffbbbb","#ff7777","#ff2222"))+
    stat_summary(geom="errorbar", width=0.5) + 
    ylab(outlabs[o]) + xlab("ssIIIa") + theme_bw()  + theme(axis.title.x = element_text(face="italic")) + 
    theme(legend.position = "none", axis.text.x = element_text(angle=30, hjust=1))
})

names(graphs2) <- outcomes
names(graphs) <- outcomes
names(graphs3) <- outcomes
names(graphs2b) <- outcomes

## This is only one Brittany is using in the main paper:
 p1 <- patchwork::wrap_plots( list(graphs3[["starch"]] + ggtitle("a"),
                        graphs3[["amylose"]] + ggtitle("b"),
                        graphs3[["RS"]] + ggtitle("c")))
 
 ggsave("multigraph1.png", p1)
 ggsave("multigraph.png", p1, width=6, height=4, dpi=600)
 ggsave("multigraph.eps", p1, width=6, height=4)
 ggsave("multigraph.emf", p1, width=6, height=4)
 

 
## Below here is organising means, sd, and counts into tables. 

means <- lapply(outcomes, function(y){
  d <- aggregate( data = allData , as.matrix(allData[,y]) ~ genotype, FUN=function(x) mean(x, na.rm=TRUE))
  names(d)[2] <- y
  d
})

counts <- lapply(outcomes, function(y){
  d <- aggregate( data = allData , as.matrix(allData[,y]) ~ genotype,FUN=function(x) sum(!is.na(x)))
  names(d)[2] <- y
  d
})

sds <- lapply(outcomes, function(y){
  d <- aggregate( data = allData , as.matrix(allData[,y]) ~ genotype, FUN=function(x) sd(x, na.rm=TRUE))
  names(d)[2] <- y
  d
})

means2 <- t(Reduce( function(x,y) merge(x,y,by="genotype", all=TRUE) , means )[,-1])
sds2 <- t(Reduce( function(x,y) merge(x,y,by="genotype", all=TRUE) , sds )[,-1])
counts2 <- as.data.frame(t(Reduce( function(x,y) merge(x,y,by="genotype", all=TRUE) , counts )[,-1]))

means_and_sds <- data.frame(Map( function(x,y) paste0(signif(x,3)," (", signif(y,2), ")") , as.data.frame(means2) , as.data.frame(sds2)))
names(means_and_sds) <- means[[1]]$genotype[order(as.character(means[[1]]$genotype))]

rownames(counts2) <- outcomes

names(counts2) <- means[[1]]$genotype[order(as.character(means[[1]]$genotype))]
counts2 <- counts2[,as.character(means[[1]]$genotype)]

rownames(means_and_sds) <- outcomes
means_and_sds <- means_and_sds[,as.character(means[[1]]$genotype)]
