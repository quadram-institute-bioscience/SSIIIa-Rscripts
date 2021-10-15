
library(readxl)
library(data.table)
library(ggplot2)
library(tidyr)

datAmylose <- read_excel("data/2021-09-23AmyloseContent.xlsx",
                         sheet="Gene Average",
                         range="A5:C161")

names(datAmylose) <- c("plant", "id", "amylose")


datAmylose <- separate(datAmylose, plant, c("Genotype", "Cross1", "Cross2"),sep = " ")
datAmylose <- separate(datAmylose, Genotype, c("strain", "mutant"),sep = "-")
setDT(datAmylose)
datAmylose[, CrossID := as.numeric(factor(Cross1)), by=.(strain, mutant)]
View(datAmylose)

(gAMy <- ggplot(datAmylose,aes(mutant, amylose)) + 
  geom_point(aes(shape=factor(CrossID), group=factor(CrossID)), position=position_dodge(width=.3)) + 
    coord_flip() + scale_y_log10() + 
  facet_grid(rows=vars(strain)) + 
  theme_bw())

gRS + gAMy + plot_layout(guides="collect")

datMerged <- merge(datAmylose, datRS)
datMergedAverage <- datMerged[, .(meanAmylose=mean(amylose), meanRS=mean(RSper100g)),by=.(mutant, CrossID, strain)
                              ]
  
## There is no correlation between amylose and starch at the individual sample level.
ggplot(datMerged, aes(amylose, RSper100g)) + geom_point(aes(shape=mutant, color=mutant)) + 
  scale_shape_manual(values=1:10) + scale_x_log10() + scale_y_log10() + theme_bw() + 
  scale_color_manual(values=c("red", "black"))+ facet_wrap(~paste(strain,Cross1))
## Within plant there is no evidence that amylose predicts resistant starch.
summary(lm(data=datMerged , log(RSper100g) ~ log(amylose) + mutant*Cross1))

## Similarly among the wild types.. but there is no variation in amylose at the WT level so..

library(ggrepel)
ggplot(datMergedAverage, aes(meanAmylose, meanRS)) + 
  geom_point(aes(color=mutant), size=2) + 
  scale_shape_manual(values=10:19) + scale_x_log10() + scale_y_log10() + theme_bw() + 
  scale_color_manual(values=c("red", "black")) +
  geom_line(aes(group=paste(strain,CrossID)))

ggplot(datMergedAverage, aes(meanAmylose, meanRS)) + 
  geom_point(aes(color=mutant)) + 
  scale_shape_manual(values=10:19) + scale_x_log10() + scale_y_log10() + theme_bw() + 
  scale_color_manual(values=c("red", "black")) + 
  geom_label_repel(aes(label=paste(strain,mutant)))


## Does the difference in amylose caused by the treatment predict the difference in starch?

datMerged$amyloseAVG <- ave(datMerged$amylose, datMerged$CrossID, datMerged$strain, FUN=mean)

summary(lmer(dat=datMerged[mutant=="abd"], log(RSper100g) ~ log(amylose) + (1|strain:CrossID)))
summary(lmer(dat=datMerged[mutant=="abd"], log(RSper100g) ~ log(amyloseAVG) + (1|strain:CrossID)))
summary(lm(dat=datMergedAverage[mutant=="abd"], log(meanRS) ~ log(meanAmylose) ))




##### WT ONLY.

summary(lmer( data=datMerged[mutant=="wt"], log(amylose) ~ (1|strain/CrossID) ))
summary(lmer( data=datMerged[mutant=="wt"], amylose ~ (1|strain/CrossID) ))

## Amylose does not vary across WT strains
summary(lm( data=datMerged[mutant=="wt"], log(amylose) ~ paste0(strain,CrossID) ))

## RS does vary across WT strains
summary(lm( data=datMerged[mutant=="wt"], log(RSper100g) ~ paste0(strain,CrossID) ))


ggplot(datRS[mutant=="wt"], aes(gsub(pattern = "x","\n",Cross1), RSper100g)) + geom_point() + 
  facet_grid(cols=vars(strain), space="free", scale="free") + scale_y_log10() + 
  labs(x="Strain", y="Resistant starch") + theme_bw()

ggplot(datMerged[mutant=="wt"], aes(gsub(pattern = "x","\n",Cross1), amylose)) + geom_point() + 
  facet_grid(cols=vars(strain), space="free", scale="free") + scale_y_log10() + 
  labs(x="Strain", y="Resistant starch") + theme_bw()

## There is no evidence that amylose varies across WT strains, 
## but a lot of evidence that resistant starch varies across WT strains.

####  Is there a correlation between the increase in amylose and the increase in RS?

dmalong <- melt(datMergedAverage, id.vars=c("strain", "mutant", "CrossID"))

dmawide <- dcast(dmalong, CrossID+strain ~ mutant + variable)

ggplot(dmawide , aes(wt_meanAmylose, abd_meanAmylose)) + geom_point() + 
  coord_fixed() + geom_abline(slope=1, intercept=0)+ 
  scale_y_log10() + scale_x_log10()

ggplot(dmawide , aes(wt_meanRS, abd_meanRS)) + geom_point() + 
  coord_fixed() + geom_abline(slope=1, intercept=0) + 
  scale_y_log10() + scale_x_log10()

summary(lm(data=dmawide , log(abd_meanRS) ~ log(abd_meanAmylose) + log(wt_meanRS) + log(wt_meanAmylose)))

dmawide
dmawide2 <- dmawide[, .(abd_meanAmylose=mean(abd_meanAmylose),
                        wt_meanAmylose=mean(wt_meanAmylose),
                        abd_meanRS=mean(abd_meanRS),
                        wt_meanRS=mean(wt_meanRS)), by=strain]
dmawide2

summary(lm(data=dmawide2 , log(abd_meanRS) ~ log(abd_meanAmylose) + log(wt_meanRS) + log(wt_meanAmylose)))
summary(lm(data=dmawide2 , log(abd_meanRS) ~ log(wt_meanRS) + log(abd_meanAmylose) ))
summary(lm(data=dmawide2 , log(abd_meanRS) ~ log(abd_meanAmylose) ))
dmawide2$difRS=dmawide2$abd_meanRS - dmawide2$wt_meanRS
dmawide2$difAmy=dmawide2$abd_meanAmylose - dmawide2$wt_meanAmylose

## For every increase in log amylose there is an increase in log RS.
summary(lm(data=dmawide2 , log(difRS) ~ log(abd_meanAmylose) ))

## If we try to correlate the average change in RS with the change in amylose at the strain level..
summary(lm(data=dmawide2 , log(difRS) ~ log(difAmy) ))
summary(lm(data=dmawide2 , difRS ~ difAmy ))

ggplot(dmawide2, aes(difRS, difAmy)) + geom_point()
