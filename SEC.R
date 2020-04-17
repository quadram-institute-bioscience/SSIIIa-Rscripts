## SEC data

## restart R before running.

library(readxl)
library(pracma)
library(tidyr)

getSEC <- function() {

sheets <- excel_sheets(path="../SEC.xlsx")
SECdata <- do.call(rbind,lapply(sheets[2:11], FUN=function(x){
       dat <- read_excel(path="../SEC.xlsx", sheet=x)
       dat <- dat[,c(7,10)]
       names(dat) <- c("wLogVh","length")
       dat$file <- x
       dat
       }))

SECdata2 <- read_excel(path="../SEC.xlsx", sheet="Combined", range="A4:BT1803", col_names = FALSE)
SECdata2names <- read_excel(path="../SEC.xlsx", sheet="Combined", range="A1:BT1", col_names = FALSE)
id <- unlist(SECdata2names[1,1+6*(0:11)])

names(SECdata2) <- c(outer(c("YA","XA","YB","XB","YT","XT"),id,paste0))
SECdata2$index <- 1:(dim(SECdata2)[1])


wt_ids <- id[1:6]
triple_ids <- id[7:12]

## Make one row per data point
SECdatlong <- pivot_longer(SECdata2, cols=-index, names_to=c("var","rep","ID"), names_pattern="(.)(.)(.*)")

## Make one row per observation (X and Y paired)
SECdatlong2 <- pivot_wider(SECdatlong , id_cols= c("ID","rep","var","index"), names_from=var, values_from=value)

# ggplot( subset(SECdatlong2, X>1 &log(X)<15&  rep=="T") , 
#         aes(x=log(X), y=Y, col=ID)) + 
#         geom_line(size=2) 
# 
# ggplot( subset(SECdatlong2,   rep=="T") , 
#         aes(x=log(X), y=Y, col=ID)) + 
#         geom_line(size=2) 
# 


## Normalise this within the allowed X range:

SECdatlong2 <- subset(SECdatlong2, X>1 &X<exp(15)&  rep=="T")
SECdatlong2$totals <- ave(SECdatlong2$Y, SECdatlong2$ID, FUN=sum)
SECdatlong2$Y2 <- SECdatlong2$Y/SECdatlong2$totals

# ggplot( subset(SECdatlong2, X>1 &log(X)<15&  rep=="T") , 
#          aes(x=log(X), y=Y2, col=ID)) + 
#          geom_line(size=2) 
# 

SECdatlong2$genotype <- ifelse( SECdatlong2$ID %in% wt_ids, "WT","TRIPLE")

# ggplot( subset(SECdatlong2, X>1 &log(X)<15&  rep=="T") ,
#         aes(x=log(X), y=Y, col=genotype, group=ID)) +
#         geom_line(size=2) + facet_wrap(~ID)

means <- aggregate(data=SECdatlong2, Y ~ X+genotype , mean)
sds  <- aggregate(data=SECdatlong2, Y ~ X+genotype , sd)
secmeans <- merge(means, sds, by=c("X","genotype"))
secmeans$uci <- secmeans$Y.x + 1.96*secmeans$Y.y
secmeans$lci <- secmeans$Y.x - 1.96*secmeans$Y.y


(sec_g1 <- ggplot( subset(SECdatlong2, X>1 & X<=100000&  rep=="T") ,
        aes(x=X, y=Y, color=genotype, group=ID)) + scale_x_log10(labels=c("1","10","37","100","1000","10000","1000000"),
                                                                 breaks=c(1,10,37,100,1000,10000,100000)) + 
        #geom_line(size=.5) + 
        theme_bw() + 
        geom_vline(xintercept=c(37,100), size=c(1,1),col="grey") + 
        theme(legend.position = c(0.8,0.8), legend.background = element_rect(linetype = "solid", color="black")) + labs(color="Genotype") + 
        geom_line(data=subset(secmeans, X<100000), aes(group=NULL, y=Y.x), width=1) + 
        labs(y="w(logVh)",x="Degree of polymerization") + 
        scale_color_manual(values=c("red","blue"), labels=c("Triple mutant", "Sibling control")))

ggsave(filename = "SEC.png", sec_g1, width=6, height=4, dpi=600)
ggsave(filename = "SEC.eps", sec_g1, width=6, height=4)
ggsave(filename = "SEC.emf", sec_g1, width=6, height=4)

        

# 
# with(subset(SECdatlong2, ID=="66-7-102"),
#         mixdist::mix(mixdat=data.frame(log(X),Y),
#              mixpar=data.frame(pi=c(0.4,0.1,0.1,0.4),
#                                mu=c(3,4,6,8),
#                                sd=c(.2,.2,.5,.5))))
# 
# short <- subset(SECdatlong2, ID=="66-2-126")
# short <- short[dim(short)[1]:1,]
# mix1 <- mixdist::mix(mixdat=data.frame(log(short$X),short$Y*10000),
#                   mixpar=data.frame(pi=c(0.4,0.3,0.3,0.4),
#                                     mu=c(3,4,6,8),
#                                     sd=c(.2,.1,1,1)))
# 
# plot(mix1)

# plot(log(short$X))

library(pracma)

## Find the heights of the peaks

SECsummary <- do.call(rbind,by(SECdatlong2, SECdatlong2$ID, FUN=function(d){
        secondpeak <- max(d[log(d$X)>5,]$Y)
        firstpeak <- max(d[log(d$X)<5,]$Y)
        peakratio <- firstpeak/secondpeak
        ID <- d$ID[1]
        genotype <- d$genotype[1]
        ## add a cutoff at log(37) and log(100)
        totalarea <- trapz(log(d$X[order(d$X)]), d$Y[order(d$X)])
        area05 <- with(subset(d,log(X)>0 & log(X)<5) , trapz(log(X[order(X)]), Y[order(X)]))
        
        area0to37 <- with(subset(d,log(X)>0 & X<37) , trapz(log(X[order(X)]), Y[order(X)]))
        area37to100 <- with(subset(d,X>=37 & X<100) , trapz(log(X[order(X)]), Y[order(X)]))
        area100plus  <- with(subset(d,X>=100 & log(X)<15) , trapz(log(X[order(X)]), Y[order(X)]))
        
        area515 <- with(subset(d,log(X)>=5 & log(X)<15) , trapz(log(X[order(X)]), Y[order(X)]))
        
        arearatio1 <- 100 * area100plus / (area0to37+area37to100+area100plus)
        arearatio2 <- 100 * area37to100 / (area0to37+area37to100) 
        
        data.frame(ID, genotype,secondpeak, firstpeak, peakratio, totalarea, area05, area515, arearatio1,arearatio2, area0to37, area37to100, area100plus) 
          
        }))
list(SECsummary, sec_g1)
}

gs <- getSEC() 
SECsummary <- gs[[1]]
sec_graph <- gs[[2]]
#sec_graph


### Percentage of the total area that are in each area

#### THESE AREN'T USED 

# ggplot(SECsummary, aes(genotype, arearatio)) + geom_point() 
# ggplot(SECsummary, aes(genotype, peakratio)) + geom_point() 
# ggplot(SECsummary, aes(arearatio, peakratio)) + geom_point(aes(shape=genotype)) 
# 
# t.test(data=SECsummary, arearatio ~ genotype)
# t.test(data=SECsummary, peakratio ~ genotype)
# 
## FInd the areas bigger or smaller than some threshold



