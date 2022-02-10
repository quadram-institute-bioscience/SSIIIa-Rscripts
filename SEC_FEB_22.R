## SEC data

## restart R and clear workspace before running!

library(readxl)
library(pracma)
library(tidyr)
library(data.table)
library(kableExtra)
library(knitr)
library(patchwork)
library(ggplot2)

#' # Load data
#'

sheets <- excel_sheets(path="data/WheatStarchSistribution8Dec21.xlsx")
SECdataA <- do.call(rbind,lapply(sheets[3:41], FUN=function(x){
       dat <- read_excel(path="data/WheatStarchSistribution8Dec21.xlsx", sheet=x)
       dat <- dat[,c(7,10)]
       names(dat) <- c("wLogVh","length")
       dat$file <- x
       dat
       }))

datesA <- read_excel(path="data/WheatStarchSistribution8Dec21.xlsx",
                     sheet=1,
                     range="AC4:BQ6",col_names = FALSE)
datesA <- as.data.table(t(datesA))
datesA[, date:=as.Date(as.numeric(V2)-2, origin="1900-01-01")]
datesA[, ID:=V1]
datesA <- datesA[, .(ID, date)]

sheets <- excel_sheets(path="data/WheatStarchSistribution9Dec21.xlsx")
SECdataB <- do.call(rbind,lapply(sheets[3:42], FUN=function(x){
        dat <- read_excel(path="data/WheatStarchSistribution9Dec21.xlsx", sheet=x)
        dat <- dat[,c(7,10)]
        names(dat) <- c("wLogVh","length")
        dat$file <- x
        dat
}))


SECdata <- rbind(SECdataA, SECdataB)


#' # Clean, code and normalise.
#' 

setDT(SECdata)
SECdat2 <- SECdata[ length>1 & length<exp(15)]
SECdat2[,totals := sum(wLogVh), file]
SECdat2[,Y := wLogVh / totals]
SECdat2[,X := length]
SECdat2[, c("ID", "REP") := tstrsplit(file, "-", fixed=TRUE)]
SECdat2[REP=="7)", REP:="7"]
SECdat2[, REP:=as.numeric(REP)]
SECdat2[, day:=ifelse(REP%%2==0,2,1)]
SECdat2[, BIOREP:=(REP+1)%/%2]
SECdat2[, c("ID", "GENOTYPE") := tstrsplit(ID, " ", fixed=TRUE)]
SECdat2[, boundaries := cut(length, c(0,37,100,1e20))]
SECdat2[, boundaries2 := cut(length, c(0,30,100,1e20))]
SECdat2[, GENOTYPE := factor(GENOTYPE, levels=c("WT", "A", "B", "D", "AB", "AD", "BD", "ABD"))]
SECdat2[ , A := grepl("A" , GENOTYPE)]
SECdat2[ , B := grepl("B" , GENOTYPE)]
SECdat2[ , D := grepl("D" , GENOTYPE)]
SECdat2[ , mutations := factor(A+B+D)]

SECdat2 <- merge(SECdat2, datesA, all.x=TRUE)

SECdat2[ , order := as.numeric(substr(ID,2,10))]

#' # Plots
#'

ggplot(SECdat2[], aes(y=Y,x=length )) + 
        geom_line(aes(group=REP, col=factor(day))) + 
        facet_grid(GENOTYPE~BIOREP) + scale_x_log10() + 
        scale_color_manual(values=c("red", "black"))


SECmeans <- SECdat2[ , .(meanY=mean(Y)), by=.(GENOTYPE, length, day)]
SECmeans2 <- SECdat2[ , .(meanY=mean(Y)), by=.(mutations, length, day)]

g1 <- ggplot(SECmeans, aes(y=meanY,x=length )) + 
        geom_line(aes(group=interaction(GENOTYPE,day), col=GENOTYPE), size=1) +
        theme_bw() + 
        scale_color_manual(values=c("black", "red", "blue", "green", "orange", "purple", "pink", "brown")) +
        scale_x_log10(labels=c("1","10","37","100","1000","10000","1000000"),
                      breaks=c(1,10,37,100,1000,10000,100000)) + 
        #geom_line(size=.5) + 
        theme_bw() + 
        geom_vline(xintercept=c(37,100), size=c(1,1),col="grey") + 
        theme(legend.position = c(0.8,.5), legend.background = element_rect(linetype = "solid", color="black")) + labs(color="Genotype") + 
        labs(y="w(logVh)",x="Degree of polymerization") 

g2 <-   ggplot(SECmeans2, aes(y=meanY,x=length )) + 
        geom_line(aes(group=interaction(mutations,day), col=mutations), size=1) +
        theme_bw() + 
        facet_wrap(~day)+

        scale_color_manual(values=c("black", "red", "blue", "green")) +
        scale_x_log10(labels=c("1","10","30","37","100","1000","10000","1000000"),
        breaks=c(1,10,30,37,100,1000,10000,100000)) +
        theme_bw() +
#        geom_vline(xintercept=c(37,100), size=c(1,1),col="grey") 
        theme(legend.position = c(0.8,.5), legend.background = element_rect(linetype = "solid", color="black")) + labs(color="Mutations") + 
        labs(y="w(logVh)",x="Degree of polymerization") 


### THIS IS THE ONE FOR PUBLICATION
g3 <-   ggplot(SECmeans2[day==1], aes(y=meanY,x=length )) + 
        geom_line(aes(group=interaction(mutations,day), col=mutations), size=1) +
        theme_bw() + 
        scale_color_manual(
                values=c("black", "red", "blue", "green"),
                labels=c(`0`="Sibling control", 
                         `1`="Single mutants", 
                         `2`="Double mutants",
                         `3`="Triple mutant"),name="Genotype class") +
        scale_x_log10(labels=c("1","10","37","100","1000","10000","1000000"),
                      breaks=c(1,10,37,100,1000,10000,100000)) +
        theme_bw() +
        geom_vline(xintercept=c(37,100), size=c(1,1),col="grey") +
        theme(legend.position = c(0.8,.5), legend.background = element_rect(linetype = "solid", color="black")) + labs(color="Mutations") + 
        labs(y="w(logVh)",x="Degree of polymerization") 
# g1 + ggtitle("Normalised by sum(wlogVh)") 

g2 + ggtitle("Normalised by sum(wlogVh)") + facet_grid(mutations~day)
# g2
ggsave(plot=g3,filename = "SEC.png", width=6, height=4, dpi="retina")
#' # Make summary statistics (peaks, areas, ratios).


summaries <- SECdat2[, .(sums = sum(Y), 
                         maxs = max(Y),
                         areas = sum((Y+shift(Y))/2*(shift(log(X))-log(X)),na.rm = TRUE)
#                         areas2 = trapz(log(X[order(X)]), Y[order(X)]),
#                         areas3 = trapz((X[order(X)]), Y[order(X)])
                         ), 
                     by=.(boundaries, file, GENOTYPE,BIOREP, day,date,order)]

summariesWide <- data.table::dcast(summaries, file + GENOTYPE + BIOREP + day+order   ~ boundaries, 
                       value.var = c("maxs", "sums", "areas"))

summariesWide[ , maxs37vs100 := `maxs_(37,100]`/ (`maxs_(37,100]`+`maxs_(0,37]`)]
summariesWide[ , maxsgl100 := `maxs_(100,1e+20]`/`maxs_(37,100]`]

summariesWide[ , sums37vs100 := `sums_(37,100]`/ (`sums_(37,100]`+`sums_(0,37]`)]
summariesWide[ , sumsgl100 := `sums_(100,1e+20]`/(`sums_(100,1e+20]`+`sums_(0,37]`+ `sums_(37,100]`)]

summariesWide[ , areas37vs100 := `areas_(37,100]`/ (`areas_(37,100]`+`areas_(0,37]`)]
summariesWide[ , areasgl100 := `areas_(100,1e+20]`/(`areas_(100,1e+20]`+`areas_(0,37]`+ `areas_(37,100]`)]



summaries2 <- data.table::melt(summariesWide[, 
                   .(GENOTYPE, file, BIOREP, day, order,
                     maxs37vs100,maxsgl100,
                     sums37vs100,sumsgl100,
                     areas37vs100,areasgl100)],
     id.vars=c("GENOTYPE", "file","day", "BIOREP","order"))

summaries2[ , A := grepl("A" , GENOTYPE)]
summaries2[ , B := grepl("B" , GENOTYPE)]
summaries2[ , D := grepl("D" , GENOTYPE)]
summaries2[ , mutations := factor(A+B+D)]
summaries[ , A := grepl("A" , GENOTYPE)]
summaries[ , B := grepl("B" , GENOTYPE)]
summaries[ , D := grepl("D" , GENOTYPE)]
summaries[ , mutations := factor(A+B+D)]


#' Summary of sum of w(logVh)

propsDay1 <- dcast(summaries[day==1,.(GENOTYPE,sums,boundaries)], 
                   GENOTYPE~boundaries, fun.aggregate = function(x) round(mean(x)*100,digits = 1), value.var="sums")
propsDay2 <- dcast(summaries[day==2,.(GENOTYPE,sums,boundaries)], 
                   GENOTYPE~boundaries, fun.aggregate = function(x) round(mean(x)*100,digits = 1), value.var="sums")

ggplot(summaries, aes(x=GENOTYPE, y=sums,fill=boundaries)) + stat_summary(geom="col", fun = "mean", position=position_stack()) + facet_wrap(~day) + 
        labs(y="Proportion")

ggplot(summaries[day==1], aes(x=factor(date), y=sums,fill=boundaries)) + 
        stat_summary(geom="col", fun = "mean", position=position_stack()) + 
        facet_wrap(~GENOTYPE, ncol=1) + 
        labs(y="Proportion")

ggplot(summaries[day==1], aes(x=factor(order), y=sums,fill=boundaries)) + 
        stat_summary(geom="col", fun = "mean", position=position_stack()) + 
        facet_wrap(~GENOTYPE,  scales = "free_x") + 
        labs(y="Proportion")


names(propsDay1) <- c("Genotype", "Short chain amylopectin", "Long chain amylopectin", "Amylose")
names(propsDay2) <- c("Genotype", "Short chain amylopectin", "Long chain amylopectin", "Amylose")

##' ## Distribution of mass Day 1

kable(propsDay1)

propsDay1

##' ## Distribution of mass Day 2

kable(propsDay2)

summaries[day==1 & boundaries=="(100,1e+20]", 
          sprintf("%0.1f (%0.1f)",
                round(mean(sums)*100,1), 
            round(100*sd(sums),1)),
          by=GENOTYPE]

##' ## All data by genotype

ggplot(summaries2[!(variable %in% c("areasgl100","areas37vs100"))], aes(x=GENOTYPE, y=value,color=factor(day))) + 
        geom_point() + 
        facet_wrap(~variable, scale="free",nrow=2,dir="v", 
                   labeller=labeller(variable=c(`maxs37vs100`="Ratio of peaks\n0-37 vs 37-100",
                                                `maxsgl100`="Ratio of peaks\namylose vs amylopectin",
                                                `sums37vs100`="Proportion long chain amylopectin\nvs total amylopectin",
                                                `sumsgl100`="Proportion amylose"))) + 
        theme_bw() + labs(y=NULL)+
        scale_color_manual(values=c("black", "red"))



##' ## All data by number of mutations

ggplot(summaries2[!(variable %in% c("areasgl100","areas37vs100"))], aes(x=mutations, y=value,color=factor(day))) + 
        geom_point() + 
        facet_wrap(~variable, scale="free",nrow=2,dir="v", 
                   labeller=labeller(variable=c(`maxs37vs100`="Ratio of peaks\n0-37 vs 37-100",
                                                `maxsgl100`="Ratio of peaks\namylose vs amylopectin",
                                                `sums37vs100`="Proportion long chain amylopectin\nvs total amylopectin",
                                                `sumsgl100`="Proportion amylose"))) + 
        theme_bw() + labs(y=NULL)+
        scale_color_manual(values=c("black", "red"))


summaries2[day==1 & variable=="maxsgl100", 
           .(meanpeakratio=sprintf("%0.2f (%0.2f)",mean(value),sd(value))), 
           by=GENOTYPE]

summaries2[day==1, 
           .(meanpeakratio=sprintf("%0.1f (%0.1f)",100*mean(value),100*sd(value))), 
           by=.(mutations, variable)]

summaries2[day==1, 
           .(meanpeakratio=sprintf("%0.1f (%0.1f)",100*mean(value),100*sd(value))), 
           by=.(GENOTYPE, variable)]


##' ## Statistical tests
##' 
lmday1 <- lm(data=summaries2[variable=="sumsgl100" & day==1], value ~ mutations)
lmday2 <- lm(data=summaries2[variable=="sumsgl100" & day==2], value ~ mutations)

lmday1b <- lm(data=summaries2[variable=="sums37vs100" & day==1], value ~ mutations)
lmday2b <- lm(data=summaries2[variable=="sums37vs100" & day==2], value ~ mutations)

lmday1c <- lm(data=summaries2[variable=="maxsgl100" & day==1], value ~ mutations)
summary(lmday1c)

## Testing whether there is an order effect:
lmday1order <- lm(data=summaries2[variable=="sumsgl100" & day==1], value ~ GENOTYPE + order)

## Testing whether there is an order effect:
lmday1table2 <- lm(data=summaries2[variable=="sumsgl100" & day==1], 100*value ~ GENOTYPE)
lmday1table2b <- lm(data=summaries2[variable=="sums37vs100" & day==1], 100*value ~ GENOTYPE)
lmday1table2c <- lm(data=summaries2[variable=="maxsgl100" & day==1], value ~ GENOTYPE)

lmday1table2 <- lm(data=summaries2[variable=="sumsgl100" & day==1], 100*value ~ (A+B+D)^3)
lmday1table2b <- lm(data=summaries2[variable=="sums37vs100" & day==1], 100*value ~ (A+B+D)^3)
lmday1table2c <- lm(data=summaries2[variable=="maxsgl100" & day==1], value ~ (A+B+D)^3)



modelsummary(list(lmday1table2c,lmday1table2, lmday1table2b), 
             estimate  = "{estimate} ({std.error}) p={p.value} {stars}",
             fmt=2, statistic=NULL, stars=TRUE)

summary(lmday1)
summary(lmday1b)

summary(lmday1table2)
summary(lmday1table2b)

summary(lmday1)
library(modelsummary)

modelsummary(models=list("Proportion amylose\nday 1"=lmday1, "Proportion amylose\nday 2"=lmday2, 
                         "Proportion amylopectin\nlong chain\nday 1"=lmday1b, 
                         "Proportion amylopectin\nlong chain\nday 2"=lmday2b), 
             estimate= "{estimate} (p={p.value}) {stars}", statistic = NULL)
