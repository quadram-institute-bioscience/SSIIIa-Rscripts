
library(ggplot2)
library(mixdist)

getMixtures <- function() {

## Mixture modelling coulter counter data.
## This file:
##  Loads data from a list of csv's or all the csv's in one directory.
##  Plots data together
##  Estimates 2-component normal mixture models using mixdist package
##     (though could be easily adapated to other mixtures as mixdist is very flexible)
##  Returns mixture parameters for each file.

## List all the csv files in the directory
## This code only works if working directory is the directory with the data in!
files <- c(dir(".",".CSV", full.names = TRUE), dir(".",".csv", full.names = TRUE))

# Apply to each file:
# load data from each file
cdata <- lapply(files, function(x) {
    #read the csv and name the columns appropriately
    d <- read.csv(file=x, header = FALSE, skip=58, nrows=400)
    names(d) <- c("bin", "diam", "volpc", "numpc", "numpg", "sa")

    #the diam variable is the lower bound of the bin, so make the upper bound as this is 
    #what mixdist needs.
    d$diamupper <- c(d$diam[-1],100)
    
    #get the width of each bin
    d$binwidth <- d$diamupper - d$diam
    
    # Record the filename
    d$file <- x
    
    # return d
    d
    })

# Collapse into a single dataset.
alldata <- do.call(rbind, cdata)


### Now I want to remove the runs where I got a 
library(ggplot2)
# Plot the counts as a line graph.
ggplot( alldata , aes(y=numpc/binwidth, x=diamupper) ) + 
  geom_line() + facet_wrap(~file) + theme_bw()

# Plot the volume contribution as a line graph
# Check for any problems in the data here by visual inspection
ggplot( alldata , aes(y=volpc/binwidth, x=diamupper) ) + 
  geom_line() + facet_wrap(~file) + theme_bw()

### 'mixdist' can estimate mixtures.  it's a better version of 'mix'
#### Read the mixdist manual for how this all works!!

# Set up the initial values based on the estimates from mixtools of the whole dataset.
alldata <- subset(alldata, diamupper<30)
mixpar1 <- data.frame(pi=c(0.2,0.8), mu=c(6,16.6), sigma=c(10,10))

## Don't set constraints on either of the components or the probabability
## You might need to change this if you get computation errors or you know
## that you want to fix certain parameters.
mixcon3 <- mixconstr(conpi="NONE",conmu="NONE", consigma = "NONE")

# Now run mixgroup and mix for each 'file' individually
# Watch for warnings in this step!
mixes3b <- by(alldata, INDICES=alldata$file, FUN=function(d){
  print(d$file[1])
  mix(d[,c("diamupper","volpc")], mixpar1, constr=mixcon3 )
})

# Create a list of density functions from each fit
densityFunctions <- lapply(mixes3b, function(m) function(x) {
  m$parameters$pi[1] * dnorm( x , mean=m$parameters$mu[1], sd=m$parameters$sigma[1]) + 
  m$parameters$pi[2] * dnorm( x , mean=m$parameters$mu[2], sd=m$parameters$sigma[2])
})

# Now use those density functions to create 'modelled' densities in each dataframe
alldataPreds <- do.call(rbind, by( data=alldata, alldata$file  , function(d) {
  d$density <- densityFunctions[[d$file[1]]](d$diamupper)
  d
}))

## Plot the volume data alongside the modelled density
ggplot(alldataPreds , aes(diamupper, density)) + facet_wrap(~file,ncol=10)    + 
  geom_line(aes(y = (volpc/binwidth) / 100), col="blue")+ geom_line(col="red")

## rename to something sensible
means <- as.data.frame(t(sapply(mixes3b, function(x) 
  c(x$parameters$pi,x$parameters$mu,x$parameters$sigma,x$chisq))))
names(means) <- c("pi1","pi2","mu1","mu2","sigma1","sigma2","chi2")

means


#####################################################################
###  ANYTHING FROM HERE DOWN IS SPECIFIC TO BRITTANY'S ANALYSIS   ###
#####################################################################

### IT WON'T WORK WITHOUT THE CROSS FILE

library(car)
library(dplyr)
library(tidyr)

### Remove anything with a chi-2 of more that 3.5.
### I judged this by eye - these were the 7 runs with the obvious errors.
### Brittany suggests these should be removed as obvious technical errors.
means <- subset(means, chi2<3.5)


means$file=rownames(means)
# Remove rownames as it makes the dataset untidy!

rownames(means) <- NULL

# Split filename out to get genotype and IDs
means <- separate(means, file, into=c("ID","genotype","rep"), sep="\\s", remove=FALSE)
means$genotype[means$genotype=="-ABD"] <-"ABD"

alldataPreds$genotype <- sapply(strsplit(alldataPreds$file, split = " ", fixed=TRUE),`[[`,2)
alldataPreds$genotype[alldataPreds$genotype=="-ABD"] <-"ABD"
alldataPreds$genotype[alldataPreds$genotype=="wt"] <-"WT"

meanvolpcs <- aggregate( data = alldataPreds , volpc ~ genotype+diam+binwidth, mean)
meandensity <- aggregate( data = alldataPreds , density ~ genotype+diam, mean)
meansforplotting <- merge(meanvolpcs, meandensity, by=c("diam","genotype"))

gCoulter <- ggplot(subset(meansforplotting, genotype%in%c("ABD","WT")), 
                   aes(x=diam, y= (volpc/binwidth) / 100, color=genotype)) + 
  geom_line() +  geom_line(aes(y=density,x=diam,group=genotype)) + theme_bw() + 
  theme(legend.position = c(0.8,0.8), legend.background = element_rect(linetype = "solid", color="black")) + 
  labs(y="Volume (%)", color="Genotype", x="Particle diameter (\u03BCm)") +
  scale_color_manual(labels=c("Triple mutant","Sibling control"), values=c("red", "blue"))

ggsave(filename = "gCoulter.emf" , gCoulter, width=6, height=4)
ggsave(filename = "gCoulter.eps" , gCoulter, width=6, height=4)
ggsave(filename = "gCoulter.png" , gCoulter, width=6, height=4, dpi=600)


means$file<-NULL
means$rep<-NULL

## Average over technical replicates then do some tidying up.
propA <- aggregate(data=means, .~genotype+ID, FUN=mean)
propA$ID <- substr(propA$ID,3,100)
propA$ID[propA$ID=="69-2-103-"]<- "69-2-103"
list(propA, gCoulter)
}

cr <- getMixtures()

coulterResults <- cr[[1]]
cr[[2]]
# 
# 
# ## Merge this with the file containing the 'cross' information for each ID.
# ## And only keep the ones we are interested in.
# propA <- merge(propA, crosses, all.x = TRUE)
# propA <- subset(propA, cross=="2074x0905x0291")
# 
# ## Now for the proportions and the component locations create a plot and test.
# 
# gpropA <- ggplot(propA, aes(y=pi2, x=genotype)) + stat_summary(geom="bar", aes(fill=mutations)) + geom_point()
# summary(lm(data=propA, pi2 ~ genotype))
# car::Anova(lm(data=propA, pi2 ~ A*B*D))
# 
# gmu1 <- ggplot(propA, aes(y=mu1, x=genotype)) + stat_summary(geom="bar", aes(fill=mutations)) + geom_point()
# summary(lm(data=propA, mu1 ~ genotype))
# car::Anova(lm(data=propA, mu1 ~ A*B*D))
# 
# gmu2 <- ggplot(propA, aes(y=mu2, x=genotype)) + stat_summary(geom="bar", aes(fill=mutations)) + geom_point()
# summary(lm(data=propA, mu2 ~ genotype))
# car::Anova(lm(data=propA, mu2 ~ A*B*D))
# 
# gpropA + theme_bw() | gmu1 + theme_bw() | gmu2 + theme_bw()
