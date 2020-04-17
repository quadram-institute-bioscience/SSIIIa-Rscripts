
library(readxl)
library(ggplot2)

sugars <- read_excel(path = "sugars.xlsx", sheet=2)
sugars$Genotype <- factor(sugars$Genotype, level=c("Control", "A", "B", "D", "AB", "AD", "BD", "ABD"))

sugars$A <- grepl("A",sugars$Genotype)>0
sugars$B <- grepl("B",sugars$Genotype)>0
sugars$D <- grepl("D",sugars$Genotype)>0
sugars$mutations <- factor(sugars$A + sugars$B +sugars$D)

sugars$Sugar <- factor(sugars$Sugar, level=c("Sucrose", "Glucose", "Fructose"))
g_sugars <- ggplot(remove_missing(sugars), aes(Genotype, Level)) + 
  geom_point(col="grey") + 
  facet_wrap(~Sugar, scales="free") + 
  theme_bw() + scale_y_log10() + 
  stat_summary(fun.y = "mean", geom="point") +
  stat_summary(geom="errorbar")

g_sugars

# aes(fill=factor(mutations)), col="black") + 
#   geom_point(color="black")+scale_fill_manual(values=c("white","#ffbbbb","#ff7777","#ff2222"))

g_sugars_bar <- ggplot(remove_missing(sugars), aes(Genotype, Level)) + 
  facet_wrap(~Sugar, scales="free") + 
  theme_bw() + 
  stat_summary(fun.y = "mean", geom="bar", aes(fill=factor(mutations)), col="black") +
  stat_summary(geom="errorbar", width=0.5) +
  geom_point(col="black") + scale_fill_manual(values=c("white","#ffbbbb","#ff7777","#ff2222")) + 
  labs(y="Concentration (\u03bcl per g)") + theme(legend.position = "none", axis.text.x = element_text(angle=30, hjust=1))

g_sugars_bar
ggsave(filename = "sugars.png", g_sugars_bar, width=8, height=6)
library(emmeans)

lm_sucrose <- lm(data= subset(sugars, Sugar=="Sucrose"), log(Level)~Genotype)
lm_glucose <- lm(data= subset(sugars, Sugar=="Glucose"), log(Level)~Genotype)
lm_fructose <- lm(data= subset(sugars, Sugar=="Fructose"), log(Level)~Genotype)

lm_sucrose_NOLOG <- lm(data= subset(sugars, Sugar=="Sucrose"), Level~  Genotype)
lm_glucose_NOLOG <- lm(data= subset(sugars, Sugar=="Glucose"), Level ~  Genotype)
lm_fructose_NOLOG <- lm(data= subset(sugars, Sugar=="Fructose"), Level~ Genotype)


lm_sucrose_LINEAR <- lm(data= subset(sugars, Sugar=="Sucrose"), log(Level)~  A+B+D)
lm_glucose_LINEAR <- lm(data= subset(sugars, Sugar=="Glucose"), log(Level)~  A+B+D)
lm_fructose_LINEAR <- lm(data= subset(sugars, Sugar=="Fructose"), log(Level)~A+B+D)


lm_sucrose_MUT <- lm(data= subset(sugars, Sugar=="Sucrose"), log(Level)~  A*B*D)
lm_glucose_MUT <- lm(data= subset(sugars, Sugar=="Glucose"), log(Level)~  A*B*D)
lm_fructose_MUT <- lm(data= subset(sugars, Sugar=="Fructose"), log(Level)~A*B*D)

lm_sucrose_MUT_NOLOG <- lm(data= subset(sugars, Sugar=="Sucrose"), Level~  A*B*D)
lm_glucose_MUT_NOLOG <- lm(data= subset(sugars, Sugar=="Glucose"), Level~  A*B*D)
lm_fructose_MUT_NOLOG <- lm(data= subset(sugars, Sugar=="Fructose"), Level ~A*B*D)

lm_fructose_MUT

lm_sucrose_N <- lm(data= subset(sugars, Sugar=="Sucrose"), (Level)~  A*B*D - A:D)
lm_glucose_N <- lm(data= subset(sugars, Sugar=="Glucose"), (Level)~  A*B*D - A:D)
lm_fructose_N <- lm(data= subset(sugars, Sugar=="Fructose"), (Level)~A*B*D - A:D)



sjPlot::tab_model(lm_sucrose, lm_glucose, lm_fructose, dv.labels = c("Sucrose", "Glucose", "Fructose"), transform = "exp")
sjPlot::tab_model(lm_sucrose_MUT, lm_glucose_MUT, lm_fructose_MUT, dv.labels = c("Sucrose", "Glucose", "Fructose"), transform = "exp")
sjPlot::tab_model(lm_sucrose_LINEAR, lm_glucose_LINEAR, lm_fructose_LINEAR, dv.labels = c("Sucrose", "Glucose", "Fructose"), transform = "exp")

sjPlot::tab_model(lm_sucrose_NOLOG, lm_glucose_NOLOG, lm_fructose_NOLOG, dv.labels = c("Sucrose", "Glucose", "Fructose"))
sjPlot::tab_model(lm_sucrose_MUT_NOLOG, lm_glucose_MUT_NOLOG, lm_fructose_MUT_NOLOG, dv.labels = c("Sucrose", "Glucose", "Fructose"))

sjPlot::tab_model(lm_sucrose_N, lm_glucose_N, lm_fructose_N, dv.labels = c("Sucrose", "Glucose", "Fructose"))


anova(lm_sucrose_MUT, lm_sucrose_LINEAR)
anova(lm_glucose_MUT, lm_glucose_LINEAR)
anova(lm_fructose_MUT, lm_fructose_LINEAR)

summary(lm_sucrose_MUT)
summary(lm_glucose_MUT)
summary(lm_fructose_MUT)

summary(lm_sucrose_LINEAR)
summary(lm_glucose_LINEAR)
summary(lm_fructose_LINEAR)


head(sugars)

sugarstab <- aggregate( data=sugars , Level ~ Sugar + Genotype , 
                        FUN=function(x) sprintf("%0.2f (%0.2f)",mean(x) , sd(x))
                        )
sugarstab <- reshape(data=sugarstab,v.names = "Level", timevar = "Sugar", idvar = "Genotype",direction="wide", )
names(sugarstab)<- c("Genotype", "Sucrose", "Glucose", "Fructose")

