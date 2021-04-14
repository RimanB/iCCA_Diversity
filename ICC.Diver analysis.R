library(tidyverse)
library("survminer")
library(haven)
library(ggcorrplot)
library(pander)
library(psych)
library("corrgram")


dat <- read_sav("/Volumes/Macintosh.D/ICC.paper/ICC with stem cell feature/OS.ICC.sav")
dat2 <- read_sav("/Volumes/Macintosh.D/ICC.paper/ICC with stem cell feature/ICC.Final analysis.sav")
head(dat2)
attach(dat2)


vars <- cbind (SBD, CLC, DPM, Rep.HGP, Size.cm, MVD.mm, Portal.tract, AVD.mm, u.AVD.mm, OS.m, DFS.m)
a1 <- pairs.panels(vars, stars = TRUE,col.stars = "red", cex.cor = 0.8, cor=TRUE, show.legend=TRUE, colors = TRUE, adjust="BH", method = "pearson", alpha = .05)

cocol<- colorRampPalette(c("#6D9EC1", "white", "#E46726"))

corPlot(vars, stars=TRUE, cex =1, keep.par=FALSE, upper = TRUE, scale = TRUE, pval = TRUE, adjust="BH")



mycorrelations <- psych::corr.test(vars, use = "pairwise",method="pearson",adjust="BH",
          alpha=.05,ci=TRUE,minlength=5)

mycorrelations$r %>%
  pander()

mycorrelations$p %>%
  pander()
write.csv(mycorrelations$p, file="iCCA_corre-p-values.csv")

library("PerformanceAnalytics")
chart.Correlation(vars, histogram=TRUE, pch=19, col = "blue")

library(RColorBrewer)
library(relaimpo)

model1 = lm(OS.m ~ BD + Pus.HGP + Size.cm + MVD.mm + AVD.mm + Portal.tract)

metrics = calc.relimp (model1, type = c("lmg"))
metrics
attach(dat)

# * p<0.05, ** p<0.01, *** p<0.001 

# Diagnosis5a in survirval analysis = Diagnosis5WHO in others

# Definition: HCC-like = 100 %  hepatoid pattern
# Bile ductular (cholangiolocarcinoma >80% CLC pattern = WHO 2019)
# Small bile duct > SBD component < 5 %, in this study actually SBD = 100 % SBD pattern 
# Mixied with DPM => DPM pattern >=5% 
# Mixed SBD-BD: 5% =< BD =< 80% 
# DFS: Disease free survival: The period from operation to first recurence or die of disease with either comimng first
# Cases no recurence and alive at last followup were censored for DFS 
# OS: Overal survival: Time from operation till death regardless of etiology 
# Cases alive at last followup were censored for OS 
# Recur = original data 
# Recur2 = Original data + death of disease - stage 4 cases or case die within 1 month
# Recur 3 = Original data + death of disease including stage 4
# reported Recur2 in paper. 
# Cases with stage IV were excluded from DFS analysis
# Cases died within one months from operation were excluded from OS analysis 
# Interpretation of the results: 



os5 <- survfit(Surv(OSY, Death) ~ Diagnosis5, data = dat)
os5

os5a <- survfit(Surv(OSY, Death) ~ Diagnosis5a, data = dat)
os5a

os.SOTH <- survfit(Surv(OSY, Death) ~ Diagnosis.SOTH, data = dat)
os.SOTH

os.D100 <- survfit(Surv(OSY, Death) ~ DiagnosisD100, data = dat)
os.D100

os.D3compo <- survfit(Surv(OSY, Death) ~ D3compo, data = dat)
os.D3compo

osst <- survfit(Surv(OSM, Death) ~ St.AJCC2018, data = dat)
osst

osIDH1 <- survfit(Surv(OSY, Death) ~ IDH1.status, data = dat)

osIDH <- survfit(Surv(OSY, Death) ~ IDHmutation, data = dat)

fds5 <- survfit(Surv(FDS.Y, Recur2) ~ Diagnosis5, data = dat)
fds5

fds5a <- survfit(Surv(FDS.Y, Recur2) ~ Diagnosis5a, data = dat)
fds5a
fds5e4 <- survfit(Surv(FDS.Me4, Recur2) ~ Diagnosis5, data = dat)

fds5e4

fdsSOTH <- survfit(Surv(FDS.Y, Recur2) ~ Diagnosis.SOTH, data = dat)

fdsSOTH

fdsD3compo <- survfit(Surv(FDS.Y, Recur2) ~ D3compo, data = dat)

fdsD3compo

ggsurvplot(os5,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           break.time.by = 1, # Specify median survival
           legend.title = "Histological subtypes",
           xlab = "Years since operation",
           ylab = "Probability of Survival",
           xlim = c(0, 15),
           main = "Survival curve",
           font.main = 14,
           font.x = c(13, "bold"),
           font.y = c(13, "bold"),
           
           font.tickslab = 10,
           # legend.labs = c("Ductular type", "Mixed type with DPM", "Mixed type", "Small bile duct type"),
           legend = "right",
           palette = "lancet") # Change ggplot2 theme

ggsurvplot(os5a,
           pval = TRUE, conf.int = F,
           risk.table = T, # Add risk table
           # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           break.time.by = 1, # Specify median survival
           legend.title = "Histological subtypes",
           xlab = "Years since operation",
           ylab = "Probability of Survival",
           xlim = c(0, 15),
           main = "Survival curve",
           font.main = 14,
           font.x = c(13, "bold"),
           font.y = c(13, "bold"),
           
           font.tickslab = 10,
           # legend.labs = c("Ductular type", "Mixed type with DPM", "Mixed type", "Small bile duct type"),
           legend = "right",
           palette = "lancet") # Change ggplot2 theme

ggsurvplot(os.SOTH,
           pval = T, conf.int = F,
           risk.table = F,  # Add risk table
            # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           break.time.by = 1, # Specify median survival
           legend.title = "Histological subtypes",
           xlab = "Years since operation",
           ylab = "Probability of Survival",
           xlim = c(0, 15),
           main = "Survival curve",
           font.main = 14,
           font.x = c(13, "bold"),
           font.y = c(13, "bold"),
           # legend.labs = c("Ductular type", "Mixed type with DPM", "Mixed type", "Small bile duct type"),
           legend = "right",
           palette = "lancet") # Change ggplot2 theme

ggsurvplot(fds5,
           pval = T, conf.int = F,
           risk.table = F, # Add risk table
           # Change risk table color by groups
           linetype = "strata",# Change line type by groups
           break.time.by = 1, # Specify median survival
           legend.title = "Histological subtypes",
           
          xlim = c(0, 10),
           
           xlab = "Years since operation",
           ylab = "Probability of Disease-free survival",
      
           main = "Survival curve",
           font.main = 14,
           font.x = c(13, "bold"),
           font.y = c(13, "bold"),
           
           font.tickslab = 10,
          
           legend.labs = levels(dat$Diagnosis5),
           legend = "right",
           palette = "lancet") # Change ggplot2 theme

ggsurvplot(fds5a,
           pval = T, conf.int = F,
           risk.table = F, # Add risk table
           # Change risk table color by groups
           linetype = "strata",# Change line type by groups
           break.time.by = 1, # Specify median survival
           legend.title = "Histological subtypes",
           
           xlim = c(0, 10),
           
           xlab = "Years since operation",
           ylab = "Probability of Disease-free survival",
           
           main = "Survival curve",
           font.main = 14,
           font.x = c(13, "bold"),
           font.y = c(13, "bold"),
           
           font.tickslab = 10,
           
           legend.labs = levels(dat$Diagnosis5),
           legend = "right",
           palette = "lancet") # Change ggplot2 theme

summary(fds5e4)$table

ggsurvplot(fdsSOTH,
           pval = T, conf.int = F,
           risk.table = F, # Add risk table
           # Change risk table color by groups
           linetype = "strata",# Change line type by groups
           break.time.by = 1, # Specify median survival
           legend.title = "Histological subtypes",
         
           xlim = c(0, 10),
           ylab = "Probability of Disease-free survival",
           xlab = "Years since operation",
          
           
           main = "Survival curve",
           font.main = 14,
           font.x = c(13, "bold"),
           font.y = c(13, "bold"),
           
           font.tickslab = 10,
      
           legend.labs = levels(dat$Diagnosis5),
           legend = "right",
           palette = "lancet") # Change ggplot2 theme



ggsurvplot(osIDH1,
           pval = T, conf.int = F,
           risk.table = F,  # Add risk table
           # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           break.time.by = 1, # Specify median survival
           legend.title = "Histological subtypes",
           xlab = "Years since operation",
           ylab = "Probability of Survival",
           xlim = c(0, 15),
           main = "Survival curve",
           font.main = 14,
           font.x = c(13, "bold"),
           font.y = c(13, "bold"),
           # legend.labs = c("Ductular type", "Mixed type with DPM", "Mixed type", "Small bile duct type"),
           legend = "right",
           palette = "lancet") # Change ggplot2 theme

ggsurvplot(osIDH,
           pval = T, conf.int = F,
           risk.table = T,  # Add risk table
           # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           break.time.by = 1, # Specify median survival
           legend.title = "Histological subtypes",
           xlab = "Years since operation",
           ylab = "Probability of Survival",
           xlim = c(0, 15),
           main = "Survival curve",
           font.main = 14,
           font.x = c(13, "bold"),
           font.y = c(13, "bold"),
           # legend.labs = c("Ductular type", "Mixed type with DPM", "Mixed type", "Small bile duct type"),
           legend = "right",
           palette = "lancet") # Change ggplot2 theme




ggsurvplot(os.D3compo,
           pval = T, conf.int = F,
           risk.table = F,  # Add risk table
           # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           break.time.by = 1, # Specify median survival
           legend.title = "Histological subtypes",
           xlab = "Years since operation",
           ylab = "Probability of Survival",
           xlim = c(0, 15),
           main = "Survival curve",
           font.main = 14,
           font.x = c(13, "bold"),
           font.y = c(13, "bold"),
           # legend.labs = c("Ductular type", "Mixed type with DPM", "Mixed type", "Small bile duct type"),
           legend = "right",
           palette = "lancet") # Change ggplot2 theme

ggsurvplot(fdsD3compo,
           pval = T, conf.int = F,
           risk.table = F,  # Add risk table
           # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           break.time.by = 1, # Specify median survival
           legend.title = "Histological subtypes",
           xlab = "Years since operation",
           ylab = "Probability of Survival",
           xlim = c(0, 15),
           main = "Survival curve",
           font.main = 14,
           font.x = c(13, "bold"),
           font.y = c(13, "bold"),
           # legend.labs = c("Ductular type", "Mixed type with DPM", "Mixed type", "Small bile duct type"),
           legend = "right",
           palette = "lancet") # Change ggplot2 theme


# Cox regression models 

head(dat)

surv_diffa <- survdiff(Surv(OSM, Death) ~ Diagnosis5a , data = dat)

surv_diffa

surv_diffSOTH <- survdiff(Surv(OSM, Death) ~ Diagnosis.SOTH , data = dat)

surv_diffSOTH

res.cox5a <- coxph(Surv(OSM, Death==1) ~ Diagnosis5a, data =  dat)
res.cox5a
summary(res.cox5a)


res.coxT <- coxph(Surv(OSM, Death) ~ T.AJCC18Group, data =  dat)
res.coxT
summary(res.coxT)


res.coxst <- coxph(Surv(OSM, Death) ~ St.AJCC2018, data =  dat)
res.coxst
summary(res.coxst)

res.coxstgroup <- coxph(Surv(OSM, Death) ~ Stagegroup, data =  dat)
res.coxstgroup
summary(res.coxstgroup)


res.coxB.ductalspread <- coxph(Surv(OSM, Death==1) ~ B.ductalspread, data =  dat)
res.coxB.ductalspread
summary(res.coxB.ductalspread)


res.coxVascularinvasion <- coxph(Surv(OSM, Death==1) ~ Vascularinvasion, data =  dat)
res.coxVascularinvasion
summary(res.coxVascularinvasion)

res.coxNeuralinvasion <- coxph(Surv(OSM, Death==1) ~ Neuralinvasion, data =  dat)
res.coxNeuralinvasion
summary(res.coxNeuralinvasion)

res.coxIntrahepaticMets <- coxph(Surv(OSY, Death==1) ~ IntrahepaticMets, data =  dat)
res.coxIntrahepaticMets
summary(res.coxIntrahepaticMets)

res.coxF <- coxph(Surv(OSM, Death==1) ~ F.stage, data =  dat)
res.coxF
summary(res.coxF)

res.coxno.TumorSt.Mt <- coxph(Surv(OSM, Death==1) ~ no.TumorSt.Mt, data =  dat)
res.coxno.TumorSt.Mt
summary(res.coxno.TumorSt.Mt)

res.coxSenosaS <- coxph(Surv(OSM, Death==1) ~ SenosaS, data =  dat)
res.coxSenosaS
summary(res.coxSenosaS)


res.coxN<- coxph(Surv(OSM, Death==1) ~ N, data =  dat)
res.coxN
summary(res.coxN)

res.coxT.diff<- coxph(Surv(OSY, Death==1) ~ T.diffgroup, data =  dat)
res.coxT.diff
summary(res.coxT.diff)

res.coxAgegroup60<- coxph(Surv(OSM, Death==1) ~ Agegroup60, data =  dat)
res.coxAgegroup60
summary(res.coxAgegroup60)

res.coxSizegroup<- coxph(Surv(OSM, Death==1) ~ Sizegroup, data =  dat)
res.coxSizegroup
summary(res.coxSizegroup)

res.coxF.stage<- coxph(Surv(OSM, Death==1) ~ F.stage, data =  dat)
res.coxF.stage
summary(res.coxF.stage)

res.coxGrow<- coxph(Surv(OSM, Death==1) ~ Growthpattern, data =  dat)
res.coxGrow
summary(res.coxGrow)


res.coxp53.Status50 <- coxph(Surv(OSM, Death==1) ~ p53.Status50, data =  dat)
res.coxp53.Status50 
summary(res.coxp53.Status50 )

res.coxT2 <- coxph(Surv(OSM, Death) ~ T.AJCC18Group2, data =  dat)
res.coxT2
summary(res.coxT2)

res.coxgender <- coxph(Surv(OSM, Death) ~ Gender , data =  dat)
res.coxgender
summary(res.coxgender)

res.coxhepatoid <- coxph(Surv(OSM, Death) ~ Hepatoid , data =  dat)
res.coxhepatoid
summary(res.coxhepatoid)
# Multiple analysis OS 

res.coxOsmul <- coxph(Surv(OSY, Death==1) ~  Diagnosis5a  + T.AJCC18Group2  + T.diff + IntrahepaticMets + no.TumorSt.Mt,  data =  dat)
res.coxOsmul

summary(res.coxOsmul)



res.coxOsmulSOTH <- coxph(Surv(OSY, Death==1) ~  Diagnosis5a  + T.AJCC18Group2  + T.diff + IntrahepaticMets,   data =  dat)
res.coxOsmulSOTH

summary(res.coxOsmulSOTH)



# Disease free survival 


res.coxfdD5a <- coxph(Surv(FDS.M, Recur2) ~ Diagnosis5a, data =  dat)

res.coxfdD5a
summary(res.coxfdD5a)

res.coxfdstg <- coxph(Surv(FDS.M, Recur2) ~ Stagegroup, data =  dat)
res.coxfdstg
summary(res.coxfdstg)

res.coxfn <- coxph(Surv(FDS.M, Recur2) ~ N, data =  dat)
res.coxfn
summary(res.coxfn)

res.coxfim <- coxph(Surv(FDS.M, Recur2) ~ IntrahepaticMets, data =  dat)
res.coxfim
summary(res.coxfim)


res.coxfstg <- coxph(Surv(FDS.M, Recur2) ~ Stage.AJCC2018, data =  dat)
res.coxfstg

res.coxfst <- coxph(Surv(FDS.M, Recur2) ~ St.AJCC2018, data =  dat)
res.coxfst 
summary(res.coxst)

res.coxfvs <- coxph(Surv(FDS.Y, Recur2==1) ~ Vascularinvasion, data =  dat)
res.coxfvs
summary(res.coxfvs)


res.coxfneu <- coxph(Surv(FDS.Y, Recur2) ~ Neuralinvasion, data =  dat)
res.coxfneu
summary(res.coxfneu)

res.coxfgro <- coxph(Surv(FDS.M, Recur2) ~ Growthpattern , data =  dat)
res.coxfgro
summary(res.coxfgro)

res.coxftdif <- coxph(Surv(FDS.M, Recur2) ~ T.diff , data =  dat)
res.coxftdif
summary(res.coxftdif)


res.coxfdT <- coxph(Surv(FDS.Y, Recur2==1) ~ T.AJCC18Group, data =  dat)
res.coxfdT
summary(res.coxfdT)

res.coxfdB.ductalspread <- coxph(Surv(FDS.Y, Recur2==1) ~ B.ductalspread, data =  dat)
res.coxfdB.ductalspread
summary(res.coxfdB.ductalspread)


res.coxfdmul <- coxph(Surv(FDS.M, Recur2) ~ Diagnosis5a  + T.AJCC18Group2  + T.diff + Neuralinvasion + IntrahepaticMets + Growthcode,  data =  dat)

res.coxfdmul 
summary(res.coxfdmul)

res.coxD.SOTH <- coxph(Surv(OSY, Death==1) ~ Diagnosis.SOTH , data =  dat)
res.coxD.SOTH



summary(res.coxD.SOTH)

res.coxD.SOTHm <- coxph(Surv(OSY, Death==1) ~ Diagnosis.SOTH + IntrahepaticMets  + T.diff + T.AJCC18Group2, data = dat)

res.coxD.SOTHm

summary(res.coxD.SOTHm)

res.coxfdmul <- coxph(Surv(FDS.M, Recur2==1) ~ Diagnosis5a  + T.AJCC18Group2  + T.diff + IntrahepaticMets + N, data = dat)

res.coxfdmul
summary (res.coxfdmul)


res.size <- coxph(Surv(OSM, Death==1) ~ Sizegroup2,  data =  dat)
res.size

res.agegroup60 <- coxph(Surv(OSM, Death==1) ~ Agegroup60,  data =  dat)
res.agegroup60 
summary (res.agegroup60)

res.agegroup65 <- coxph(Surv(OSM, Death==1) ~ Agegroup65,  data =  dat)
res.agegroup65

res.Fstage <- coxph(Surv(FDS.Y, Recur2==1) ~ F.stage,  data =  dat)
res.Fstage

summary(res.Fstage)


res.coxfno.TumorSt.Mt <- coxph(Surv(FDS.M, Recur2==1) ~ no.TumorSt.Mt, data = dat)

res.coxfno.TumorSt.Mt
summary (res.coxfno.TumorSt.Mt)

res.coxfdT2 <- coxph(Surv(FDS.Y, Recur2) ~ T.AJCC18Group2, data =  dat)
res.coxfdT2
summary(res.coxfdT2)

res.coxfdgender <- coxph(Surv(FDS.Y, Recur2) ~ Gender , data =  dat)
res.coxfdgender
summary(res.coxfdgender)

res.fdagegroup60 <- coxph(Surv(FDS.Y, Recur2==1) ~ Agegroup60,  data =  dat)
res.fdagegroup60 
summary (res.fdagegroup60)


res.coxfdhepatoid <- coxph(Surv(FDS.Y, Death) ~ Hepatoid , data =  dat)
res.coxfdhepatoid
summary(res.coxfdhepatoid)
# Vascular analysis 

library(readxl)
vas <- read_excel("/Volumes/Macintosh.D/ICC.paper/ICC with stem cell feature/Vascular.xlsx", 
                  col_types = c("numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "numeric",
                                "text",    "numeric", "numeric", 
                                "numeric","numeric"))

# Note Diangosis5 = Diagnosis5a WHO 2019

View(vas)
head(vas)
attach(vas)
library(tidyverse)

library(ggplot2)
library(RColorBrewer)
library("ggpubr")
library(ggsignif)
display.brewer.all(colorblindFriendly = TRUE)
ordered(vas$Diagnosis.5,
                         levels = c("sbd", "mixed", "mdpm", "clc", "HCC-like"))

dat3 <-  gather(vas, "CD34.bl.mm", "CD34.IM.mm", "Cd34.tc.mm", "MDV.CD34mm", "CD34.TC.p.mm", key = "location", value = "CD34")


dat4 <- gather(vas, "Ki67.im", "Ki67.tc", key = "Ki.location", value = "Ki67")

dat5 <- gather(vas, "AVD.lb.mm",  "AVD.im.mm", "u.ADV.im.mm",  "ADV.TC.mm", "u.ADV.TC.mm", "AVD.mm", "u.AVD.mm", key = "adv.location", value = "adv")
               
dat6 <- gather(vas, "AVD.lb.mm",  "AVD.im.mm", "u.ADV.im.mm",  "ADV.TC.mm", "u.ADV.TC.mm", "AVD.mm", "u.AVD.mm", key = "HGPs.location", value = "HGPs")



GCOLOR <-  c("#82d142", "darkorange1", "deepskyblue", "blueviolet", "#247f79")

GCOLOR2 <- c("#6605f4", "#5500f9", "#34d1f0", "#a4f4f3", "#fbfea3")

g34 <-  ggplot(dat3, mapping = aes(x = location, y = CD34, fill = Diagnosis.5), ylim(0, 300)) 

cd34plot <- g34 + geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual (values = GCOLOR) +
  
 geom_point(position=position_jitterdodge(jitter.width=0.02, dodge.width = 0.8), size = ,
                                ) + 
  
  scale_color_manual (values = GCOLOR) +
  
  theme_classic() + theme(legend.position = "top")+ 
  
  stat_compare_means(aes(group = Diagnosis.5), label = "p.format") + 
  
  scale_x_discrete(limits=c("CD34.bl.mm", "MDV.CD34mm", "CD34.IM.mm", "Cd34.tc.mm")) +
  labs(title="Vascular density",x="Locations", y = "Number of vessels/mm2") + 
  
  
  geom_signif(comparisons = list(c("CD34.bl.mm", "MDV.CD34mm"),
                                 c("CD34.IM.mm", "Cd34.tc.mm"))) +
  theme( plot.title=element_text(size = 16, face = "bold"),
         text=element_text(size = 14),
         axis.text.x=element_text(colour="black", size = 12),
         axis.text.y=element_text(colour="black", size = 12))  

  
  


cd34plot



g67 <-  ggplot(dat4, mapping = aes(x = Ki.location, y = Ki67, fill = Diagnosis.5), ylim(0, 150)) 


ki67plot <-  g67 + geom_boxplot(outlier.size = 0.8) + 
  geom_point(position=position_jitterdodge(jitter.width=0.02, dodge.width = 0.8), size = 2,
                                  aes(fill= Diagnosis.5)) + 
  scale_fill_manual (values = GCOLOR) +

  scale_colour_manual(values = GCOLOR) +
  theme_classic() + theme(legend.position = "top")  + 
  stat_compare_means(aes(group = Diagnosis.5), label = "p.format") + 
  scale_x_discrete(limits=c("Ki67.im", "Ki67.tc")) +
  labs(title="Ki67 proliferation index",x="Locations", y = "Percentage of tumor cells possitive for Ki67") + 
  geom_signif(comparisons = list(c("Ki67.im", "Ki67.tc"))) +

geom_signif(y_position = c(120, 120, 120)) +
  theme( plot.title=element_text(size = 14, face = "bold"),
         text=element_text(size = 12),
         axis.text.x=element_text(colour="black", size = 12),
         axis.text.y=element_text(colour="black", size = 12))
ki67plot


gadv <-  ggplot(dat5, mapping = aes(x = adv.location, y = adv, fill = Diagnosis.5), ylim(0, 300)) 

advplot <- gadv + geom_boxplot(outlier.size = 0.8) + 
  scale_fill_manual (values = GCOLOR) +
  geom_point(position=position_jitterdodge(jitter.width=0.02, dodge.width = 0.8), size = 1.5,
                                              ) + theme_classic() + 
  theme(legend.position = "top")  + 
  
  stat_compare_means(aes(group = Diagnosis.5), label = "p.format") + 
  scale_x_discrete(limits=c("AVD.lb.mm", "AVD.mm", "AVD.im.mm",  "ADV.TC.mm","u.ADV.im.mm", "u.ADV.TC.mm")) +
  labs(title="Arterial density",x="Locations", y = "Number of arteries/mm2") + 
  geom_signif(comparisons = list(c("AVD.lb.mm", "AVD.mm"),
                                 c("AVD.im.mm",  "ADV.TC.mm"),
                                 c("u.ADV.im.mm", "u.ADV.TC.mm"))) + 
  
  geom_signif(y_position = c(40, 40, 40)) + 
  theme( plot.title=element_text(size = 20, face = "bold"),
         text=element_text(size = 16),
         axis.text.x=element_text(colour="black", size = 12),
         axis.text.y=element_text(colour="black", size = 12)) 

advplot

dat3$Fibro <- transmute(dat3, Fibro = as.character(dat3$Fibrosis))


g34f <-  ggplot(dat3, mapping = aes(x = location, y = CD34,  color = Fibro$Fibro)) 

  
cd34fplot <- g34f + geom_boxplot(outlier.shape = NA) + 
  
  geom_point(position=position_jitterdodge(jitter.width=0.02, dodge.width = 0.8), size = 1,
             aes(shape = Fibro$Fibro)) + 
  
  
  theme_classic() + theme(legend.position = "top")+ 
  
  stat_compare_means(aes(group = Fibro$Fibro), label = "p.format") + 
  
  scale_x_discrete(limits=c("CD34.bl.mm", "MDV.CD34mm", "CD34.IM.mm", "Cd34.tc.mm")) +
  labs(title="Vascular density",x="Locations", y = "Number of vessels/mm2") + 
  
  geom_signif(comparisons = list(c("CD34.bl.mm", "MDV.CD34mm"),
                                 c("CD34.IM.mm", "Cd34.tc.mm")))

cd34fplot

dat5$Fibro <- transmute(dat5, Fibro = as.character(dat5$Fibrosis))

gadvf <-  ggplot(dat5, mapping = aes(x = adv.location, y = adv, color = Fibro$Fibro)) 

advfplot <- gadvf + geom_boxplot(outlier.size = 0.8) + 
  
  geom_point(position=position_jitterdodge(jitter.width=0.02, dodge.width = 0.8), size = 0.8,
             aes(shape= Fibro$Fibro)) + theme_classic() + 
  theme(legend.position = "top")  + 
  
  stat_compare_means(aes(group = Fibro$Fibro), label = "p.format") + 
  scale_x_discrete(limits=c("AVD.lb.mm", "AVD.mm", "u.AVD.mm")) +
  labs(title="Arterial density",x="Locations", y = "Number of arteries/mm2") + 
  geom_signif(comparisons = list(c("AVD.lb.mm", "AVD.mm"))) + 
  
  geom_signif(y_position = c(40, 40, 40)) + 

theme( plot.title=element_text(size = 20),
      text=element_text(size = 16),
      axis.text.x=element_text(colour="black", size = 14),
      axis.text.y=element_text(colour="black", size = 14))

advfplot

 
wilcox.test(vas$CD34.IM, vas$Cd34.tc,
                     p.adjust.method = "BH")

wilcox.test(vas$CD34.bl, vas$MDV.CD34,
            p.adjust.method = "BH")

wilcox.test(dat3$Ki67.im, dat3$Ki67.tc,
            p.adjust.method = "BH")

wilcox.test(vas$ADV.TC.mm, vas$AVD.im.mm,
            p.adjust.method = "BH")
wilcox.test(vas$AVD.lb.mm, vas$AVD.mm,
            p.adjust.method = "BH")







library(haven)
PGC1a <- read_sav("Documents/PGC1a.sav")
View(PGC1a)
head(PGC1a)

osPGC1a <- survfit(Surv(OSM, Death) ~ PGC1a.status, data = PGC1a)

osPGC1a

dfsPGC1a <- survfit(Surv(DFSM, Recurence) ~ PGC1a.status, data = PGC1a)
dfsPGC1a


ggsurvplot(osPGC1a,
           pval = T, conf.int = F,
           risk.table = T, # Add risk table
           # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           break.time.by = 12, # Specify median survival
           legend.title = "PGC1a",
           legend.labs = levels(PGC1a$PGC1a.status),
           xlab = "Time (month)",
           ylab = "Overall survival",
           legend = "right",
           palette = "lancet") # Change ggplot2 theme

ggsurvplot(dfsPGC1a,
           pval = T, conf.int = F,
           risk.table = T, # Add risk table
           # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           break.time.by = 12, # Specify median survival
           legend.title = "PGC1a",
           legend.labs = levels(PGC1a$PGC1a.status),
           xlab = "Time (month)",
           ylab = "Probability of Disease-free survival",
           legend = "right",
           palette = "lancet") # Change ggplot2 theme


