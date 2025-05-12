##############################
# Aerobic/Anaerobe Analysis
#
#
##############################
library(tidyr)
library(ggplot2)
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


patient_tissue_o2_summary <-UNC_patients_meta_no_anti_aero_table_relevant_20200309

o2_data_long <- gather(patient_tissue_o2_summary, classification,fraction, Aerobic:Obligate.Anaerobic, factor_key=TRUE)
o2_data_long.NI <- subset(o2_data_long,o2_data_long$Pathology =="NI")

o2_data_long <- o2_data_long.NI

o2_data_long$DiseaseTissue <- paste(o2_data_long$Disease,o2_data_long$Tissue,sep=" ")
o2_data_long$DiseaseNo <- paste(o2_data_long$Disease,o2_data_long$Tissue,o2_data_long$PatientNo,sep=" ")



#### Subsets ####
o2_data_long.UNC.ileum <- subset(o2_data_long,o2_data_long$Tissue =="ileum")
o2_data_long.UNC.colon <- subset(o2_data_long,o2_data_long$Tissue =="colon")
o2_data_long.UNC.colon.IBD.nonIBD <- subset(o2_data_long.UNC.colon,o2_data_long.UNC.colon$Disease %in% c("CD","UC","nonIBD"))
o2_data_long.UNC.colon.IBD <- subset(o2_data_long.UNC.colon.IBD.nonIBD,o2_data_long.UNC.colon.IBD.nonIBD$Disease %in% c("CD","UC"))
o2_data_long.UNC.UC <- subset(o2_data_long,o2_data_long$Disease == "UC")
o2_data_long.UNC.nonIBD <- subset(o2_data_long,o2_data_long$Disease == "nonIBD")

o2_data_long.UNC.UC.colon <- subset(o2_data_long.UNC.UC,o2_data_long.UNC.UC$Tissue == "colon")

#### TESTING ####
# Tissue - nonIBD
o2_data_long.UNC.nonIBD.OAerobic <- subset(o2_data_long.UNC.nonIBD,o2_data_long.UNC.nonIBD$classification=="Obligate.Aerobic")
o2_data_long.UNC.nonIBD.Aerobic <- subset(o2_data_long.UNC.nonIBD,o2_data_long.UNC.nonIBD$classification=="Aerobic")
o2_data_long.UNC.nonIBD.FAerobic <- subset(o2_data_long.UNC.nonIBD,o2_data_long.UNC.nonIBD$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.nonIBD.Anaerobic <- subset(o2_data_long.UNC.nonIBD,o2_data_long.UNC.nonIBD$classification=="Anaerobic")
o2_data_long.UNC.nonIBD.OAnaerobic <- subset(o2_data_long.UNC.nonIBD,o2_data_long.UNC.nonIBD$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Tissue,o2_data_long.UNC.nonIBD.OAerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.nonIBD.Aerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.nonIBD.FAerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.nonIBD.Anaerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.nonIBD.OAnaerobic) 

# Tissue - UC
o2_data_long.UNC.UC.OAerobic <- subset(o2_data_long.UNC.UC,o2_data_long.UNC.UC$classification=="Obligate.Aerobic")
o2_data_long.UNC.UC.Aerobic <- subset(o2_data_long.UNC.UC,o2_data_long.UNC.UC$classification=="Aerobic")
o2_data_long.UNC.UC.FAerobic <- subset(o2_data_long.UNC.UC,o2_data_long.UNC.UC$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.UC.Anaerobic <- subset(o2_data_long.UNC.UC,o2_data_long.UNC.UC$classification=="Anaerobic")
o2_data_long.UNC.UC.OAnaerobic <- subset(o2_data_long.UNC.UC,o2_data_long.UNC.UC$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Tissue,o2_data_long.UNC.UC.OAerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.UC.Aerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.UC.FAerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.UC.Anaerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.UC.OAnaerobic) #boarder

# Ileum - Disease
o2_data_long.UNC.ileum.OAerobic <- subset(o2_data_long.UNC.ileum,o2_data_long.UNC.ileum$classification=="Obligate.Aerobic")
o2_data_long.UNC.ileum.Aerobic <- subset(o2_data_long.UNC.ileum,o2_data_long.UNC.ileum$classification=="Aerobic")
o2_data_long.UNC.ileum.FAerobic <- subset(o2_data_long.UNC.ileum,o2_data_long.UNC.ileum$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.ileum.Anaerobic <- subset(o2_data_long.UNC.ileum,o2_data_long.UNC.ileum$classification=="Anaerobic")
o2_data_long.UNC.ileum.OAnaerobic <- subset(o2_data_long.UNC.ileum,o2_data_long.UNC.ileum$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Disease,o2_data_long.UNC.ileum.OAerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.ileum.Aerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.ileum.FAerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.ileum.Anaerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.ileum.OAnaerobic) 

# Colon - Disease: nonIBD vs UC
o2_data_long.UNC.colon.UC.nonIBD <- subset(o2_data_long.UNC.colon.IBD.nonIBD,o2_data_long.UNC.colon.IBD.nonIBD$Disease %in% c("UC","nonIBD"))

o2_data_long.UNC.colon.OAerobic <- subset(o2_data_long.UNC.colon.UC.nonIBD,o2_data_long.UNC.colon.UC.nonIBD$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.Aerobic <- subset(o2_data_long.UNC.colon.UC.nonIBD,o2_data_long.UNC.colon.UC.nonIBD$classification=="Aerobic")
o2_data_long.UNC.colon.FAerobic <- subset(o2_data_long.UNC.colon.UC.nonIBD,o2_data_long.UNC.colon.UC.nonIBD$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.Anaerobic <- subset(o2_data_long.UNC.colon.UC.nonIBD,o2_data_long.UNC.colon.UC.nonIBD$classification=="Anaerobic")
o2_data_long.UNC.colon.OAnaerobic <- subset(o2_data_long.UNC.colon.UC.nonIBD,o2_data_long.UNC.colon.UC.nonIBD$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Disease,o2_data_long.UNC.colon.OAerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.Aerobic) #SIG
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.FAerobic) #SIG
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.Anaerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.OAnaerobic) #SIG

# Colon - Disease: nonIBD vs CD
o2_data_long.UNC.colon.CD.nonIBD <- subset(o2_data_long.UNC.colon.IBD.nonIBD,o2_data_long.UNC.colon.IBD.nonIBD$Disease %in% c("CD","nonIBD"))

o2_data_long.UNC.colon.OAerobic <- subset(o2_data_long.UNC.colon.CD.nonIBD,o2_data_long.UNC.colon.CD.nonIBD$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.Aerobic <- subset(o2_data_long.UNC.colon.CD.nonIBD,o2_data_long.UNC.colon.CD.nonIBD$classification=="Aerobic")
o2_data_long.UNC.colon.FAerobic <- subset(o2_data_long.UNC.colon.CD.nonIBD,o2_data_long.UNC.colon.CD.nonIBD$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.Anaerobic <- subset(o2_data_long.UNC.colon.CD.nonIBD,o2_data_long.UNC.colon.CD.nonIBD$classification=="Anaerobic")
o2_data_long.UNC.colon.OAnaerobic <- subset(o2_data_long.UNC.colon.CD.nonIBD,o2_data_long.UNC.colon.CD.nonIBD$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Disease,o2_data_long.UNC.colon.OAerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.Aerobic) 
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.FAerobic) #Boarder
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.Anaerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.OAnaerobic) 

# Colon - Disease: UC vs CD
o2_data_long.UNC.colon.IBD.OAerobic <- subset(o2_data_long.UNC.colon.IBD,o2_data_long.UNC.colon.IBD$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.IBD.Aerobic <- subset(o2_data_long.UNC.colon.IBD,o2_data_long.UNC.colon.IBD$classification=="Aerobic")
o2_data_long.UNC.colon.IBD.FAerobic <- subset(o2_data_long.UNC.colon.IBD,o2_data_long.UNC.colon.IBD$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.IBD.Anaerobic <- subset(o2_data_long.UNC.colon.IBD,o2_data_long.UNC.colon.IBD$classification=="Anaerobic")
o2_data_long.UNC.colon.IBD.OAnaerobic <- subset(o2_data_long.UNC.colon.IBD,o2_data_long.UNC.colon.IBD$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Disease,o2_data_long.UNC.colon.IBD.OAerobic) #boarder
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.IBD.Aerobic) #SIG
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.IBD.FAerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.IBD.Anaerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.IBD.OAnaerobic) #boarder



# Tissue - UC (E2v E3) colon --- FIX#####
o2_data_long.UNC.UC.colon.OAerobic <- subset(o2_data_long.UNC.UC.colon,o2_data_long.UNC.UC.colon$classification=="Obligate.Aerobic")
o2_data_long.UNC.UC.colon.Aerobic <- subset(o2_data_long.UNC.UC.colon,o2_data_long.UNC.UC.colon$classification=="Aerobic")
o2_data_long.UNC.UC.colon.FAerobic <- subset(o2_data_long.UNC.UC.colon,o2_data_long.UNC.UC.colon$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.UC.colon.Anaerobic <- subset(o2_data_long.UNC.UC.colon,o2_data_long.UNC.UC.colon$classification=="Anaerobic")
o2_data_long.UNC.UC.colon.OAnaerobic <- subset(o2_data_long.UNC.UC.colon,o2_data_long.UNC.UC.colon$classification=="Obligate.Anaerobic")

wilcox.test(fraction~TypeUC,o2_data_long.UNC.UC.colon.OAerobic)
wilcox.test(fraction~TypeUC,o2_data_long.UNC.UC.colon.Aerobic)
wilcox.test(fraction~TypeUC,o2_data_long.UNC.UC.colon.FAerobic)
wilcox.test(fraction~TypeUC,o2_data_long.UNC.UC.colon.Anaerobic)
wilcox.test(fraction~TypeUC,o2_data_long.UNC.UC.colon.OAnaerobic) 


####Tissue - UC colon - aero level###
o2_data_long.UNC.UC.colon$AeroLvl <- ifelse(o2_data_long.UNC.UC.colon$PatientNo %in% c("127","41500","42700","43000","46","52"),"Low","Health")

o2_data_long.UNC.UC.colon.OAerobic <- subset(o2_data_long.UNC.UC.colon,o2_data_long.UNC.UC.colon$classification=="Obligate.Aerobic")
o2_data_long.UNC.UC.colon.Aerobic <- subset(o2_data_long.UNC.UC.colon,o2_data_long.UNC.UC.colon$classification=="Aerobic")
o2_data_long.UNC.UC.colon.FAerobic <- subset(o2_data_long.UNC.UC.colon,o2_data_long.UNC.UC.colon$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.UC.colon.Anaerobic <- subset(o2_data_long.UNC.UC.colon,o2_data_long.UNC.UC.colon$classification=="Anaerobic")
o2_data_long.UNC.UC.colon.OAnaerobic <- subset(o2_data_long.UNC.UC.colon,o2_data_long.UNC.UC.colon$classification=="Obligate.Anaerobic")

wilcox.test(fraction~AeroLvl,o2_data_long.UNC.UC.colon.OAerobic)
wilcox.test(fraction~AeroLvl,o2_data_long.UNC.UC.colon.Aerobic)
wilcox.test(fraction~AeroLvl,o2_data_long.UNC.UC.colon.FAerobic)
wilcox.test(fraction~AeroLvl,o2_data_long.UNC.UC.colon.Anaerobic)
wilcox.test(fraction~AeroLvl,o2_data_long.UNC.UC.colon.OAnaerobic) 

#### MEDICATION USAGE - COLON FIX####
#### Colon - 5-ASA
o2_data_long.UNC.colon.UC.asa <- subset(o2_data_long.UNC.UC.colon,o2_data_long.UNC.UC.colon$X5_ASA %in% c("0","1"))

o2_data_long.UNC.colon.UC.asa.OAerobic <- subset(o2_data_long.UNC.colon.UC.asa,o2_data_long.UNC.colon.UC.asa$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.UC.asa.Aerobic <- subset(o2_data_long.UNC.colon.UC.asa,o2_data_long.UNC.colon.UC.asa$classification=="Aerobic")
o2_data_long.UNC.colon.UC.asa.FAerobic <- subset(o2_data_long.UNC.colon.UC.asa,o2_data_long.UNC.colon.UC.asa$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.UC.asa.Anaerobic <- subset(o2_data_long.UNC.colon.UC.asa,o2_data_long.UNC.colon.UC.asa$classification=="Anaerobic")
o2_data_long.UNC.colon.UC.asa.OAnaerobic <- subset(o2_data_long.UNC.colon.UC.asa,o2_data_long.UNC.colon.UC.asa$classification=="Obligate.Anaerobic")

wilcox.test(fraction~X5_ASA,o2_data_long.UNC.colon.UC.asa.OAerobic)
wilcox.test(fraction~X5_ASA,o2_data_long.UNC.colon.UC.asa.Aerobic)
wilcox.test(fraction~X5_ASA,o2_data_long.UNC.colon.UC.asa.FAerobic)
wilcox.test(fraction~X5_ASA,o2_data_long.UNC.colon.UC.asa.Anaerobic)
wilcox.test(fraction~X5_ASA,o2_data_long.UNC.colon.UC.asa.OAnaerobic)

#### Colon - Immunomodulator
o2_data_long.UNC.colon.UC.imm <- subset(o2_data_long.UNC.UC.colon,o2_data_long.UNC.UC.colon$Immunomodulators %in% c("0","1"))

o2_data_long.UNC.colon.UC.imm.OAerobic <- subset(o2_data_long.UNC.colon.UC.imm,o2_data_long.UNC.colon.UC.imm$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.UC.imm.Aerobic <- subset(o2_data_long.UNC.colon.UC.imm,o2_data_long.UNC.colon.UC.imm$classification=="Aerobic")
o2_data_long.UNC.colon.UC.imm.FAerobic <- subset(o2_data_long.UNC.colon.UC.imm,o2_data_long.UNC.colon.UC.imm$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.UC.imm.Anaerobic <- subset(o2_data_long.UNC.colon.UC.imm,o2_data_long.UNC.colon.UC.imm$classification=="Anaerobic")
o2_data_long.UNC.colon.UC.imm.OAnaerobic <- subset(o2_data_long.UNC.colon.UC.imm,o2_data_long.UNC.colon.UC.imm$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Immunomodulators,o2_data_long.UNC.colon.UC.imm.OAerobic)
wilcox.test(fraction~Immunomodulators,o2_data_long.UNC.colon.UC.imm.Aerobic)
wilcox.test(fraction~Immunomodulators,o2_data_long.UNC.colon.UC.imm.FAerobic)
wilcox.test(fraction~Immunomodulators,o2_data_long.UNC.colon.UC.imm.Anaerobic)
wilcox.test(fraction~Immunomodulators,o2_data_long.UNC.colon.UC.imm.OAnaerobic)

#### Colon - Steroids
o2_data_long.UNC.colon.UC.str <- subset(o2_data_long.UNC.UC.colon,o2_data_long.UNC.UC.colon$Steroids %in% c("0","1"))

o2_data_long.UNC.colon.UC.str.OAerobic <- subset(o2_data_long.UNC.colon.UC.str,o2_data_long.UNC.colon.UC.str$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.UC.str.Aerobic <- subset(o2_data_long.UNC.colon.UC.str,o2_data_long.UNC.colon.UC.str$classification=="Aerobic")
o2_data_long.UNC.colon.UC.str.FAerobic <- subset(o2_data_long.UNC.colon.UC.str,o2_data_long.UNC.colon.UC.str$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.UC.str.Anaerobic <- subset(o2_data_long.UNC.colon.UC.str,o2_data_long.UNC.colon.UC.str$classification=="Anaerobic")
o2_data_long.UNC.colon.UC.str.OAnaerobic <- subset(o2_data_long.UNC.colon.UC.str,o2_data_long.UNC.colon.UC.str$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Steroids,o2_data_long.UNC.colon.UC.str.OAerobic)
wilcox.test(fraction~Steroids,o2_data_long.UNC.colon.UC.str.Aerobic)
wilcox.test(fraction~Steroids,o2_data_long.UNC.colon.UC.str.FAerobic)
wilcox.test(fraction~Steroids,o2_data_long.UNC.colon.UC.str.Anaerobic)
wilcox.test(fraction~Steroids,o2_data_long.UNC.colon.UC.str.OAnaerobic)

#### Colon - anti-TNF
o2_data_long.UNC.colon.UC.tnf <- subset(o2_data_long.UNC.UC.colon,o2_data_long.UNC.UC.colon$TNF %in% c("0","1"))

o2_data_long.UNC.colon.UC.tnf.OAerobic <- subset(o2_data_long.UNC.colon.UC.tnf,o2_data_long.UNC.colon.UC.tnf$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.UC.tnf.Aerobic <- subset(o2_data_long.UNC.colon.UC.tnf,o2_data_long.UNC.colon.UC.tnf$classification=="Aerobic")
o2_data_long.UNC.colon.UC.tnf.FAerobic <- subset(o2_data_long.UNC.colon.UC.tnf,o2_data_long.UNC.colon.UC.tnf$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.UC.tnf.Anaerobic <- subset(o2_data_long.UNC.colon.UC.tnf,o2_data_long.UNC.colon.UC.tnf$classification=="Anaerobic")
o2_data_long.UNC.colon.UC.tnf.OAnaerobic <- subset(o2_data_long.UNC.colon.UC.tnf,o2_data_long.UNC.colon.UC.tnf$classification=="Obligate.Anaerobic")

wilcox.test(fraction~TNF,o2_data_long.UNC.colon.UC.tnf.OAerobic)
wilcox.test(fraction~TNF,o2_data_long.UNC.colon.UC.tnf.Aerobic)
wilcox.test(fraction~TNF,o2_data_long.UNC.colon.UC.tnf.FAerobic)
wilcox.test(fraction~TNF,o2_data_long.UNC.colon.UC.tnf.Anaerobic)
wilcox.test(fraction~TNF,o2_data_long.UNC.colon.UC.tnf.OAnaerobic)

#### Colon - Probiotics
o2_data_long.UNC.colon.UC.pro <- subset(o2_data_long.UNC.UC.colon,o2_data_long.UNC.UC.colon$Probiotic %in% c("0","1"))

o2_data_long.UNC.colon.UC.pro.OAerobic <- subset(o2_data_long.UNC.colon.UC.pro,o2_data_long.UNC.colon.UC.pro$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.UC.pro.Aerobic <- subset(o2_data_long.UNC.colon.UC.pro,o2_data_long.UNC.colon.UC.pro$classification=="Aerobic")
o2_data_long.UNC.colon.UC.pro.FAerobic <- subset(o2_data_long.UNC.colon.UC.pro,o2_data_long.UNC.colon.UC.pro$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.UC.pro.Anaerobic <- subset(o2_data_long.UNC.colon.UC.pro,o2_data_long.UNC.colon.UC.pro$classification=="Anaerobic")
o2_data_long.UNC.colon.UC.pro.OAnaerobic <- subset(o2_data_long.UNC.colon.UC.pro,o2_data_long.UNC.colon.UC.pro$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Probiotic,o2_data_long.UNC.colon.UC.pro.OAerobic)
wilcox.test(fraction~Probiotic,o2_data_long.UNC.colon.UC.pro.Aerobic)
wilcox.test(fraction~Probiotic,o2_data_long.UNC.colon.UC.pro.FAerobic)
wilcox.test(fraction~Probiotic,o2_data_long.UNC.colon.UC.pro.Anaerobic)
wilcox.test(fraction~Probiotic,o2_data_long.UNC.colon.UC.pro.OAnaerobic)

#### FIGURES ####
theme_classic()


# Colon
o2_data_long.UNC.colon.IBD.nonIBD$classification_f <- factor(o2_data_long.UNC.colon.IBD.nonIBD$classification, levels=c('Obligate.Anaerobic','Anaerobic','Facultative.Anaerobic...Microaerophilic','Aerobic','Obligate.Aerobic'))
o2_data_long.UNC.colon.IBD.nonIBD$AeroLvl <- ifelse(o2_data_long.UNC.colon.IBD.nonIBD$Disease == "UC",ifelse(o2_data_long.UNC.colon.IBD.nonIBD$PatientNo %in% c("127","41500","42700","43000","46","52"),"Low","Health"),"")
o2_data_long.UNC.colon.IBD.nonIBD$DiseaseAero <- paste(o2_data_long.UNC.colon.IBD.nonIBD$Disease,o2_data_long.UNC.colon.IBD.nonIBD$AeroLvl,sep = " ")

o2.UNC.colon <- ggplot(o2_data_long.UNC.colon.IBD.nonIBD, aes(x=Disease, y=fraction, group=DiseaseTissue)) + 
  geom_boxplot(aes(fill=DiseaseTissue),outlier.size=3) + facet_wrap(~ classification_f, ncol=3)+
  scale_fill_manual(values=c('#cb5998','#00b0f0','#cc5a49')) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.x = unit(0,"line")) 
o2.UNC.colon+theme_classic()

o2.UNC.colon.aero <- ggplot(o2_data_long.UNC.colon.IBD.nonIBD, aes(x=Disease, y=fraction, group=DiseaseAero)) + 
  geom_boxplot(aes(fill=DiseaseAero),outlier.size=3) + facet_wrap(~ classification_f, ncol=3)+
  scale_fill_manual(values=c('#cb5998','#00b0f0','#cc5a49','#d36211')) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.x = unit(0,"line")) 
o2.UNC.colon.aero+theme_classic()


####Scatterplots - Supplementary####
#nonIBD
o2.UNC.nonIBD.OAerobic <- ggplot(o2_data_long.UNC.nonIBD.OAerobic , aes(x=DiseaseNo, y=fraction,color=Tissue)) + 
  geom_point(size=3) +
  labs(y = "Relative Abundance")+
  scale_color_manual(values=c('#00b0f0','#8b6cc9'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.nonIBD.Aerobic <- ggplot(o2_data_long.UNC.nonIBD.Aerobic , aes(x=DiseaseNo, y=fraction,color=Tissue)) + 
  geom_point(size=3) +
  labs(y = "Relative Abundance")+
  scale_color_manual(values=c('#00b0f0','#8b6cc9'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.nonIBD.FAerobic <- ggplot(o2_data_long.UNC.nonIBD.FAerobic , aes(x=DiseaseNo, y=fraction,color=Tissue)) + 
  geom_point(size=3) +
  labs(y = "Relative Abundance")+
  scale_color_manual(values=c('#00b0f0','#8b6cc9'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.nonIBD.Anaerobic <- ggplot(o2_data_long.UNC.nonIBD.Anaerobic , aes(x=DiseaseNo, y=fraction,color=Tissue)) + 
  geom_point(size=3) +
  labs(y = "Relative Abundance")+
  scale_color_manual(values=c('#00b0f0','#8b6cc9'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.nonIBD.OAnaerobic <- ggplot(o2_data_long.UNC.nonIBD.OAnaerobic , aes(x=DiseaseNo, y=fraction,color=Tissue)) + 
  geom_point(size=3) +
  labs(y = "Relative Abundance")+
  scale_color_manual(values=c('#00b0f0','#8b6cc9'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.nonIBD.OAerobic
o2.UNC.nonIBD.Aerobic
o2.UNC.nonIBD.FAerobic
o2.UNC.nonIBD.Anaerobic
o2.UNC.nonIBD.OAnaerobic

nonIBD_scatter_figure <- ggarrange(o2.UNC.nonIBD.OAerobic + rremove("xlab"),o2.UNC.nonIBD.Aerobic + rremove("xlab") , o2.UNC.nonIBD.FAerobic,o2.UNC.nonIBD.Anaerobic + rremove("xlab"),o2.UNC.nonIBD.OAnaerobic + rremove("xlab"),
                                   common.legend = TRUE,
                                   labels = c("A", "B", "C","D","E"),
                                   ncol = 3, nrow = 2)
annotate_figure(nonIBD_scatter_figure,
                bottom = text_grob("Aerotolerance Profiles in nonIBD Mucosa", color = "black", face = "bold", size = 14),
                fig.lab = "Figure S2", fig.lab.face = "bold",fig.lab.pos = "bottom.left"
)

#CD - Tissue
o2.UNC.CD.OAerobic <- ggplot(o2_data_long.UNC.CD.OAerobic , aes(x=DiseaseNo, y=fraction,color=Tissue)) + 
  geom_point(size=3) +
  labs(y = "Relative Abundance")+
  scale_color_manual(values=c('#cb5998','#cc5b48'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.CD.Aerobic <- ggplot(o2_data_long.UNC.CD.Aerobic , aes(x=DiseaseNo, y=fraction,color=Tissue)) + 
  geom_point(size=3) +
  labs(y = "Relative Abundance")+
  scale_color_manual(values=c('#cb5998','#cc5b48'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.CD.FAerobic <- ggplot(o2_data_long.UNC.CD.FAerobic , aes(x=DiseaseNo, y=fraction,color=Tissue)) + 
  geom_point(size=3) +
  labs(y = "Relative Abundance")+
  scale_color_manual(values=c('#cb5998','#cc5b48'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.CD.Anaerobic <- ggplot(o2_data_long.UNC.CD.Anaerobic , aes(x=DiseaseNo, y=fraction,color=Tissue)) + 
  geom_point(size=3) +
  labs(y = "Relative Abundance")+
  scale_color_manual(values=c('#cb5998','#cc5b48'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.CD.OAnaerobic <- ggplot(o2_data_long.UNC.CD.OAnaerobic , aes(x=DiseaseNo, y=fraction,color=Tissue)) + 
  geom_point(size=3) +
  labs(y = "Relative Abundance")+
  scale_color_manual(values=c('#cb5998','#cc5b48'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.CD.OAerobic
o2.UNC.CD.Aerobic
o2.UNC.CD.FAerobic
o2.UNC.CD.Anaerobic
o2.UNC.CD.OAnaerobic

CD_scatter_figure <- ggarrange(o2.UNC.CD.OAerobic + rremove("xlab"),o2.UNC.CD.Aerobic + rremove("xlab") , o2.UNC.CD.FAerobic,o2.UNC.CD.Anaerobic + rremove("xlab"),o2.UNC.CD.OAnaerobic + rremove("xlab"),
                               common.legend = TRUE,
                               labels = c("A", "B", "C","D","E"),
                               ncol = 3, nrow = 2)
annotate_figure(CD_scatter_figure,
                bottom = text_grob("Aerotolerance Profiles in CD Mucosa", color = "black", face = "bold", size = 14),
                fig.lab = "Figure S4", fig.lab.face = "bold",fig.lab.pos = "bottom.left"
)
#Colon - Disease
o2_data_long.UNC.colon.OAnaerobic$DiseaseNo <- as.factor(o2_data_long.UNC.colon.OAnaerobic$PatientNo)

o2.UNC.colon.OAerobic <- ggplot(o2_data_long.UNC.colon.OAerobic , aes(x=DiseaseNo, y=fraction,color=Disease)) + 
  geom_point(size=5) +
  labs(y = "Relative Abundance")+
  scale_color_manual(values=c('#cb5998','#00b0f0'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
o2_data_long.UNC.colon.Aerobic$DiseaseNo <- as.factor(o2_data_long.UNC.colon.Aerobic$PatientNo)
o2.UNC.colon.Aerobic <- ggplot(o2_data_long.UNC.colon.Aerobic , aes(x=DiseaseNo, y=fraction,color=Disease)) + 
  geom_point(size=5) +
  labs(y = "Relative Abundance")+
  scale_color_manual(values=c('#cb5998','#00b0f0'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.colon.FAerobic <- ggplot(o2_data_long.UNC.colon.FAerobic , aes(x=DiseaseNo, y=fraction,color=Disease)) + 
  geom_point(size=5) +
  labs(y = "Relative Abundance")+
  scale_color_manual(values=c('#cb5998','#00b0f0'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.colon.Anaerobic <- ggplot(o2_data_long.UNC.colon.Anaerobic , aes(x=DiseaseNo, y=fraction,color=Disease)) + 
  geom_point(size=5) +
  labs(y = "Relative Abundance")+
  scale_color_manual(values=c('#cb5998','#00b0f0'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.colon.OAnaerobic <- ggplot(o2_data_long.UNC.colon.OAnaerobic , aes(x=DiseaseNo, y=fraction,color=Disease)) + 
  geom_point(size=5) +
  labs(y = "Relative Abundance")+
  scale_color_manual(values=c('#cb5998','#00b0f0'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.colon.OAerobic
o2.UNC.colon.Aerobic
o2.UNC.colon.FAerobic
o2.UNC.colon.Anaerobic
o2.UNC.colon.OAnaerobic

colon_scatter_figure <- ggarrange(o2.UNC.colon.OAerobic + rremove("xlab"),o2.UNC.colon.Aerobic + rremove("xlab") , o2.UNC.colon.FAerobic,o2.UNC.colon.Anaerobic + rremove("xlab"),o2.UNC.colon.OAnaerobic + rremove("xlab"),
                                  common.legend = TRUE,
                                  labels = c("A", "B", "C","D","E"),
                                  ncol = 3, nrow = 2)
annotate_figure(colon_scatter_figure,
                bottom = text_grob("Aerotolerance Profiles in Colonic Mucosa", color = "black", face = "bold", size = 14),
                fig.lab = "Figure S1", fig.lab.face = "bold",fig.lab.pos = "bottom.left"
)

#Ileum - Disease
o2.UNC.ileum.OAerobic <- ggplot(o2_data_long.UNC.ileum.OAerobic , aes(x=DiseaseNo, y=fraction,color=Disease)) + 
  geom_point(size=5) +
  labs(y = "Relative Abundance",caption = "Figure S2")+
  scale_color_manual(values=c('#cc5b48','#8b6cc9'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.ileum.Aerobic <- ggplot(o2_data_long.UNC.ileum.Aerobic , aes(x=DiseaseNo, y=fraction,color=Disease)) + 
  geom_point(size=5) +
  scale_color_manual(values=c('#cc5b48','#8b6cc9'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.ileum.FAerobic <- ggplot(o2_data_long.UNC.ileum.FAerobic , aes(x=DiseaseNo, y=fraction,color=Disease)) + 
  geom_point(size=5) +
  scale_color_manual(values=c('#cc5b48','#8b6cc9'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.ileum.Anaerobic <- ggplot(o2_data_long.UNC.ileum.Anaerobic , aes(x=DiseaseNo, y=fraction,color=Disease)) + 
  geom_point(size=5) +
  scale_color_manual(values=c('#cc5b48','#8b6cc9'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

o2.UNC.ileum.OAnaerobic <- ggplot(o2_data_long.UNC.ileum.OAnaerobic , aes(x=DiseaseNo, y=fraction,color=Disease)) + 
  geom_point(size=5) +
  scale_color_manual(values=c('#cc5b48','#8b6cc9'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
o2.UNC.ileum.OAerobic


ileum_scatter_figure <- ggarrange(o2.UNC.ileum.OAerobic + rremove("xlab"),o2.UNC.ileum.Aerobic + rremove("xlab") , o2.UNC.ileum.FAerobic,o2.UNC.ileum.Anaerobic + rremove("xlab"),o2.UNC.ileum.OAnaerobic + rremove("xlab"),
                                  common.legend = TRUE,
                                  labels = c("A", "B", "C","D","E"),
                                  ncol = 3, nrow = 2)
annotate_figure(ileum_scatter_figure,
                bottom = text_grob("Aerotolerance Profiles in Ileal Mucosa", color = "black", face = "bold", size = 14),
                fig.lab = "Figure S3", fig.lab.face = "bold",fig.lab.pos = "bottom.left"
)
