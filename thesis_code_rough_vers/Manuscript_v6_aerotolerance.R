##############################
# Aerobic/Anaerobe Analysis
#
#
##############################
library(tidyr)
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


patient_tissue_o2_summary <- UNC_patients_meta_no_anti_relevant_nonIBD_CD_aero_table_20190720

o2_data_long <- gather(patient_tissue_o2_summary, classification,fraction, Aerobic:Obligate.Anaerobic, factor_key=TRUE)
o2_data_long.NI <- subset(o2_data_long,o2_data_long$Pathology =="NI")

o2_data_long <- o2_data_long.NI

o2_data_long$DiseaseTissue <- paste(o2_data_long$Disease,o2_data_long$Tissue,sep=" ")
o2_data_long$DiseaseNo <- paste(o2_data_long$Disease,o2_data_long$Tissue,o2_data_long$PatientNo,sep=" ")



#### Subsets ####
o2_data_long.UNC.CD <- subset(o2_data_long, o2_data_long$Disease == "CD")
o2_data_long.UNC.nonIBD <- subset(o2_data_long, o2_data_long$Disease == "nonIBD")

o2_data_long.UNC.ileum <- subset(o2_data_long,o2_data_long$Tissue =="ileum")
o2_data_long.UNC.colon <- subset(o2_data_long,o2_data_long$Tissue =="colon")


o2_data_long.UNC.colon.CD <- subset(o2_data_long.UNC.colon,o2_data_long.UNC.colon$Disease == "CD")
o2_data_long.UNC.ileum.CD <- subset(o2_data_long.UNC.ileum, o2_data_long.UNC.ileum$Disease == "CD")

o2_data_long.UNC.colon.nonIBD <- subset(o2_data_long.UNC.colon,o2_data_long.UNC.colon$Disease == "nonIBD")
o2_data_long.UNC.ileum.nonIBD <- subset(o2_data_long.UNC.ileum,o2_data_long.UNC.ileum$Disease == "nonIBD")

o2_data_long.UNC.ileum.CD.recur <- subset(o2_data_long.UNC.ileum.CD,o2_data_long.UNC.ileum.CD$Recurrence %in% c("no recurrence","recurrence"))
o2_data_long.UNC.colon.CD.recur <- subset(o2_data_long.UNC.colon.CD,o2_data_long.UNC.colon.CD$Recurrence %in% c("no recurrence","recurrence"))

o2_data_long.UNC.colon.CD.subtype <- subset(o2_data_long.UNC.colon.CD,o2_data_long.UNC.colon.CD$Subtype %in% c("IL","CL"))
o2_data_long.UNC.ileum.CD.subtype <- subset(o2_data_long.UNC.ileum.CD,o2_data_long.UNC.ileum.CD$Subtype %in% c("IL","CL"))


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

# Tissue - CD
o2_data_long.UNC.CD.OAerobic <- subset(o2_data_long.UNC.CD,o2_data_long.UNC.CD$classification=="Obligate.Aerobic")
o2_data_long.UNC.CD.Aerobic <- subset(o2_data_long.UNC.CD,o2_data_long.UNC.CD$classification=="Aerobic")
o2_data_long.UNC.CD.FAerobic <- subset(o2_data_long.UNC.CD,o2_data_long.UNC.CD$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.CD.Anaerobic <- subset(o2_data_long.UNC.CD,o2_data_long.UNC.CD$classification=="Anaerobic")
o2_data_long.UNC.CD.OAnaerobic <- subset(o2_data_long.UNC.CD,o2_data_long.UNC.CD$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.OAerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.Aerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.FAerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.Anaerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.OAnaerobic) 

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
wilcox.test(fraction~Disease,o2_data_long.UNC.ileum.OAnaerobic) #SIGNIF

# Colon - Disease

o2_data_long.UNC.colon.OAerobic <- subset(o2_data_long.UNC.colon,o2_data_long.UNC.colon$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.Aerobic <- subset(o2_data_long.UNC.colon,o2_data_long.UNC.colon$classification=="Aerobic")
o2_data_long.UNC.colon.FAerobic <- subset(o2_data_long.UNC.colon,o2_data_long.UNC.colon$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.Anaerobic <- subset(o2_data_long.UNC.colon,o2_data_long.UNC.colon$classification=="Anaerobic")
o2_data_long.UNC.colon.OAnaerobic <- subset(o2_data_long.UNC.colon,o2_data_long.UNC.colon$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Disease,o2_data_long.UNC.colon.OAerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.Aerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.FAerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.Anaerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.OAnaerobic) 

# Tissue - CD (B2, B3)
o2_data_long.UNC.CD.B23.OAerobic <- subset(o2_data_long.UNC.CD.OAerobic,o2_data_long.UNC.CD.OAerobic$Behavior %in% c("B2","B3"))
o2_data_long.UNC.CD.B23.Aerobic <- subset(o2_data_long.UNC.CD.Aerobic,o2_data_long.UNC.CD.Aerobic$Behavior %in% c("B2","B3"))
o2_data_long.UNC.CD.B23.FAerobic <- subset(o2_data_long.UNC.CD.FAerobic,o2_data_long.UNC.CD.FAerobic$Behavior %in% c("B2","B3"))
o2_data_long.UNC.CD.B23.Anaerobic <- subset(o2_data_long.UNC.CD.Anaerobic,o2_data_long.UNC.CD.Anaerobic$Behavior %in% c("B2","B3"))
o2_data_long.UNC.CD.B23.OAnaerobic <- subset(o2_data_long.UNC.CD.OAnaerobic,o2_data_long.UNC.CD.OAnaerobic$Behavior %in% c("B2","B3"))

wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.B23.OAerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.B23.Aerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.B23.FAerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.B23.Anaerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.B23.OAnaerobic) 

#### Aerotolerance - Disease Behavior - Colon ####
o2_data_long.UNC.colon.CD.B12 <- subset(o2_data_long.UNC.colon.CD,o2_data_long.UNC.colon.CD$Behavior %in% c("B1","B2"))
o2_data_long.UNC.colon.CD.B13 <- subset(o2_data_long.UNC.colon.CD,o2_data_long.UNC.colon.CD$Behavior %in% c("B1","B3"))
o2_data_long.UNC.colon.CD.B23 <- subset(o2_data_long.UNC.colon.CD,o2_data_long.UNC.colon.CD$Behavior %in% c("B2","B3"))

#### CD Colon - B1 vs B2
o2_data_long.UNC.colon.CD.B12.OAerobic <- subset(o2_data_long.UNC.colon.CD.B12,o2_data_long.UNC.colon.CD.B12$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.CD.B12.Aerobic <- subset(o2_data_long.UNC.colon.CD.B12,o2_data_long.UNC.colon.CD.B12$classification=="Aerobic")
o2_data_long.UNC.colon.CD.B12.FAerobic <- subset(o2_data_long.UNC.colon.CD.B12,o2_data_long.UNC.colon.CD.B12$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.CD.B12.Anaerobic <- subset(o2_data_long.UNC.colon.CD.B12,o2_data_long.UNC.colon.CD.B12$classification=="Anaerobic")
o2_data_long.UNC.colon.CD.B12.OAnaerobic <- subset(o2_data_long.UNC.colon.CD.B12,o2_data_long.UNC.colon.CD.B12$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.CD.B12.OAerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.CD.B12.Aerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.CD.B12.FAerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.CD.B12.Anaerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.CD.B12.OAnaerobic) 

#### CD Colon - B1 vs B3
o2_data_long.UNC.colon.CD.B13.OAerobic <- subset(o2_data_long.UNC.colon.CD.B13,o2_data_long.UNC.colon.CD.B13$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.CD.B13.Aerobic <- subset(o2_data_long.UNC.colon.CD.B13,o2_data_long.UNC.colon.CD.B13$classification=="Aerobic")
o2_data_long.UNC.colon.CD.B13.FAerobic <- subset(o2_data_long.UNC.colon.CD.B13,o2_data_long.UNC.colon.CD.B13$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.CD.B13.Anaerobic <- subset(o2_data_long.UNC.colon.CD.B13,o2_data_long.UNC.colon.CD.B13$classification=="Anaerobic")
o2_data_long.UNC.colon.CD.B13.OAnaerobic <- subset(o2_data_long.UNC.colon.CD.B13,o2_data_long.UNC.colon.CD.B13$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.CD.B13.OAerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.CD.B13.Aerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.CD.B13.FAerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.CD.B13.Anaerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.CD.B13.OAnaerobic)

#### CD Colon - B2 vs B3
o2_data_long.UNC.colon.CD.B23.OAerobic <- subset(o2_data_long.UNC.colon.CD.B23,o2_data_long.UNC.colon.CD.B23$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.CD.B23.Aerobic <- subset(o2_data_long.UNC.colon.CD.B23,o2_data_long.UNC.colon.CD.B23$classification=="Aerobic")
o2_data_long.UNC.colon.CD.B23.FAerobic <- subset(o2_data_long.UNC.colon.CD.B23,o2_data_long.UNC.colon.CD.B23$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.CD.B23.Anaerobic <- subset(o2_data_long.UNC.colon.CD.B23,o2_data_long.UNC.colon.CD.B23$classification=="Anaerobic")
o2_data_long.UNC.colon.CD.B23.OAnaerobic <- subset(o2_data_long.UNC.colon.CD.B23,o2_data_long.UNC.colon.CD.B23$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.CD.B23.OAerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.CD.B23.Aerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.CD.B23.FAerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.CD.B23.Anaerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.CD.B23.OAnaerobic) 

#### Colon - nonIBD vs CD B1
o2_data_long.UNC.colon$CatB <- paste(o2_data_long.UNC.colon$Disease,o2_data_long.UNC.colon$Behavior,sep=" ")
o2_data_long.UNC.colon.nonIBD.CD.B1 <- subset(o2_data_long.UNC.colon,o2_data_long.UNC.colon$CatB %in% c("nonIBD ", "CD B1"))

o2_data_long.UNC.colon.nonIBD.CD.B1.OAerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.B1,o2_data_long.UNC.colon.nonIBD.CD.B1$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.nonIBD.CD.B1.Aerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.B1,o2_data_long.UNC.colon.nonIBD.CD.B1$classification=="Aerobic")
o2_data_long.UNC.colon.nonIBD.CD.B1.FAerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.B1,o2_data_long.UNC.colon.nonIBD.CD.B1$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.nonIBD.CD.B1.Anaerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.B1,o2_data_long.UNC.colon.nonIBD.CD.B1$classification=="Anaerobic")
o2_data_long.UNC.colon.nonIBD.CD.B1.OAnaerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.B1,o2_data_long.UNC.colon.nonIBD.CD.B1$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.nonIBD.CD.B1.OAerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.nonIBD.CD.B1.Aerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.nonIBD.CD.B1.FAerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.nonIBD.CD.B1.Anaerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.nonIBD.CD.B1.OAnaerobic) 

#### Colon - nonIBD vs CD B2
o2_data_long.UNC.colon.nonIBD.CD.B2 <- subset(o2_data_long.UNC.colon,o2_data_long.UNC.colon$CatB %in% c("nonIBD ", "CD B2"))

o2_data_long.UNC.colon.nonIBD.CD.B2.OAerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.B2,o2_data_long.UNC.colon.nonIBD.CD.B2$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.nonIBD.CD.B2.Aerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.B2,o2_data_long.UNC.colon.nonIBD.CD.B2$classification=="Aerobic")
o2_data_long.UNC.colon.nonIBD.CD.B2.FAerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.B2,o2_data_long.UNC.colon.nonIBD.CD.B2$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.nonIBD.CD.B2.Anaerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.B2,o2_data_long.UNC.colon.nonIBD.CD.B2$classification=="Anaerobic")
o2_data_long.UNC.colon.nonIBD.CD.B2.OAnaerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.B2,o2_data_long.UNC.colon.nonIBD.CD.B2$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.nonIBD.CD.B2.OAerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.nonIBD.CD.B2.Aerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.nonIBD.CD.B2.FAerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.nonIBD.CD.B2.Anaerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.nonIBD.CD.B2.OAnaerobic)

#### Colon - nonIBD vs CD B3
o2_data_long.UNC.colon.nonIBD.CD.B3 <- subset(o2_data_long.UNC.colon,o2_data_long.UNC.colon$CatB %in% c("nonIBD ", "CD B3"))

o2_data_long.UNC.colon.nonIBD.CD.B3.OAerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.B3,o2_data_long.UNC.colon.nonIBD.CD.B3$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.nonIBD.CD.B3.Aerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.B3,o2_data_long.UNC.colon.nonIBD.CD.B3$classification=="Aerobic")
o2_data_long.UNC.colon.nonIBD.CD.B3.FAerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.B3,o2_data_long.UNC.colon.nonIBD.CD.B3$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.nonIBD.CD.B3.Anaerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.B3,o2_data_long.UNC.colon.nonIBD.CD.B3$classification=="Anaerobic")
o2_data_long.UNC.colon.nonIBD.CD.B3.OAnaerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.B3,o2_data_long.UNC.colon.nonIBD.CD.B3$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.nonIBD.CD.B3.OAerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.nonIBD.CD.B3.Aerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.nonIBD.CD.B3.FAerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.nonIBD.CD.B3.Anaerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.colon.nonIBD.CD.B3.OAnaerobic)

#### Aerotolerance - Disease Behavior - Ileum ####
o2_data_long.UNC.ileum.CD.B23 <- subset(o2_data_long.UNC.ileum.CD,o2_data_long.UNC.ileum.CD$Behavior %in% c("B2","B3"))

#### CD Ileum - B2 vs B3
o2_data_long.UNC.ileum.CD.B23.OAerobic <- subset(o2_data_long.UNC.ileum.CD.B23,o2_data_long.UNC.ileum.CD.B23$classification=="Obligate.Aerobic")
o2_data_long.UNC.ileum.CD.B23.Aerobic <- subset(o2_data_long.UNC.ileum.CD.B23,o2_data_long.UNC.ileum.CD.B23$classification=="Aerobic")
o2_data_long.UNC.ileum.CD.B23.FAerobic <- subset(o2_data_long.UNC.ileum.CD.B23,o2_data_long.UNC.ileum.CD.B23$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.ileum.CD.B23.Anaerobic <- subset(o2_data_long.UNC.ileum.CD.B23,o2_data_long.UNC.ileum.CD.B23$classification=="Anaerobic")
o2_data_long.UNC.ileum.CD.B23.OAnaerobic <- subset(o2_data_long.UNC.ileum.CD.B23,o2_data_long.UNC.ileum.CD.B23$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Behavior,o2_data_long.UNC.ileum.CD.B23.OAerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.ileum.CD.B23.Aerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.ileum.CD.B23.FAerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.ileum.CD.B23.Anaerobic)
wilcox.test(fraction~Behavior,o2_data_long.UNC.ileum.CD.B23.OAnaerobic) 

#### Ileum - nonIBD vs CD B2
o2_data_long.UNC.ileum$CatB <- paste(o2_data_long.UNC.ileum$Disease,o2_data_long.UNC.ileum$Behavior,sep=" ")
o2_data_long.UNC.ileum.nonIBD.CD.B2 <- subset(o2_data_long.UNC.ileum,o2_data_long.UNC.ileum$CatB %in% c("nonIBD ", "CD B2"))
o2_data_long.UNC.ileum.nonIBD.CD.B3 <- subset(o2_data_long.UNC.ileum,o2_data_long.UNC.ileum$CatB %in% c("nonIBD ", "CD B3"))

o2_data_long.UNC.ileum.nonIBD.CD.B2.OAerobic <- subset(o2_data_long.UNC.ileum.nonIBD.CD.B2,o2_data_long.UNC.ileum.nonIBD.CD.B2$classification=="Obligate.Aerobic")
o2_data_long.UNC.ileum.nonIBD.CD.B2.Aerobic <- subset(o2_data_long.UNC.ileum.nonIBD.CD.B2,o2_data_long.UNC.ileum.nonIBD.CD.B2$classification=="Aerobic")
o2_data_long.UNC.ileum.nonIBD.CD.B2.FAerobic <- subset(o2_data_long.UNC.ileum.nonIBD.CD.B2,o2_data_long.UNC.ileum.nonIBD.CD.B2$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.ileum.nonIBD.CD.B2.Anaerobic <- subset(o2_data_long.UNC.ileum.nonIBD.CD.B2,o2_data_long.UNC.ileum.nonIBD.CD.B2$classification=="Anaerobic")
o2_data_long.UNC.ileum.nonIBD.CD.B2.OAnaerobic <- subset(o2_data_long.UNC.ileum.nonIBD.CD.B2,o2_data_long.UNC.ileum.nonIBD.CD.B2$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Disease,o2_data_long.UNC.ileum.nonIBD.CD.B2.OAerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.ileum.nonIBD.CD.B2.Aerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.ileum.nonIBD.CD.B2.FAerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.ileum.nonIBD.CD.B2.Anaerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.ileum.nonIBD.CD.B2.OAnaerobic)

#### Ileum - nonIBD vs CD B3

o2_data_long.UNC.ileum.nonIBD.CD.B3.OAerobic <- subset(o2_data_long.UNC.ileum.nonIBD.CD.B3,o2_data_long.UNC.ileum.nonIBD.CD.B3$classification=="Obligate.Aerobic")
o2_data_long.UNC.ileum.nonIBD.CD.B3.Aerobic <- subset(o2_data_long.UNC.ileum.nonIBD.CD.B3,o2_data_long.UNC.ileum.nonIBD.CD.B3$classification=="Aerobic")
o2_data_long.UNC.ileum.nonIBD.CD.B3.FAerobic <- subset(o2_data_long.UNC.ileum.nonIBD.CD.B3,o2_data_long.UNC.ileum.nonIBD.CD.B3$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.ileum.nonIBD.CD.B3.Anaerobic <- subset(o2_data_long.UNC.ileum.nonIBD.CD.B3,o2_data_long.UNC.ileum.nonIBD.CD.B3$classification=="Anaerobic")
o2_data_long.UNC.ileum.nonIBD.CD.B3.OAnaerobic <- subset(o2_data_long.UNC.ileum.nonIBD.CD.B3,o2_data_long.UNC.ileum.nonIBD.CD.B3$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Disease,o2_data_long.UNC.ileum.nonIBD.CD.B3.OAerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.ileum.nonIBD.CD.B3.Aerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.ileum.nonIBD.CD.B3.FAerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.ileum.nonIBD.CD.B3.Anaerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.ileum.nonIBD.CD.B3.OAnaerobic)

#### Disease Behavior - Tissue ####
o2_data_long.UNC.CD.B2 <- subset(o2_data_long.UNC.CD,o2_data_long.UNC.CD$Behavior=="B2")
o2_data_long.UNC.CD.B3 <- subset(o2_data_long.UNC.CD,o2_data_long.UNC.CD$Behavior=="B3")

#### CD B2 - Tissue
o2_data_long.UNC.CD.B2.OAerobic <- subset(o2_data_long.UNC.CD.B2,o2_data_long.UNC.CD.B2$classification=="Obligate.Aerobic")
o2_data_long.UNC.CD.B2.Aerobic <- subset(o2_data_long.UNC.CD.B2,o2_data_long.UNC.CD.B2$classification=="Aerobic")
o2_data_long.UNC.CD.B2.FAerobic <- subset(o2_data_long.UNC.CD.B2,o2_data_long.UNC.CD.B2$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.CD.B2.Anaerobic <- subset(o2_data_long.UNC.CD.B2,o2_data_long.UNC.CD.B2$classification=="Anaerobic")
o2_data_long.UNC.CD.B2.OAnaerobic <- subset(o2_data_long.UNC.CD.B2,o2_data_long.UNC.CD.B2$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.B2.OAerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.B2.Aerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.B2.FAerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.B2.Anaerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.B2.OAnaerobic)

#### CD B3 - Tissue
o2_data_long.UNC.CD.B3.OAerobic <- subset(o2_data_long.UNC.CD.B3,o2_data_long.UNC.CD.B3$classification=="Obligate.Aerobic")
o2_data_long.UNC.CD.B3.Aerobic <- subset(o2_data_long.UNC.CD.B3,o2_data_long.UNC.CD.B3$classification=="Aerobic")
o2_data_long.UNC.CD.B3.FAerobic <- subset(o2_data_long.UNC.CD.B3,o2_data_long.UNC.CD.B3$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.CD.B3.Anaerobic <- subset(o2_data_long.UNC.CD.B3,o2_data_long.UNC.CD.B3$classification=="Anaerobic")
o2_data_long.UNC.CD.B3.OAnaerobic <- subset(o2_data_long.UNC.CD.B3,o2_data_long.UNC.CD.B3$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.B3.OAerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.B3.Aerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.B3.FAerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.B3.Anaerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.B3.OAnaerobic)


#### Aerotolerance - Disease Location ####
o2_data_long.UNC.colon.CD.L23 <- subset(o2_data_long.UNC.colon.CD, o2_data_long.UNC.colon.CD$Location %in% c("L2","L3"))
o2_data_long.UNC.colon$CatL <- paste(o2_data_long.UNC.colon$Disease,o2_data_long.UNC.colon$Location,sep = " ")
o2_data_long.UNC.colon.nonIBD.CD.L2 <- subset(o2_data_long.UNC.colon,o2_data_long.UNC.colon$CatL %in% c("nonIBD ","CD L2"))
o2_data_long.UNC.colon.nonIBD.CD.L3 <- subset(o2_data_long.UNC.colon,o2_data_long.UNC.colon$CatL %in% c("nonIBD ","CD L3"))

#### CD Colon - L2 vs L3
o2_data_long.UNC.colon.CD.L23.OAerobic <- subset(o2_data_long.UNC.colon.CD.L23,o2_data_long.UNC.colon.CD.L23$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.CD.L23.Aerobic <- subset(o2_data_long.UNC.colon.CD.L23,o2_data_long.UNC.colon.CD.L23$classification=="Aerobic")
o2_data_long.UNC.colon.CD.L23.FAerobic <- subset(o2_data_long.UNC.colon.CD.L23,o2_data_long.UNC.colon.CD.L23$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.CD.L23.Anaerobic <- subset(o2_data_long.UNC.colon.CD.L23,o2_data_long.UNC.colon.CD.L23$classification=="Anaerobic")
o2_data_long.UNC.colon.CD.L23.OAnaerobic <- subset(o2_data_long.UNC.colon.CD.L23,o2_data_long.UNC.colon.CD.L23$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Location,o2_data_long.UNC.colon.CD.L23.OAerobic)
wilcox.test(fraction~Location,o2_data_long.UNC.colon.CD.L23.Aerobic)
wilcox.test(fraction~Location,o2_data_long.UNC.colon.CD.L23.FAerobic)
wilcox.test(fraction~Location,o2_data_long.UNC.colon.CD.L23.Anaerobic)
wilcox.test(fraction~Location,o2_data_long.UNC.colon.CD.L23.OAnaerobic)

#### Colon - nonIBD vs CD L2
o2_data_long.UNC.colon.nonIBD.CD.L2.OAerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.L2,o2_data_long.UNC.colon.nonIBD.CD.L2$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.nonIBD.CD.L2.Aerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.L2,o2_data_long.UNC.colon.nonIBD.CD.L2$classification=="Aerobic")
o2_data_long.UNC.colon.nonIBD.CD.L2.FAerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.L2,o2_data_long.UNC.colon.nonIBD.CD.L2$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.nonIBD.CD.L2.Anaerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.L2,o2_data_long.UNC.colon.nonIBD.CD.L2$classification=="Anaerobic")
o2_data_long.UNC.colon.nonIBD.CD.L2.OAnaerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.L2,o2_data_long.UNC.colon.nonIBD.CD.L2$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Disease,o2_data_long.UNC.colon.nonIBD.CD.L2.OAerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.nonIBD.CD.L2.Aerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.nonIBD.CD.L2.FAerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.nonIBD.CD.L2.Anaerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.nonIBD.CD.L2.OAnaerobic)

#### Colon - nonIBD vs CD L3
o2_data_long.UNC.colon.nonIBD.CD.L3.OAerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.L3,o2_data_long.UNC.colon.nonIBD.CD.L3$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.nonIBD.CD.L3.Aerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.L3,o2_data_long.UNC.colon.nonIBD.CD.L3$classification=="Aerobic")
o2_data_long.UNC.colon.nonIBD.CD.L3.FAerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.L3,o2_data_long.UNC.colon.nonIBD.CD.L3$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.nonIBD.CD.L3.Anaerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.L3,o2_data_long.UNC.colon.nonIBD.CD.L3$classification=="Anaerobic")
o2_data_long.UNC.colon.nonIBD.CD.L3.OAnaerobic <- subset(o2_data_long.UNC.colon.nonIBD.CD.L3,o2_data_long.UNC.colon.nonIBD.CD.L3$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Disease,o2_data_long.UNC.colon.nonIBD.CD.L3.OAerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.nonIBD.CD.L3.Aerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.nonIBD.CD.L3.FAerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.nonIBD.CD.L3.Anaerobic)
wilcox.test(fraction~Disease,o2_data_long.UNC.colon.nonIBD.CD.L3.OAnaerobic)

#### CD L3 - Tissue
o2_data_long.UNC.CD.L3 <- subset(o2_data_long.UNC.CD, o2_data_long.UNC.CD$Location == "L3")

o2_data_long.UNC.CD.L3.OAerobic <- subset(o2_data_long.UNC.CD.L3,o2_data_long.UNC.CD.L3$classification=="Obligate.Aerobic")
o2_data_long.UNC.CD.L3.Aerobic <- subset(o2_data_long.UNC.CD.L3,o2_data_long.UNC.CD.L3$classification=="Aerobic")
o2_data_long.UNC.CD.L3.FAerobic <- subset(o2_data_long.UNC.CD.L3,o2_data_long.UNC.CD.L3$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.CD.L3.Anaerobic <- subset(o2_data_long.UNC.CD.L3,o2_data_long.UNC.CD.L3$classification=="Anaerobic")
o2_data_long.UNC.CD.L3.OAnaerobic <- subset(o2_data_long.UNC.CD.L3,o2_data_long.UNC.CD.L3$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.L3.OAerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.L3.Aerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.L3.FAerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.L3.Anaerobic)
wilcox.test(fraction~Tissue,o2_data_long.UNC.CD.L3.OAnaerobic)

#### MEDICATION USAGE - COLON ####
#### Colon - 5-ASA
o2_data_long.UNC.colon.CD.asa <- subset(o2_data_long.UNC.colon.CD,o2_data_long.UNC.colon.CD$X5_ASA %in% c("0","1"))

o2_data_long.UNC.colon.CD.asa.OAerobic <- subset(o2_data_long.UNC.colon.CD.asa,o2_data_long.UNC.colon.CD.asa$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.CD.asa.Aerobic <- subset(o2_data_long.UNC.colon.CD.asa,o2_data_long.UNC.colon.CD.asa$classification=="Aerobic")
o2_data_long.UNC.colon.CD.asa.FAerobic <- subset(o2_data_long.UNC.colon.CD.asa,o2_data_long.UNC.colon.CD.asa$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.CD.asa.Anaerobic <- subset(o2_data_long.UNC.colon.CD.asa,o2_data_long.UNC.colon.CD.asa$classification=="Anaerobic")
o2_data_long.UNC.colon.CD.asa.OAnaerobic <- subset(o2_data_long.UNC.colon.CD.asa,o2_data_long.UNC.colon.CD.asa$classification=="Obligate.Anaerobic")

wilcox.test(fraction~X5_ASA,o2_data_long.UNC.colon.CD.asa.OAerobic)
wilcox.test(fraction~X5_ASA,o2_data_long.UNC.colon.CD.asa.Aerobic)
wilcox.test(fraction~X5_ASA,o2_data_long.UNC.colon.CD.asa.FAerobic)
wilcox.test(fraction~X5_ASA,o2_data_long.UNC.colon.CD.asa.Anaerobic)
wilcox.test(fraction~X5_ASA,o2_data_long.UNC.colon.CD.asa.OAnaerobic)

#### Colon - Immunomodulator
o2_data_long.UNC.colon.CD.imm <- subset(o2_data_long.UNC.colon.CD,o2_data_long.UNC.colon.CD$Immunomodulators %in% c("0","1"))

o2_data_long.UNC.colon.CD.imm.OAerobic <- subset(o2_data_long.UNC.colon.CD.imm,o2_data_long.UNC.colon.CD.imm$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.CD.imm.Aerobic <- subset(o2_data_long.UNC.colon.CD.imm,o2_data_long.UNC.colon.CD.imm$classification=="Aerobic")
o2_data_long.UNC.colon.CD.imm.FAerobic <- subset(o2_data_long.UNC.colon.CD.imm,o2_data_long.UNC.colon.CD.imm$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.CD.imm.Anaerobic <- subset(o2_data_long.UNC.colon.CD.imm,o2_data_long.UNC.colon.CD.imm$classification=="Anaerobic")
o2_data_long.UNC.colon.CD.imm.OAnaerobic <- subset(o2_data_long.UNC.colon.CD.imm,o2_data_long.UNC.colon.CD.imm$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Immunomodulators,o2_data_long.UNC.colon.CD.imm.OAerobic)
wilcox.test(fraction~Immunomodulators,o2_data_long.UNC.colon.CD.imm.Aerobic)
wilcox.test(fraction~Immunomodulators,o2_data_long.UNC.colon.CD.imm.FAerobic)
wilcox.test(fraction~Immunomodulators,o2_data_long.UNC.colon.CD.imm.Anaerobic)
wilcox.test(fraction~Immunomodulators,o2_data_long.UNC.colon.CD.imm.OAnaerobic)

#### Colon - Steroids
o2_data_long.UNC.colon.CD.str <- subset(o2_data_long.UNC.colon.CD,o2_data_long.UNC.colon.CD$Steroids %in% c("0","1"))

o2_data_long.UNC.colon.CD.str.OAerobic <- subset(o2_data_long.UNC.colon.CD.str,o2_data_long.UNC.colon.CD.str$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.CD.str.Aerobic <- subset(o2_data_long.UNC.colon.CD.str,o2_data_long.UNC.colon.CD.str$classification=="Aerobic")
o2_data_long.UNC.colon.CD.str.FAerobic <- subset(o2_data_long.UNC.colon.CD.str,o2_data_long.UNC.colon.CD.str$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.CD.str.Anaerobic <- subset(o2_data_long.UNC.colon.CD.str,o2_data_long.UNC.colon.CD.str$classification=="Anaerobic")
o2_data_long.UNC.colon.CD.str.OAnaerobic <- subset(o2_data_long.UNC.colon.CD.str,o2_data_long.UNC.colon.CD.str$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Steroids,o2_data_long.UNC.colon.CD.str.OAerobic)
wilcox.test(fraction~Steroids,o2_data_long.UNC.colon.CD.str.Aerobic)
wilcox.test(fraction~Steroids,o2_data_long.UNC.colon.CD.str.FAerobic)
wilcox.test(fraction~Steroids,o2_data_long.UNC.colon.CD.str.Anaerobic)
wilcox.test(fraction~Steroids,o2_data_long.UNC.colon.CD.str.OAnaerobic)

#### Colon - anti-TNF
o2_data_long.UNC.colon.CD.tnf <- subset(o2_data_long.UNC.colon.CD,o2_data_long.UNC.colon.CD$TNF %in% c("0","1"))

o2_data_long.UNC.colon.CD.tnf.OAerobic <- subset(o2_data_long.UNC.colon.CD.tnf,o2_data_long.UNC.colon.CD.tnf$classification=="Obligate.Aerobic")
o2_data_long.UNC.colon.CD.tnf.Aerobic <- subset(o2_data_long.UNC.colon.CD.tnf,o2_data_long.UNC.colon.CD.tnf$classification=="Aerobic")
o2_data_long.UNC.colon.CD.tnf.FAerobic <- subset(o2_data_long.UNC.colon.CD.tnf,o2_data_long.UNC.colon.CD.tnf$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.colon.CD.tnf.Anaerobic <- subset(o2_data_long.UNC.colon.CD.tnf,o2_data_long.UNC.colon.CD.tnf$classification=="Anaerobic")
o2_data_long.UNC.colon.CD.tnf.OAnaerobic <- subset(o2_data_long.UNC.colon.CD.tnf,o2_data_long.UNC.colon.CD.tnf$classification=="Obligate.Anaerobic")

wilcox.test(fraction~TNF,o2_data_long.UNC.colon.CD.tnf.OAerobic)
wilcox.test(fraction~TNF,o2_data_long.UNC.colon.CD.tnf.Aerobic)
wilcox.test(fraction~TNF,o2_data_long.UNC.colon.CD.tnf.FAerobic)
wilcox.test(fraction~TNF,o2_data_long.UNC.colon.CD.tnf.Anaerobic)
wilcox.test(fraction~TNF,o2_data_long.UNC.colon.CD.tnf.OAnaerobic)

#### MEDICATION USAGE - ILEUM ####
#### ileum - 5-ASA
o2_data_long.UNC.ileum.CD.asa <- subset(o2_data_long.UNC.ileum.CD,o2_data_long.UNC.ileum.CD$X5_ASA %in% c("0","1"))

o2_data_long.UNC.ileum.CD.asa.OAerobic <- subset(o2_data_long.UNC.ileum.CD.asa,o2_data_long.UNC.ileum.CD.asa$classification=="Obligate.Aerobic")
o2_data_long.UNC.ileum.CD.asa.Aerobic <- subset(o2_data_long.UNC.ileum.CD.asa,o2_data_long.UNC.ileum.CD.asa$classification=="Aerobic")
o2_data_long.UNC.ileum.CD.asa.FAerobic <- subset(o2_data_long.UNC.ileum.CD.asa,o2_data_long.UNC.ileum.CD.asa$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.ileum.CD.asa.Anaerobic <- subset(o2_data_long.UNC.ileum.CD.asa,o2_data_long.UNC.ileum.CD.asa$classification=="Anaerobic")
o2_data_long.UNC.ileum.CD.asa.OAnaerobic <- subset(o2_data_long.UNC.ileum.CD.asa,o2_data_long.UNC.ileum.CD.asa$classification=="Obligate.Anaerobic")

wilcox.test(fraction~X5_ASA,o2_data_long.UNC.ileum.CD.asa.OAerobic)
wilcox.test(fraction~X5_ASA,o2_data_long.UNC.ileum.CD.asa.Aerobic)
wilcox.test(fraction~X5_ASA,o2_data_long.UNC.ileum.CD.asa.FAerobic)
wilcox.test(fraction~X5_ASA,o2_data_long.UNC.ileum.CD.asa.Anaerobic)
wilcox.test(fraction~X5_ASA,o2_data_long.UNC.ileum.CD.asa.OAnaerobic)

#### Ileum - Immunomodulator
o2_data_long.UNC.ileum.CD.imm <- subset(o2_data_long.UNC.ileum.CD,o2_data_long.UNC.ileum.CD$Immunomodulators %in% c("0","1"))

o2_data_long.UNC.ileum.CD.imm.OAerobic <- subset(o2_data_long.UNC.ileum.CD.imm,o2_data_long.UNC.ileum.CD.imm$classification=="Obligate.Aerobic")
o2_data_long.UNC.ileum.CD.imm.Aerobic <- subset(o2_data_long.UNC.ileum.CD.imm,o2_data_long.UNC.ileum.CD.imm$classification=="Aerobic")
o2_data_long.UNC.ileum.CD.imm.FAerobic <- subset(o2_data_long.UNC.ileum.CD.imm,o2_data_long.UNC.ileum.CD.imm$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.ileum.CD.imm.Anaerobic <- subset(o2_data_long.UNC.ileum.CD.imm,o2_data_long.UNC.ileum.CD.imm$classification=="Anaerobic")
o2_data_long.UNC.ileum.CD.imm.OAnaerobic <- subset(o2_data_long.UNC.ileum.CD.imm,o2_data_long.UNC.ileum.CD.imm$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Immunomodulators,o2_data_long.UNC.ileum.CD.imm.OAerobic)
wilcox.test(fraction~Immunomodulators,o2_data_long.UNC.ileum.CD.imm.Aerobic)
wilcox.test(fraction~Immunomodulators,o2_data_long.UNC.ileum.CD.imm.FAerobic)
wilcox.test(fraction~Immunomodulators,o2_data_long.UNC.ileum.CD.imm.Anaerobic)
wilcox.test(fraction~Immunomodulators,o2_data_long.UNC.ileum.CD.imm.OAnaerobic)

#### Ileum - Steroids
o2_data_long.UNC.ileum.CD.str <- subset(o2_data_long.UNC.ileum.CD,o2_data_long.UNC.ileum.CD$Steroids %in% c("0","1"))

o2_data_long.UNC.ileum.CD.str.OAerobic <- subset(o2_data_long.UNC.ileum.CD.str,o2_data_long.UNC.ileum.CD.str$classification=="Obligate.Aerobic")
o2_data_long.UNC.ileum.CD.str.Aerobic <- subset(o2_data_long.UNC.ileum.CD.str,o2_data_long.UNC.ileum.CD.str$classification=="Aerobic")
o2_data_long.UNC.ileum.CD.str.FAerobic <- subset(o2_data_long.UNC.ileum.CD.str,o2_data_long.UNC.ileum.CD.str$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.ileum.CD.str.Anaerobic <- subset(o2_data_long.UNC.ileum.CD.str,o2_data_long.UNC.ileum.CD.str$classification=="Anaerobic")
o2_data_long.UNC.ileum.CD.str.OAnaerobic <- subset(o2_data_long.UNC.ileum.CD.str,o2_data_long.UNC.ileum.CD.str$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Steroids,o2_data_long.UNC.ileum.CD.str.OAerobic)
wilcox.test(fraction~Steroids,o2_data_long.UNC.ileum.CD.str.Aerobic)
wilcox.test(fraction~Steroids,o2_data_long.UNC.ileum.CD.str.FAerobic)
wilcox.test(fraction~Steroids,o2_data_long.UNC.ileum.CD.str.Anaerobic)
wilcox.test(fraction~Steroids,o2_data_long.UNC.ileum.CD.str.OAnaerobic)

#### Ileum - anti-TNF
o2_data_long.UNC.ileum.CD.tnf <- subset(o2_data_long.UNC.ileum.CD,o2_data_long.UNC.ileum.CD$TNF %in% c("0","1"))

o2_data_long.UNC.ileum.CD.tnf.OAerobic <- subset(o2_data_long.UNC.ileum.CD.tnf,o2_data_long.UNC.ileum.CD.tnf$classification=="Obligate.Aerobic")
o2_data_long.UNC.ileum.CD.tnf.Aerobic <- subset(o2_data_long.UNC.ileum.CD.tnf,o2_data_long.UNC.ileum.CD.tnf$classification=="Aerobic")
o2_data_long.UNC.ileum.CD.tnf.FAerobic <- subset(o2_data_long.UNC.ileum.CD.tnf,o2_data_long.UNC.ileum.CD.tnf$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.ileum.CD.tnf.Anaerobic <- subset(o2_data_long.UNC.ileum.CD.tnf,o2_data_long.UNC.ileum.CD.tnf$classification=="Anaerobic")
o2_data_long.UNC.ileum.CD.tnf.OAnaerobic <- subset(o2_data_long.UNC.ileum.CD.tnf,o2_data_long.UNC.ileum.CD.tnf$classification=="Obligate.Anaerobic")

wilcox.test(fraction~TNF,o2_data_long.UNC.ileum.CD.tnf.OAerobic)
wilcox.test(fraction~TNF,o2_data_long.UNC.ileum.CD.tnf.Aerobic)
wilcox.test(fraction~TNF,o2_data_long.UNC.ileum.CD.tnf.FAerobic)
wilcox.test(fraction~TNF,o2_data_long.UNC.ileum.CD.tnf.Anaerobic)
wilcox.test(fraction~TNF,o2_data_long.UNC.ileum.CD.tnf.OAnaerobic)





#### Recurrence - Ileum ####
o2_data_long.UNC.ileum.CD.recur.OAerobic <- subset(o2_data_long.UNC.ileum.CD.recur,o2_data_long.UNC.ileum.CD.recur$classification=="Obligate.Aerobic")
o2_data_long.UNC.ileum.CD.recur.Aerobic <- subset(o2_data_long.UNC.ileum.CD.recur,o2_data_long.UNC.ileum.CD.recur$classification=="Aerobic")
o2_data_long.UNC.ileum.CD.recur.FAerobic <- subset(o2_data_long.UNC.ileum.CD.recur,o2_data_long.UNC.ileum.CD.recur$classification=="Facultative.Anaerobic...Microaerophilic")
o2_data_long.UNC.ileum.CD.recur.Anaerobic <- subset(o2_data_long.UNC.ileum.CD.recur,o2_data_long.UNC.ileum.CD.recur$classification=="Anaerobic")
o2_data_long.UNC.ileum.CD.recur.OAnaerobic <- subset(o2_data_long.UNC.ileum.CD.recur,o2_data_long.UNC.ileum.CD.recur$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Recurrence,o2_data_long.UNC.ileum.CD.recur.OAerobic) 
wilcox.test(fraction~Recurrence,o2_data_long.UNC.ileum.CD.recur.Aerobic)
wilcox.test(fraction~Recurrence,o2_data_long.UNC.ileum.CD.recur.FAerobic)
wilcox.test(fraction~Recurrence,o2_data_long.UNC.ileum.CD.recur.Anaerobic)
wilcox.test(fraction~Recurrence,o2_data_long.UNC.ileum.CD.recur.OAnaerobic) 


#### FIGURES ####
theme_classic()
# nonIBD
o2_data_long.UNC.nonIBD$classification_f <- factor(o2_data_long.UNC.nonIBD$classification, levels=c('Obligate.Anaerobic','Anaerobic','Facultative.Anaerobic...Microaerophilic','Aerobic','Obligate.Aerobic'))
o2.UNC.nonIBD <- ggplot(o2_data_long.UNC.nonIBD, aes(x=Tissue, y=fraction, group=DiseaseTissue)) + 
  geom_boxplot(aes(fill=DiseaseTissue),outlier.size=3) + facet_wrap(~ classification_f, ncol=3) +
  stat_compare_means(aes(group = DiseaseTissue),label = "p.signif", label.x.npc = "center", method = "wilcox.test") +
  scale_fill_manual(values=c('#00b0f0','#8b6cc9')) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.x = unit(0,"line"))
o2.UNC.nonIBD+theme_classic()

# CD
o2_data_long.UNC.CD$classification_f <- factor(o2_data_long.UNC.CD$classification, levels=c('Obligate.Anaerobic','Anaerobic','Facultative.Anaerobic...Microaerophilic','Aerobic','Obligate.Aerobic'))
o2.UNC.CD <- ggplot(o2_data_long.UNC.CD , aes(x=DiseaseTissue, y=fraction, group=DiseaseTissue)) + 
  geom_boxplot(aes(fill=DiseaseTissue),outlier.size=3) + facet_wrap(~ classification_f, ncol=3) +
  stat_compare_means(aes(group = DiseaseTissue),label = "p.signif", label.x.npc = "center", method = "wilcox.test") +
  scale_fill_manual(values=c('#cb5998','#cc5b48')) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.x = unit(0,"line"))
o2.UNC.CD+theme_classic()

# Colon
o2_data_long.UNC.colon$classification_f <- factor(o2_data_long.UNC.colon$classification, levels=c('Obligate.Anaerobic','Anaerobic','Facultative.Anaerobic...Microaerophilic','Aerobic','Obligate.Aerobic'))

o2.UNC.colon <- ggplot(o2_data_long.UNC.colon, aes(x=Disease, y=fraction, group=DiseaseTissue)) + 
  geom_boxplot(aes(fill=DiseaseTissue),outlier.size=3) + facet_wrap(~ classification_f, ncol=3)+
  stat_compare_means(aes(group = DiseaseTissue),label = "p.signif",label.x.npc = "center", method = "wilcox.test") +
  scale_fill_manual(values=c('#cb5998','#00b0f0')) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.x = unit(0,"line")) 
o2.UNC.colon+theme_classic()

# Ileum 
o2_data_long.UNC.ileum$classification_f <- factor(o2_data_long.UNC.ileum$classification, levels=c('Obligate.Anaerobic','Anaerobic','Facultative.Anaerobic...Microaerophilic','Aerobic','Obligate.Aerobic'))
o2.UNC.ileum <- ggplot(o2_data_long.UNC.ileum, aes(x=Disease, y=fraction, group=DiseaseTissue)) + 
  geom_boxplot(aes(fill=DiseaseTissue),outlier.size=3) + facet_wrap(~ classification_f, ncol=3) +
  stat_compare_means(aes(group = DiseaseTissue),label = "p.signif", label.x.npc = "center", method = "wilcox.test") +
  scale_fill_manual(values=c('#cc5b48','#8b6cc9')) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.spacing.x = unit(0,"line"))
o2.UNC.ileum+theme_classic()

aero_figure <- ggarrange(o2.UNC.nonIBD + rremove("x.text"),o2.UNC.colon + rremove("x.text") , o2.UNC.ileum + rremove("x.text"),o2.UNC.CD + rremove("x.text"), 
                         labels = c("A", "B", "C","D"),
                         ncol = 2, nrow = 2)
annotate_figure(aero_figure,
                bottom = text_grob("Oxygen Shifts in Intestinal Mucosa in CD", color = "black", face = "bold", size = 14),
                fig.lab = "Figure 2", fig.lab.face = "bold",fig.lab.pos = "bottom.left"
)

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
o2.UNC.colon.OAerobic <- ggplot(o2_data_long.UNC.colon.OAerobic , aes(x=DiseaseNo, y=fraction,color=Disease)) + 
  geom_point(size=5) +
  labs(y = "Relative Abundance")+
  scale_color_manual(values=c('#cb5998','#00b0f0'))+
  coord_cartesian(ylim= c(0,1))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

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
#### WashU ####
washu_unc_data_long <- gather(UNC_WashU_CD_meta_aero_table_no_anti_ileum_only_20190720, classification,fraction, Aerobic:Obligate.Anaerobic, factor_key=TRUE)
washu_unc_data_long.NI <- subset(washu_unc_data_long, washu_unc_data_long$Pathology =="NI")

washu_unc_data_long <- washu_unc_data_long.NI

washu_unc_data_long$DiseaseTissue <- paste(washu_unc_data_long$Disease,washu_unc_data_long$Tissue,sep=" ")

washu_unc_data_recur_only <- subset(washu_unc_data_long,washu_unc_data_long$Recurrence=="recurrence")
washu_unc_data_no_recur_only <- subset(washu_unc_data_long,washu_unc_data_long$Recurrence=="no recurrence")

washu_unc_data_long_recurrence <- subset(washu_unc_data_long,washu_unc_data_long$Recurrence %in% c("recurrence","no recurrence"))

washu_recurrence <- subset(washu_unc_data_long_recurrence,washu_unc_data_long_recurrence$Cohort=="WashU")
unc_recurrence <- subset(washu_unc_data_long_recurrence,washu_unc_data_long_recurrence$Cohort=="UNC")

####Compare UNC vs WashU ####
washu_unc_data_long.OAerobic <- subset(washu_unc_data_long,washu_unc_data_long$classification=="Obligate.Aerobic")
washu_unc_data_long.Aerobic <- subset(washu_unc_data_long,washu_unc_data_long$classification=="Aerobic")
washu_unc_data_long.FAerobic <- subset(washu_unc_data_long,washu_unc_data_long$classification=="Facultative.Anaerobic...Microaerophilic")
washu_unc_data_long.Anaerobic <- subset(washu_unc_data_long,washu_unc_data_long$classification=="Anaerobic")
washu_unc_data_long.OAnaerobic <- subset(washu_unc_data_long,washu_unc_data_long$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Cohort,washu_unc_data_long.OAerobic) 
wilcox.test(fraction~Cohort,washu_unc_data_long.Aerobic)
wilcox.test(fraction~Cohort,washu_unc_data_long.FAerobic)
wilcox.test(fraction~Cohort,washu_unc_data_long.Anaerobic)
wilcox.test(fraction~Cohort,washu_unc_data_long.OAnaerobic) 

####Compare UNC vs WashU : recurrence ####
washu_unc_data_recur_only.OAerobic <- subset(washu_unc_data_recur_only,washu_unc_data_recur_only$classification=="Obligate.Aerobic")
washu_unc_data_recur_only.Aerobic <- subset(washu_unc_data_recur_only,washu_unc_data_recur_only$classification=="Aerobic")
washu_unc_data_recur_only.FAerobic <- subset(washu_unc_data_recur_only,washu_unc_data_recur_only$classification=="Facultative.Anaerobic...Microaerophilic")
washu_unc_data_recur_only.Anaerobic <- subset(washu_unc_data_recur_only,washu_unc_data_recur_only$classification=="Anaerobic")
washu_unc_data_recur_only.OAnaerobic <- subset(washu_unc_data_recur_only,washu_unc_data_recur_only$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Cohort,washu_unc_data_recur_only.OAerobic) 
wilcox.test(fraction~Cohort,washu_unc_data_recur_only.Aerobic)
wilcox.test(fraction~Cohort,washu_unc_data_recur_only.FAerobic)
wilcox.test(fraction~Cohort,washu_unc_data_recur_only.Anaerobic)
wilcox.test(fraction~Cohort,washu_unc_data_recur_only.OAnaerobic) 

####Compare UNC vs WashU : no recurrence ####
washu_unc_data_no_recur_only.OAerobic <- subset(washu_unc_data_no_recur_only,washu_unc_data_no_recur_only$classification=="Obligate.Aerobic")
washu_unc_data_no_recur_only.Aerobic <- subset(washu_unc_data_no_recur_only,washu_unc_data_no_recur_only$classification=="Aerobic")
washu_unc_data_no_recur_only.FAerobic <- subset(washu_unc_data_no_recur_only,washu_unc_data_no_recur_only$classification=="Facultative.Anaerobic...Microaerophilic")
washu_unc_data_no_recur_only.Anaerobic <- subset(washu_unc_data_no_recur_only,washu_unc_data_no_recur_only$classification=="Anaerobic")
washu_unc_data_no_recur_only.OAnaerobic <- subset(washu_unc_data_no_recur_only,washu_unc_data_no_recur_only$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Cohort,washu_unc_data_no_recur_only.OAerobic) 
wilcox.test(fraction~Cohort,washu_unc_data_no_recur_only.Aerobic)
wilcox.test(fraction~Cohort,washu_unc_data_no_recur_only.FAerobic)
wilcox.test(fraction~Cohort,washu_unc_data_no_recur_only.Anaerobic)
wilcox.test(fraction~Cohort,washu_unc_data_no_recur_only.OAnaerobic) 

####WashU - recurrence ####
washu_recurrence.OAerobic <- subset(washu_recurrence ,washu_recurrence$classification=="Obligate.Aerobic")
washu_recurrence.Aerobic <- subset(washu_recurrence ,washu_recurrence$classification=="Aerobic")
washu_recurrence.FAerobic <- subset(washu_recurrence ,washu_recurrence$classification=="Facultative.Anaerobic...Microaerophilic")
washu_recurrence.Anaerobic <- subset(washu_recurrence ,washu_recurrence$classification=="Anaerobic")
washu_recurrence.OAnaerobic <- subset(washu_recurrence, washu_recurrence$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Recurrence,washu_recurrence.OAerobic) 
wilcox.test(fraction~Recurrence,washu_recurrence.Aerobic)
wilcox.test(fraction~Recurrence,washu_recurrence.FAerobic)
wilcox.test(fraction~Recurrence,washu_recurrence.Anaerobic)
wilcox.test(fraction~Recurrence,washu_recurrence.OAnaerobic) 

####UNC - recurrence ####
unc_recurrence.OAerobic <- subset(unc_recurrence ,unc_recurrence$classification=="Obligate.Aerobic")
unc_recurrence.Aerobic <- subset(unc_recurrence ,unc_recurrence$classification=="Aerobic")
unc_recurrence.FAerobic <- subset(unc_recurrence ,unc_recurrence$classification=="Facultative.Anaerobic...Microaerophilic")
unc_recurrence.Anaerobic <- subset(unc_recurrence ,unc_recurrence$classification=="Anaerobic")
unc_recurrence.OAnaerobic <- subset(unc_recurrence, unc_recurrence$classification=="Obligate.Anaerobic")

wilcox.test(fraction~Recurrence,unc_recurrence.OAerobic) 
wilcox.test(fraction~Recurrence,unc_recurrence.Aerobic)
wilcox.test(fraction~Recurrence,unc_recurrence.FAerobic)
wilcox.test(fraction~Recurrence,unc_recurrence.Anaerobic)
wilcox.test(fraction~Recurrence,unc_recurrence.OAnaerobic) 
