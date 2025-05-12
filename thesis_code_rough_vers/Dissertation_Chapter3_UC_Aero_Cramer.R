##############################
# Aerobic/Anaerobe Analysis
# Cramer Test
#
##############################
library(cramer)
library(MASS)
patient_tissue_o2_summary <- UNC_patients_meta_no_anti_aero_table_relevant_20200309

UNC.NI <- subset(patient_tissue_o2_summary,patient_tissue_o2_summary$Pathology=="NI")

UNC.nonIBD <- subset(UNC.NI,UNC.NI$Disease == "nonIBD")
UNC.nonIBD.colon <- subset(UNC.nonIBD,UNC.nonIBD$Tissue=="colon")

UNC.CD <- subset(UNC.NI,UNC.NI$Disease=="CD")
UNC.CD.colon <- subset(UNC.CD,UNC.CD$Tissue=="colon")

UNC.UC <- subset(UNC.NI,UNC.NI$Disease=="UC")
UNC.UC.colon <- subset(UNC.UC,UNC.UC$Tissue=="colon")

UNC.UC.colon$AeroLvl <- ifelse(UNC.UC.colon$PatientNo %in% c("127","41500","42700","43000","46","52"),"Low","Health")

UNC.UC.colon.health <- subset(UNC.UC.colon,UNC.UC.colon$AeroLvl=="Health")
UNC.UC.colon.low <- subset(UNC.UC.colon,UNC.UC.colon$AeroLvl=="Low")
### Get Aerotolereance Categories
UNC.nonIBD.colon.aero <- subset(UNC.nonIBD.colon,select=c("Aerobic","Anaerobic","Facultative.Anaerobic...Microaerophilic","Obligate.Aerobic","Obligate.Anaerobic"))

UNC.CD.colon.aero <- subset(UNC.CD.colon,select=c("Aerobic","Anaerobic","Facultative.Anaerobic...Microaerophilic","Obligate.Aerobic","Obligate.Anaerobic"))

UNC.UC.colon.aero <- subset(UNC.UC.colon,select=c("Aerobic","Anaerobic","Facultative.Anaerobic...Microaerophilic","Obligate.Aerobic","Obligate.Anaerobic"))

UNC.UC.colon.h.aero <- subset(UNC.UC.colon.health,select=c("Aerobic","Anaerobic","Facultative.Anaerobic...Microaerophilic","Obligate.Aerobic","Obligate.Anaerobic"))
UNC.UC.colon.l.aero <- subset(UNC.UC.colon.low,select=c("Aerobic","Anaerobic","Facultative.Anaerobic...Microaerophilic","Obligate.Aerobic","Obligate.Anaerobic"))

#Remove row names
rownames(UNC.nonIBD.colon.aero) <- c()
rownames(UNC.CD.colon.aero) <- c()
rownames(UNC.UC.colon.aero) <- c()

rownames(UNC.UC.colon.h.aero) <- c()
rownames(UNC.UC.colon.l.aero) <- c()

#Remove column names
colnames(UNC.nonIBD.colon.aero) <- c()
colnames(UNC.CD.colon.aero) <- c()
colnames(UNC.UC.colon.aero) <- c()

colnames(UNC.UC.colon.h.aero) <- c()
colnames(UNC.UC.colon.l.aero) <- c()

#Convert into a data matrix
UNC.CD.colon.aero.mat <- as.matrix(UNC.CD.colon.aero)
UNC.nonIBD.colon.aero.mat <- as.matrix(UNC.nonIBD.colon.aero)
UNC.UC.colon.aero.mat <- as.matrix(UNC.UC.colon.aero)

UNC.UC.colon.H.aero.mat <- as.matrix(UNC.UC.colon.h.aero)
UNC.UC.colon.L.aero.mat <- as.matrix(UNC.UC.colon.l.aero)

#Tests!!

#Colon disease
cramer.test(UNC.nonIBD.colon.aero.mat,UNC.CD.colon.aero.mat)
cramer.test(UNC.nonIBD.colon.aero.mat,UNC.UC.colon.aero.mat)
cramer.test(UNC.UC.colon.aero.mat,UNC.CD.colon.aero.mat)

cramer.test(UNC.UC.colon.H.aero.mat,UNC.UC.colon.L.aero.mat)

cramer.test(UNC.UC.colon.H.aero.mat,UNC.CD.colon.aero.mat)
cramer.test(UNC.UC.colon.L.aero.mat,UNC.CD.colon.aero.mat)

cramer.test(UNC.UC.colon.H.aero.mat,UNC.nonIBD.colon.aero.mat)
cramer.test(UNC.UC.colon.L.aero.mat,UNC.nonIBD.colon.aero.mat)
