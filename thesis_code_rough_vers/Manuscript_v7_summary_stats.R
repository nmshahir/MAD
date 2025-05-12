####Manuscript - Summary Statistics####

#### Data Clean Bit for UNC Only####
ps.patients.g <- readRDS("ps_patients_genus_20190713.rds")
sample_data(ps.patients.g)$Age <- as.numeric(levels(sample_data(ps.patients.g)$Age))[sample_data(ps.patients.g)$Age]

#Get all samples above 8000 reads
ps.patients.g.10K <- prune_samples(sample_sums(ps.patients.g)>8000,ps.patients.g)
ps.patients.g.10K <- prune_taxa(taxa_sums(ps.patients.g.10K)>0,ps.patients.g.10K)

ps.temp <- ps.patients.g.10K

#Get all samples that are either Ascending, Terminal, Transverse/Left, or unknown
ps.temp.2 <- subset_samples(ps.temp,!(sample_data(ps.temp)$Region %in% c("Cecal","Descending","fecal","Rectal/Sigmoid","Sigmoid","Transverse/Left")))
ps.temp.2 <- prune_taxa(taxa_sums(ps.temp.2)>0,ps.temp.2)

ps.patients.g <- ps.temp.2

#Get Only UNC samples
ps.UNC.patients <- subset_samples(ps.patients.g,sample_data(ps.patients.g)$Cohort == "UNC")
ps.UNC.patients <- prune_taxa(taxa_sums(ps.UNC.patients)>0,ps.UNC.patients)

#Remove All Samples with Antibiotics
ps.UNC.patients.no.anti <- subset_samples(ps.UNC.patients,sample_data(ps.UNC.patients)$Antibiotics == 0)
ps.UNC.patients.no.anti <- prune_taxa(taxa_sums(ps.UNC.patients.no.anti)>0,ps.UNC.patients.no.anti)

#704 RSVs seen specifically in antibioitics in patients with antibiotics
ps.UNC.patients <- ps.UNC.patients.no.anti

# Remove patient 25
ps.UNC.patients.no.25 <- subset_samples(ps.UNC.patients,sample_data(ps.UNC.patients)$PatientNo != "25")
ps.UNC.patients.no.25 <- prune_taxa(taxa_sums(ps.UNC.patients.no.25)>0,ps.UNC.patients.no.25)

ps.UNC.patients <- ps.UNC.patients.no.25

ps.UNC.patients.NI.I <- subset_samples(ps.UNC.patients,sample_data(ps.UNC.patients)$Pathology %in% c("NI","I"))
ps.UNC.patients.NI.I <- prune_taxa(taxa_sums(ps.UNC.patients.NI.I)>0,ps.UNC.patients.NI.I)

ps.UNC.patients <- ps.UNC.patients.NI.I

ps.WashU.patients <- subset_samples(ps.patients.g,sample_data(ps.patients.g)$Cohort == "WashU")
ps.WashU.patients <- prune_taxa(taxa_sums(ps.WashU.patients)>0,ps.WashU.patients)

ps.WashU.patients.NI <- subset_samples(ps.WashU.patients,sample_data(ps.WashU.patients)$Pathology == "NI")
ps.WashU.patients.NI <- prune_taxa(taxa_sums(ps.WashU.patients.NI)>0,ps.WashU.patients.NI)

ps.WashU.patients <- ps.WashU.patients.NI

ps.WashU.no.anti <- subset_samples(ps.WashU.patients,sample_data(ps.WashU.patients)$Antibiotics == "0")
ps.WashU.no.anti <- prune_taxa(taxa_sums(ps.WashU.no.anti)>0,ps.WashU.no.anti)

ps.WashU.patients <- ps.WashU.no.anti

ps.L3 <- subset_samples(ps.UNC.patients,sample_data(ps.UNC.patients)$Location == "L3")
ps.L3 <- prune_taxa(taxa_sums(ps.L3)>0,ps.L3)

ps.L3.B2 <- subset_samples(ps.L3,sample_data(ps.L3)$Behavior== "B2")
ps.L3.B2 <- prune_taxa(taxa_sums(ps.L3.B2)>0,ps.L3.B2)

ps.L3.B2.colon <- subset_samples(ps.L3.B2,sample_data(ps.L3.B2)$AnatomSite == "colon")
ps.L3.B2.colon <- prune_taxa(taxa_sums(ps.L3.B2.colon)>0, ps.L3.B2.colon)
#### Number of Total UNC Patients 
#Merge Samples if any
ps.UNC <- merge_samples(ps.UNC.patients,"PatientNo")
sample_data(ps.UNC)$DiseaseStatus <- as.integer(sample_data(ps.UNC)$DiseaseStatus)
sample_data(ps.UNC)$DiseaseStatus <- as.factor(sample_data(ps.UNC)$DiseaseStatus)
levels(sample_data(ps.UNC)$DiseaseStatus) <- levels(sample_data(ps.UNC.patients)$DiseaseStatus)
#CD - 46
#nonIBD - 23

#Tissue
sample_data(ps.UNC.patients)$DiseaseTissue <- paste(sample_data(ps.UNC.patients)$PatientNo,sample_data(ps.UNC.patients)$Disease, sample_data(ps.UNC.patients)$AnatomSite,sep=" ")
ps.UNC <- merge_samples(ps.UNC.patients,"DiseaseTissue")
sample_data(ps.UNC)$DiseaseStatus <- as.integer(sample_data(ps.UNC)$DiseaseStatus)
sample_data(ps.UNC)$DiseaseStatus <- as.factor(sample_data(ps.UNC)$DiseaseStatus)
levels(sample_data(ps.UNC)$DiseaseStatus) <- levels(sample_data(ps.UNC.patients)$DiseaseStatus)

sample_data(ps.UNC)$AnatomSite <- as.integer(sample_data(ps.UNC)$AnatomSite)
sample_data(ps.UNC)$AnatomSite <- as.factor(sample_data(ps.UNC)$AnatomSite)
levels(sample_data(ps.UNC)$AnatomSite) <- levels(sample_data(ps.UNC.patients)$AnatomSite)

#Tissue + Inflammation
sample_data(ps.UNC.patients)$DiseaseTissue <- paste(sample_data(ps.UNC.patients)$PatientNo,sample_data(ps.UNC.patients)$Disease, sample_data(ps.UNC.patients)$AnatomSite,sample_data(ps.UNC.patients)$Pathology,sep=" ")
ps.UNC <- merge_samples(ps.UNC.patients,"DiseaseTissue")
sample_data(ps.UNC)$DiseaseStatus <- as.integer(sample_data(ps.UNC)$DiseaseStatus)
sample_data(ps.UNC)$DiseaseStatus <- as.factor(sample_data(ps.UNC)$DiseaseStatus)
levels(sample_data(ps.UNC)$DiseaseStatus) <- levels(sample_data(ps.UNC.patients)$DiseaseStatus)

sample_data(ps.UNC)$AnatomSite <- as.integer(sample_data(ps.UNC)$AnatomSite)
sample_data(ps.UNC)$AnatomSite <- as.factor(sample_data(ps.UNC)$AnatomSite)
levels(sample_data(ps.UNC)$AnatomSite) <- levels(sample_data(ps.UNC.patients)$AnatomSite)

sample_data(ps.UNC)$Pathology <- as.integer(sample_data(ps.UNC)$Pathology)
sample_data(ps.UNC)$Pathology<- as.factor(sample_data(ps.UNC)$Pathology)
levels(sample_data(ps.UNC)$Pathology) <- levels(sample_data(ps.UNC.patients)$Pathology)

sample_data(ps.UNC)$PatientNo <- as.integer(sample_data(ps.UNC)$PatientNo)
sample_data(ps.UNC)$PatientNo <- as.factor(sample_data(ps.UNC)$PatientNo)
levels(sample_data(ps.UNC)$PatientNo) <- levels(sample_data(ps.UNC.patients)$PatientNo)

sample_data(ps.UNC)$Behavior <- as.integer(sample_data(ps.UNC)$Behavior)
sample_data(ps.UNC)$Behavior   <- as.factor(sample_data(ps.UNC)$Behavior)
levels(sample_data(ps.UNC)$Behavior  ) <- levels(sample_data(ps.UNC.patients)$Behavior  )

sample_data(ps.UNC)$Location <- as.integer(sample_data(ps.UNC)$Location)
sample_data(ps.UNC)$Location   <- as.factor(sample_data(ps.UNC)$Location)
levels(sample_data(ps.UNC)$Location  ) <- levels(sample_data(ps.UNC.patients)$Location  )

sample_data(ps.UNC)$Sex <- as.integer(sample_data(ps.UNC)$Sex)
sample_data(ps.UNC)$Sex   <- as.factor(sample_data(ps.UNC)$Sex)
levels(sample_data(ps.UNC)$Sex) <- levels(sample_data(ps.UNC.patients)$Sex)

sample_data(ps.UNC)$Race <- as.integer(sample_data(ps.UNC)$Race)
sample_data(ps.UNC)$Race   <- as.factor(sample_data(ps.UNC)$Race)
levels(sample_data(ps.UNC)$Race) <- levels(sample_data(ps.UNC.patients)$Race)

sample_data(ps.UNC)$SmokingStatus <- as.integer(sample_data(ps.UNC)$SmokingStatus)
sample_data(ps.UNC)$SmokingStatus   <- as.factor(sample_data(ps.UNC)$SmokingStatus)
levels(sample_data(ps.UNC)$SmokingStatus) <- levels(as.factor(sample_data(ps.UNC.patients)$SmokingStatus))

sample_data(ps.UNC)$Antibiotics <- as.integer(sample_data(ps.UNC)$Antibiotics)
sample_data(ps.UNC)$Antibiotics   <- as.factor(sample_data(ps.UNC)$Antibiotics)
levels(sample_data(ps.UNC)$Antibiotics) <- levels(as.factor(sample_data(ps.UNC.patients)$Antibiotics))

sample_data(ps.UNC)$X5_ASA <- as.integer(sample_data(ps.UNC)$X5_ASA)
sample_data(ps.UNC)$X5_ASA   <- as.factor(sample_data(ps.UNC)$X5_ASA)
levels(sample_data(ps.UNC)$X5_ASA) <- levels(as.factor(sample_data(ps.UNC.patients)$X5_ASA))

sample_data(ps.UNC)$Steroids <- as.integer(sample_data(ps.UNC)$Steroids)
sample_data(ps.UNC)$Steroids   <- as.factor(sample_data(ps.UNC)$Steroids)
levels(sample_data(ps.UNC)$Steroids) <- levels(as.factor(sample_data(ps.UNC.patients)$Steroids))

sample_data(ps.UNC)$Immunomodulators <- as.integer(sample_data(ps.UNC)$Immunomodulators)
sample_data(ps.UNC)$Immunomodulators   <- as.factor(sample_data(ps.UNC)$Immunomodulators)
levels(sample_data(ps.UNC)$Immunomodulators) <- levels(as.factor(sample_data(ps.UNC.patients)$Immunomodulators))

sample_data(ps.UNC)$TNF <- as.integer(sample_data(ps.UNC)$TNF)
sample_data(ps.UNC)$TNF  <- as.factor(sample_data(ps.UNC)$TNF)
levels(sample_data(ps.UNC)$TNF) <- levels(as.factor(sample_data(ps.UNC.patients)$TNF))

sample_data(ps.UNC)$Cdiff <- as.integer(sample_data(ps.UNC)$Cdiff)
sample_data(ps.UNC)$Cdiff  <- as.factor(sample_data(ps.UNC)$Cdiff)
levels(sample_data(ps.UNC)$Cdiff) <- levels(as.factor(sample_data(ps.UNC.patients)$Cdiff))

sample_data(ps.UNC)$Ethnicity <- as.integer(sample_data(ps.UNC)$Ethnicity)
sample_data(ps.UNC)$Ethnicity  <- as.factor(sample_data(ps.UNC)$Ethnicity)
levels(sample_data(ps.UNC)$Ethnicity) <- levels(as.factor(sample_data(ps.UNC.patients)$Ethnicity))


sample_data(ps.UNC)$Perianal <- as.integer(sample_data(ps.UNC)$Perianal)
sample_data(ps.UNC)$Perianal <- as.factor(sample_data(ps.UNC)$Perianal)
levels(sample_data(ps.UNC)$Perianal) <- levels(as.factor(sample_data(ps.UNC.patients)$Perianal))

sample_data(ps.UNC)$Subtype <- as.integer(sample_data(ps.UNC)$Subtype)
sample_data(ps.UNC)$Subtype <- as.factor(sample_data(ps.UNC)$Subtype)
levels(sample_data(ps.UNC)$Subtype) <- levels(as.factor(sample_data(ps.UNC.patients)$Subtype))

sample_data(ps.UNC)$Rutgeerts_Score <- as.integer(sample_data(ps.UNC)$Rutgeerts_Score )
sample_data(ps.UNC)$Rutgeerts_Score <- as.factor(sample_data(ps.UNC)$Rutgeerts_Score )
levels(sample_data(ps.UNC)$Rutgeerts_Score) <- levels(as.factor(sample_data(ps.UNC.patients)$Rutgeerts_Score))

sample_data(ps.UNC)$Recurrence <- ifelse(sample_data(ps.UNC)$Rutgeerts_Score %in% c("i0","i1","i0-i1","i1-i2"),"no recurrence",
                                              ifelse(sample_data(ps.UNC)$Rutgeerts_Score %in%c("i2-i3","i3","i4"),"recurrence"," ")
)

sample_data(ps.UNC)$Age <- as.integer(sample_data(ps.UNC)$Age)
ps.UNC.CD.nonIBD <- subset_samples(ps.UNC,sample_data(ps.UNC)$DiseaseStatus %in% c("CD","nonIBD"))
ps.UNC.CD.nonIBD <- prune_taxa(taxa_sums(ps.UNC.CD.nonIBD)>0,ps.UNC.CD.nonIBD)


ps.UNC.NI <- subset_samples(ps.UNC.CD.nonIBD,sample_data(ps.UNC.CD.nonIBD)$Pathology == "NI")
ps.UNC.NI <- prune_taxa(taxa_sums(ps.UNC.NI)>0,ps.UNC.NI)

sample_data(ps.UNC.NI)$subgrp <- paste(sample_data(ps.UNC.NI)$DiseaseStatus,sample_data(ps.UNC.NI)$AnatomSite)

sample_data(ps.UNC.NI)$Recurrence2 <- ifelse(sample_data(ps.UNC.NI)$Rutgeerts_Score %in% c("i0","i1","i0-i1","i1-i2"),"no recurrence",
                                               ifelse(sample_data(ps.UNC.NI)$Rutgeerts_Score %in%c("i2","i2+","i2-i3","i3","i4"),"recurrence"," ")
)

#### Number of Total WashU Patients
#Merge Samples if duplicates
sample_data(ps.WashU.patients)$NewID <- paste(sample_data(ps.WashU.patients)$PatientNo,sample_data(ps.WashU.patients)$DiseaseStatus, sample_data(ps.WashU.patients)$Pathology,sep=" ")
ps.WashU <- merge_samples(ps.WashU.patients,"NewID")

sample_data(ps.WashU)$DiseaseStatus <- as.integer(sample_data(ps.WashU)$DiseaseStatus)
sample_data(ps.WashU)$DiseaseStatus <- as.factor(sample_data(ps.WashU)$DiseaseStatus)
levels(sample_data(ps.WashU)$DiseaseStatus) <- levels(sample_data(ps.WashU.patients)$DiseaseStatus)

sample_data(ps.WashU)$AnatomSite <- as.integer(sample_data(ps.WashU)$AnatomSite)
sample_data(ps.WashU)$AnatomSite <- as.factor(sample_data(ps.WashU)$AnatomSite)
levels(sample_data(ps.WashU)$AnatomSite) <- levels(sample_data(ps.WashU.patients)$AnatomSite)

sample_data(ps.WashU)$AnatomSite <- as.integer(sample_data(ps.WashU)$AnatomSite)
sample_data(ps.WashU)$AnatomSite <- as.factor(sample_data(ps.WashU)$AnatomSite)
levels(sample_data(ps.WashU)$AnatomSite) <- levels(sample_data(ps.WashU.patients)$AnatomSite)

sample_data(ps.WashU)$Race <- as.integer(sample_data(ps.WashU)$Race)
sample_data(ps.WashU)$Race  <- as.factor(sample_data(ps.WashU)$Race)
levels(sample_data(ps.WashU)$Race) <- levels(sample_data(ps.WashU.patients)$Race)

sample_data(ps.WashU)$Ethnicity <- as.integer(sample_data(ps.WashU)$Ethnicity)
sample_data(ps.WashU)$Ethnicity  <- as.factor(sample_data(ps.WashU)$Ethnicity)
levels(sample_data(ps.WashU)$Ethnicity) <- levels(sample_data(ps.WashU.patients)$Ethnicity)

sample_data(ps.WashU)$Sex <- as.integer(sample_data(ps.WashU)$Sex)
sample_data(ps.WashU)$Sex  <- as.factor(sample_data(ps.WashU)$Sex)
levels(sample_data(ps.WashU)$Sex) <- levels(sample_data(ps.WashU.patients)$Sex)

sample_data(ps.WashU)$Age <- as.integer(sample_data(ps.WashU)$Age)

sample_data(ps.WashU)$Behavior <- as.integer(sample_data(ps.WashU)$Behavior)
sample_data(ps.WashU)$Behavior  <- as.factor(sample_data(ps.WashU)$Behavior)
levels(sample_data(ps.WashU)$Behavior) <- levels(sample_data(ps.WashU.patients)$Behavior)

sample_data(ps.WashU)$Location <- as.integer(sample_data(ps.WashU)$Location)
sample_data(ps.WashU)$Location  <- as.factor(sample_data(ps.WashU)$Location)
levels(sample_data(ps.WashU)$Location) <- levels(sample_data(ps.WashU.patients)$Location)

sample_data(ps.WashU)$SmokingStatus <- as.integer(sample_data(ps.WashU)$SmokingStatus)
sample_data(ps.WashU)$SmokingStatus  <- as.factor(sample_data(ps.WashU)$SmokingStatus)
levels(sample_data(ps.WashU)$SmokingStatus) <- levels(sample_data(ps.WashU.patients)$SmokingStatus)

sample_data(ps.WashU)$X5_ASA <- as.integer(sample_data(ps.WashU)$X5_ASA)
sample_data(ps.WashU)$X5_ASA  <- as.factor(sample_data(ps.WashU)$X5_ASA)
levels(sample_data(ps.WashU)$X5_ASA) <- levels(sample_data(ps.WashU.patients)$X5_ASA)

sample_data(ps.WashU)$Immunomodulators <- as.integer(sample_data(ps.WashU)$Immunomodulators)
sample_data(ps.WashU)$Immunomodulators  <- as.factor(sample_data(ps.WashU)$Immunomodulators)
levels(sample_data(ps.WashU)$Immunomodulators) <- levels(sample_data(ps.WashU.patients)$Immunomodulators)

sample_data(ps.WashU)$Steroids <- as.integer(sample_data(ps.WashU)$Steroids)
sample_data(ps.WashU)$Steroids  <- as.factor(sample_data(ps.WashU)$Steroids)
levels(sample_data(ps.WashU)$Steroids) <- levels(sample_data(ps.WashU.patients)$Steroids)

sample_data(ps.WashU)$TNF <- as.integer(sample_data(ps.WashU)$TNF)
sample_data(ps.WashU)$TNF  <- as.factor(sample_data(ps.WashU)$TNF)
levels(sample_data(ps.WashU)$TNF) <- levels(sample_data(ps.WashU.patients)$TNF)

sample_data(ps.WashU)$Rutgeerts_Score <- as.integer(sample_data(ps.WashU)$Rutgeerts_Score)
sample_data(ps.WashU)$Rutgeerts_Score  <- as.factor(sample_data(ps.WashU)$Rutgeerts_Score)
levels(sample_data(ps.WashU)$Rutgeerts_Score) <- levels(sample_data(ps.WashU.patients)$Rutgeerts_Score)

sample_data(ps.WashU)$Recurrence2 <- ifelse(sample_data(ps.WashU)$Rutgeerts_Score %in% c("i0","i1","i0-i1","i1-i2"),"no recurrence",
                                             ifelse(sample_data(ps.WashU)$Rutgeerts_Score %in%c("i2","i2+","i2-i3","i3","i4"),"recurrence"," ")
)
#CD - 79

#Total CD samples: 46+79 = 125
#Total nonIBD samples:  23
#### Number of UNC Samples Used####
Microbiome_Map_20190713 <- read.delim("C:/Users/nmshahir/Dropbox/IBD_microbiota/Analyses/Metadata/Microbiome_Map_20190713.txt")
UNC_Map_Only <- subset(Microbiome_Map_20190713,Microbiome_Map_20190713$Cohort == "UNC")
UNC_Map_NI_I_Only <- subset(UNC_Map_Only,UNC_Map_Only$Pathology %in% c("NI","I"))
UNC_CD_nonIBD_Only <- subset(UNC_Map_NI_I_Only,UNC_Map_NI_I_Only$DiseaseStatus %in% c("CD","nonIBD"))


WashU_Map_Only <- subset(Microbiome_Map_20190713,Microbiome_Map_20190713$Cohort == "WashU")
WashU_Map_NI_Only <- subset(WashU_Map_Only,WashU_Map_Only$Pathology %in% c("NI"))

patient.num.UNC <- sample_data(ps.UNC.patients)$PatientNo
patient.num.WashU <- sample_data(ps.WashU.patients)$PatientNo

UNC_Manuscript_samples_used <- subset(UNC_CD_nonIBD_Only,UNC_CD_nonIBD_Only$PatientNo %in% patient.num.UNC)
WashU_Manuscript_samples_used <- subset(WashU_Map_NI_Only,WashU_Map_NI_Only$PatientNo %in% patient.num.WashU)

#113 Samples from UNC used
#99 Samples from WashU used

####Average Seq Depth ####

ps.UNC.CD.nonIBD <- subset_samples(ps.UNC,sample_data(ps.UNC)$DiseaseStatus %in% c("CD","nonIBD"))
ps.UNC.CD.nonIBD <- prune_taxa(taxa_sums(ps.UNC.CD.nonIBD)>0,ps.UNC.CD.nonIBD)

total_seqs <- sum(sample_sums(ps.UNC.CD.nonIBD)) + sum(sample_sums(ps.WashU))
ave_seq <- total_seqs/212
