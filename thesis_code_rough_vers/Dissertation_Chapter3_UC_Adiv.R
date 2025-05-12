####Manuscript - Alpha Diversity Analysis####
#ps.patients.g <- readRDS("ps_patients_genus_20200303.rds")
ps.patients.g <- readRDS("ps_patients_genus_20190713.rds")
ps.patients.g <- readRDS("ps_patients_species_20190713.rds")

sample_data(ps.patients.g)$Age <- as.numeric(levels(sample_data(ps.patients.g)$Age))[sample_data(ps.patients.g)$Age]

#Get all samples above 8000 reads
ps.patients.g.10K <- prune_samples(sample_sums(ps.patients.g)>8000,ps.patients.g)
ps.patients.g.10K <- prune_taxa(taxa_sums(ps.patients.g.10K)>0,ps.patients.g.10K)

ps.temp <- ps.patients.g.10K

old_meta <- sample_data(ps.temp)
write.csv(old_meta,"old_metadata.csv")
newmap = import_qiime_sample_data("Updated_UC_metadata.txt")

sample_data(ps.temp) <- newmap

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

#### Make Analysis Subsets ####
#### Make Subsets ####

#ileal disease
ileum.patients <- subset_samples(ps.UNC.patients,sample_data(ps.UNC.patients)$AnatomSite == "ileum")
nonIBD.UC.ileum <- subset_samples(ileum.patients, sample_data(ileum.patients)$DiseaseStatus %in% c("nonIBD","UC"))
nonIBD.UC.NI.ileum <- subset_samples(nonIBD.UC.ileum,sample_data(nonIBD.UC.ileum)$Pathology == "NI")
nonIBD.UC.NI.ileum <- prune_taxa(taxa_sums(nonIBD.UC.NI.ileum)>0, nonIBD.UC.NI.ileum)

#colon disease
colon.patients <- subset_samples(ps.UNC.patients,sample_data(ps.UNC.patients)$AnatomSite == "colon")

nonIBD.IBD.colon <- subset_samples(colon.patients,sample_data(colon.patients)$DiseaseStatus %in% c("nonIBD","CD","UC"))
nonIBD.IBD.colon <- prune_taxa(taxa_sums(nonIBD.IBD.colon) >0, nonIBD.IBD.colon)
nonIBD.IBD.NI.colon <- subset_samples(nonIBD.IBD.colon,sample_data(nonIBD.IBD.colon)$Pathology == "NI")
nonIBD.IBD.NI.colon <- prune_taxa(taxa_sums(nonIBD.IBD.NI.colon)>0,nonIBD.IBD.NI.colon)

nonIBD.UC.colon <- subset_samples(colon.patients, sample_data(colon.patients)$DiseaseStatus %in% c("nonIBD","UC"))
nonIBD.UC.NI.colon <- subset_samples(nonIBD.UC.colon,sample_data(nonIBD.UC.colon)$Pathology == "NI")
nonIBD.UC.NI.colon <- prune_taxa(taxa_sums(nonIBD.UC.NI.colon)>0, nonIBD.UC.NI.colon)

# UC Subgroups
UC.patients <- subset_samples(ps.UNC.patients,sample_data(ps.UNC.patients)$DiseaseStatus == "UC")
UC.patients <- prune_taxa(taxa_sums(UC.patients)>0,UC.patients)

UC.NI.patients <- subset_samples(UC.patients,sample_data(UC.patients)$Pathology == "NI")
UC.NI.patients <- prune_taxa(taxa_sums(UC.NI.patients)>0,UC.NI.patients)

UC.colon <- subset_samples(UC.patients,sample_data(UC.patients)$AnatomSite == "colon")
UC.colon <- prune_taxa(taxa_sums(UC.colon)>0,UC.colon)

UC.NI.colon <- subset_samples(nonIBD.UC.NI.colon, sample_data(nonIBD.UC.NI.colon)$DiseaseStatus == "UC")
UC.NI.colon <- prune_taxa(taxa_sums(UC.NI.colon) > 0, UC.NI.colon)

UC.NI.ileum <- subset_samples(nonIBD.UC.NI.ileum, sample_data(nonIBD.UC.NI.ileum)$DiseaseStatus == "UC")
UC.NI.ileum <- prune_taxa(taxa_sums(UC.NI.ileum) > 0, UC.NI.ileum)

#matched inflamed and non-inflamed samples
UC.colon.m <- subset_samples(UC.colon, sample_data(UC.colon)$PatientNo %in% c("40","41200","42300","42700","43000","46","52"))
UC.colon.m <- prune_taxa(taxa_sums(UC.colon.m)>0,UC.colon.m)

# UC and CD colon
CD.UC.colon <- subset_samples(colon.patients,sample_data(colon.patients)$DiseaseStatus %in% c("CD","UC"))
CD.UC.NI.colon <- subset_samples(CD.UC.colon,sample_data(CD.UC.colon)$Pathology == "NI")
CD.UC.NI.colon <- prune_taxa(taxa_sums(CD.UC.NI.colon)>0, CD.UC.NI.colon)

# CD and nonIBD
nonIBD.CD.colon <- subset_samples(colon.patients,sample_data(colon.patients)$DiseaseStatus %in% c("CD","nonIBD"))
nonIBD.CD.NI.colon <- subset_samples(nonIBD.CD.colon,sample_data(nonIBD.CD.colon)$Pathology == "NI")
nonIBD.CD.NI.colon <- prune_taxa(taxa_sums(nonIBD.CD.NI.colon)>0,nonIBD.CD.NI.colon)

nonIBD.NI.colon <- subset_samples(nonIBD.CD.NI.colon, sample_data(nonIBD.CD.NI.colon)$DiseaseStatus == "nonIBD")
nonIBD.NI.colon <- prune_taxa(taxa_sums(nonIBD.NI.colon)>0,nonIBD.NI.colon)

CD.NI.colon <- subset_samples(nonIBD.CD.NI.colon, sample_data(nonIBD.CD.NI.colon)$DiseaseStatus == "CD")
CD.NI.colon <- prune_taxa(taxa_sums(CD.NI.colon)>0,CD.NI.colon)


UC.only <- getThresholdTaxa(UC.NI.colon,0)
CD.only <- getThresholdTaxa(CD.NI.colon,0)
nonIBD.only <- getThresholdTaxa(nonIBD.NI.colon,0)


##### PHYLUM OVERVIEW ####
phyla_counts_tab <- otu_table(tax_glom(nonIBD.IBD.NI.colon, taxrank="Phylum")) 
phyla_tax_vec <- as.vector(tax_table(tax_glom(nonIBD.IBD.NI.colon, taxrank="Phylum"))[,2]) 
rownames(phyla_counts_tab) <- phyla_tax_vec
#Alpha Diversity Tests
####nonIBD and CD colon ####
data <- nonIBD.UC.NI.colon
divnet.nonIBD.UC.colon <- divnet(data,
                                 X= "DiseaseStatus",
                                 ncores = 4)
saveRDS(divnet.nonIBD.UC.colon,"colon_nonIBD_v_UC.rds")

####nonIBD and CD ileum ####
data <- nonIBD.UC.NI.ileum
divnet.nonIBD.UC.ileum <- divnet(data,
                                 X= "DiseaseStatus",
                                 ncores = 4)
saveRDS(divnet.nonIBD.UC.ileum,"ileum_nonIBD_v_UC.rds")

####UC colon - inflammation ####
data <- UC.colon.m
divnet.UC.colon.inf <- divnet(data,
                                 X= "Pathology",
                                 ncores = 4)
saveRDS(divnet.UC.colon.inf,"UC_colon_inf.rds")

#Get UC Classification Data from Ashley/Caroline
data <- UC.NI.colon
divnet.UC.colon.E23 <- divnet(data,
                              X= "ClassificationUC",
                              ncores = 4)
saveRDS(divnet.UC.colon.E23,"UC_NI_colon_E23.rds")

#### UC colon - 5-ASA ####
data <- UC.NI.colon
divnet.UC.colon.asa <- divnet(data,
                              X= "X5_ASA",
                              ncores = 4)
saveRDS(divnet.UC.colon.asa,"UC_NI_colon_5_ASA.rds")

#### UC colon - Immunomodulators ####
data <- UC.NI.colon
divnet.UC.colon.imm <- divnet(data,
                              X= "Immunomodulators",
                              ncores = 4)
saveRDS(divnet.UC.colon.imm,"UC_NI_colon_IMM.rds")

#### UC colon - Probiotics ####
data <- UC.NI.colon
divnet.UC.colon.pro <- divnet(data,
                              X= "Probiotic",
                              ncores = 4)
saveRDS(divnet.UC.colon.pro,"UC_NI_colon_PRO.rds")

#### UC colon - Steroids ####
data <- UC.NI.colon
divnet.UC.colon.str <- divnet(data,
                              X= "Steroids",
                              ncores = 4)
saveRDS(divnet.UC.colon.str,"UC_NI_colon_STR.rds")

#### UC colon - TNF ####
data <- UC.NI.colon
divnet.UC.colon.tnf <- divnet(data,
                              X= "TNF",
                              ncores = 4)
saveRDS(divnet.UC.colon.tnf,"UC_NI_colon_TNF.rds")

#### UC Colon - Inflammation All ####
data <- UC.colon
data.inf <- subset_samples(UC.colon,sample_data(UC.colon)$Pathology %in% c("NI","I"))
data.inf <- prune_taxa(taxa_sums(data.inf)>0,data.inf)
divnet.UC.colon.inf.all <- divnet(data.inf,
                                  X = "Pathology",
                                  ncores = 4)
saveRDS(divnet.UC.colon.inf.all,"UC_colon_all_inf.rds")

#### Colon - UC vs CD ####
data <- CD.UC.NI.colon
divnet.CD.UC.colon <- divnet(data,
                             X = "DiseaseStatus",
                             ncores = 4)
saveRDS(divnet.CD.UC.colon,"colon_CD_UC.rds")

#### Colon - UC E2 vs nonIBD ####
data <- nonIBD.UC.NI.colon
sample_data(data)$DiseasePheno <- paste(sample_data(data)$DiseaseStatus,sample_data(data)$ClassificationUC,sep = " ")
data.E2 <- subset_samples(data,sample_data(data)$DiseasePheno %in% c("UC E2", "nonIBD "))
data.E2 <- prune_taxa(taxa_sums(data.E2)>0,data.E2)

divnet.nonIBD.UC.E2 <- divnet(data.E2,
                              X="DiseaseStatus",
                              ncores = 4)
saveRDS(divnet.nonIBD.UC.E2,"colon_nonIBD_v_UC_E2.rds")

#### Colon - UC E3 vs nonIBD ####
data <- nonIBD.UC.NI.colon
sample_data(data)$DiseasePheno <- paste(sample_data(data)$DiseaseStatus,sample_data(data)$ClassificationUC,sep = " ")
data.E3 <- subset_samples(data,sample_data(data)$DiseasePheno %in% c("UC E3", "nonIBD "))
data.E3 <- prune_taxa(taxa_sums(data.E3)>0,data.E3)

divnet.nonIBD.UC.E3 <- divnet(data.E3,
                              X="DiseaseStatus",
                              ncores = 4)
saveRDS(divnet.nonIBD.UC.E3,"colon_nonIBD_v_UC_E3.rds")

#### Colon - UC Subtypes ####
sample_data(UC.NI.colon)$AeroLvl <- ifelse(sample_data(UC.NI.colon)$PatientNo %in% c("127","41500","42700","43000","46","52"),"Low","Health")
#REUPLOAD
data <- UC.NI.colon
divnet.UC.NI.colon.aero <- divnet(data,
                              X="AeroLvl",
                              ncores = 4)
saveRDS(divnet.UC.NI.colon.aero,"UC_NI_colon_aero_split.rds")

#### Colon - nonIBD vs UC Low ####
sample_data(nonIBD.UC.NI.colon)$AeroLvl <- ifelse(sample_data(nonIBD.UC.NI.colon)$DiseaseStatus=="nonIBD","",ifelse(sample_data(nonIBD.UC.NI.colon)$PatientNo %in% c("127","41500","42700","43000","46","52"),"Low","Health"))
sample_data(nonIBD.UC.NI.colon)$AeroDisease <- paste(sample_data(nonIBD.UC.NI.colon)$DiseaseStatus,sample_data(nonIBD.UC.NI.colon)$AeroLvl, sep = " ")

data <- nonIBD.UC.NI.colon
data.nonIBD.Low <- subset_samples(data,sample_data(data)$AeroDisease %in% c("nonIBD ","UC Low"))
data.nonIBD.Low <- prune_taxa(taxa_sums(data.nonIBD.Low)>0,data.nonIBD.Low)
data <- data.nonIBD.Low
divnet.nonIBD.UC.low <- divnet(data,
                           X="DiseaseStatus",
                           ncores = 4)
saveRDS(divnet.nonIBD.UC.low,"nonIBD_v_UC_low_colon.rds")

#### Colon - nonIBD vs UC High ####
data <- nonIBD.UC.NI.colon
data.nonIBD.H <- subset_samples(data,sample_data(data)$AeroDisease %in% c("nonIBD ","UC Health"))
data.nonIBD.H <- prune_taxa(taxa_sums(data.nonIBD.H)>0,data.nonIBD.H)
data <- data.nonIBD.H
divnet.nonIBD.UC.health <- divnet(data,
                               X="DiseaseStatus",
                               ncores = 4)
saveRDS(divnet.nonIBD.UC.health,"nonIBD_v_UC_health_colon.rds")

#### Colon - CD vs UC Low ####
sample_data(CD.UC.NI.colon)$AeroLvl <- ifelse(sample_data(CD.UC.NI.colon)$DiseaseStatus=="CD","",ifelse(sample_data(CD.UC.NI.colon)$PatientNo %in% c("127","41500","42700","43000","46","52"),"Low","Health"))
sample_data(CD.UC.NI.colon)$AeroDisease <- paste(sample_data(CD.UC.NI.colon)$DiseaseStatus,sample_data(CD.UC.NI.colon)$AeroLvl, sep = " ")

data <- CD.UC.NI.colon
data.CD.Low <- subset_samples(data,sample_data(data)$AeroDisease %in% c("CD ","UC Low"))
data.CD.Low <- prune_taxa(taxa_sums(data.CD.Low)>0,data.CD.Low)
data <- data.CD.Low
divnet.CD.UC.low <- divnet(data,
                           X="DiseaseStatus",
                           ncores = 4)
saveRDS(divnet.CD.UC.low,"CD_v_UC_low_colon.rds")

#### Colon - CD vs UC High ####
data <- CD.UC.NI.colon
data.CD.H <- subset_samples(data,sample_data(data)$AeroDisease %in% c("CD ","UC Health"))
data.CD.H <- prune_taxa(taxa_sums(data.CD.H)>0,data.CD.H)
data <- data.CD.H
divnet.CD.UC.h <- divnet(data,
                         X="DiseaseStatus",
                         ncores = 4)
saveRDS(divnet.CD.UC.h,"CD_v_UC_health_colon.rds")