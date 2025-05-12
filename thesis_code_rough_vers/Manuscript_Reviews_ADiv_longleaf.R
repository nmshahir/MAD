library(phyloseq)
library(magrittr)
library(ggplot2)
####Manuscript - Alpha Diversity Analysis####
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


#### Make Analysis Subsets ####
#### Make Subsets ####

#ileal disease
ileum.patients <- subset_samples(ps.UNC.patients,sample_data(ps.UNC.patients)$AnatomSite == "ileum")
nonIBD.CD.ileum <- subset_samples(ileum.patients, sample_data(ileum.patients)$DiseaseStatus %in% c("nonIBD","CD"))
nonIBD.CD.NI.ileum <- subset_samples(nonIBD.CD.ileum,sample_data(nonIBD.CD.ileum)$Pathology == "NI")
nonIBD.CD.NI.ileum <- prune_taxa(taxa_sums(nonIBD.CD.NI.ileum)>0, nonIBD.CD.NI.ileum)

#colon disease
colon.patients <- subset_samples(ps.UNC.patients,sample_data(ps.UNC.patients)$AnatomSite == "colon")
nonIBD.CD.colon <- subset_samples(colon.patients, sample_data(colon.patients)$DiseaseStatus %in% c("nonIBD","CD"))
nonIBD.CD.NI.colon <- subset_samples(nonIBD.CD.colon,sample_data(nonIBD.CD.colon)$Pathology == "NI")
nonIBD.CD.NI.colon <- prune_taxa(taxa_sums(nonIBD.CD.NI.colon)>0, nonIBD.CD.NI.colon)

# CD Subgroups
CD.patients <- subset_samples(ps.UNC.patients,sample_data(ps.UNC.patients)$DiseaseStatus == "CD")
CD.patients <- prune_taxa(taxa_sums(CD.patients)>0,CD.patients)

CD.NI.patients <- subset_samples(CD.patients,sample_data(CD.patients)$Pathology == "NI")
CD.NI.patients <- prune_taxa(taxa_sums(CD.NI.patients)>0,CD.NI.patients)

CD.ileum <- subset_samples(CD.patients,sample_data(CD.patients)$AnatomSite == "ileum")
CD.ileum <- prune_taxa(taxa_sums(CD.ileum) > 0, CD.ileum)

#CD.ileum.m <- subset_samples(CD.ileum,sample_data(CD.ileum)$PatientNo %in% c("1004098","129","137","138","144","15","19","40100","41","41700","42100","42400","42600","44","44000","47","51","56","57","59","60","70","77"))
#CD.ileum.m <- prune_taxa(taxa_sums(CD.ileum.m) > 0, CD.ileum.m)

CD.ileum.m <- subset_samples(CD.ileum,sample_data(CD.ileum)$PatientNo %in% c("129","137","138","144","15","40100","41","41700","42100","42200","42400","42600","44000","51","57","77"))
CD.ileum.m <- prune_taxa(taxa_sums(CD.ileum.m) > 0, CD.ileum.m)

CD.colon <- subset_samples(CD.patients,sample_data(CD.patients)$AnatomSite == "colon")
CD.colon <- prune_taxa(taxa_sums(CD.colon) > 0, CD.colon)

CD.colon.m <- subset_samples(CD.colon,sample_data(CD.colon)$PatientNo %in% c("40200","40300","40500","40700","40800","42000","42900","63"))
CD.colon.m <- prune_taxa(taxa_sums(CD.colon.m) >0, CD.colon.m)

CD.NI.colon <- subset_samples(nonIBD.CD.NI.colon, sample_data(nonIBD.CD.NI.colon)$DiseaseStatus == "CD")
CD.NI.colon <- prune_taxa(taxa_sums(CD.NI.colon) > 0, CD.NI.colon)

CD.NI.ileum <- subset_samples(nonIBD.CD.NI.ileum, sample_data(nonIBD.CD.NI.ileum)$DiseaseStatus == "CD")
CD.NI.ileum <- prune_taxa(taxa_sums(CD.NI.ileum) > 0, CD.NI.ileum)

# nonIBD Subgroups
nonIBD.patients <- subset_samples(ps.UNC.patients,sample_data(ps.UNC.patients)$DiseaseStatus == "nonIBD")
nonIBD.patients <- prune_taxa(taxa_sums(nonIBD.patients)>0,nonIBD.patients)

nonIBD.NI.patients <- subset_samples(nonIBD.patients,sample_data(nonIBD.patients)$Pathology == "NI")
nonIBD.NI.patients <- prune_taxa(taxa_sums(nonIBD.NI.patients)>0,nonIBD.NI.patients)

nonIBD.NI.colon <- subset_samples(nonIBD.CD.NI.colon,sample_data(nonIBD.CD.NI.colon)$DiseaseStatus == "nonIBD")
nonIBD.NI.colon <- prune_taxa(taxa_sums(nonIBD.NI.colon)>0,nonIBD.NI.colon)

nonIBD.NI.ileum <- subset_samples(nonIBD.CD.NI.ileum,sample_data(nonIBD.CD.NI.ileum)$DiseaseStatus == "nonIBD")
nonIBD.NI.ileum <- prune_taxa(taxa_sums(nonIBD.NI.ileum)>0,nonIBD.NI.ileum)

####nonIBD and CD colon ####
data <- nonIBD.CD.NI.colon
divnet.nonIBD.CD.colon <- divnet(data,
                                 X= "DiseaseStatus",
                                 ncores = 24)
saveRDS(divnet.nonIBD.CD.colon,"colon_nonIBDvCD_all.rds")

####nonIBD and CD ileum ####
data <- nonIBD.CD.NI.ileum
divnet.nonIBD.CD.ileum <- divnet(data,
                                 X= "DiseaseStatus",
                                 ncores = 24)
saveRDS(divnet.nonIBD.CD.ileum,"ileum_disease.rds")

####nonIBD tissue ####
data <- nonIBD.NI.patients
divnet.nonIBD.tissue <- divnet(data,
                               X= "AnatomSite",
                               ncores = 24)
saveRDS(divnet.nonIBD.tissue,"nonIBD_tissue.rds")

####CD tissue ####
data <- CD.NI.patients
divnet.CD.tissue <- divnet(data,
                           X= "AnatomSite",
                           ncores = 24)
saveRDS(divnet.CD.tissue,"CD_tissue.rds")

#Colon CD: B1 vs B2
data <- CD.NI.colon
data <- subset_samples(data,sample_data(data)$Behavior %in% c("B1","B2"))
data <- prune_taxa(taxa_sums(data) > 0, data)
start.time.b12 <- Sys.time()
divnet.CD.B12.colon <- divnet(data,
                              X= "Behavior",
                              ncores = 24)
end.time.b12 <- Sys.time()

#Colon CD: B1 vs B3
data <- CD.NI.colon
data <- subset_samples(data,sample_data(data)$Behavior %in% c("B1","B3"))
data <- prune_taxa(taxa_sums(data) > 0, data)
start.time.b13 <- Sys.time()
divnet.CD.B13.colon <- divnet(data,
                              X= "Behavior",
                              ncores = 24)
end.time.b13 <- Sys.time()

#Colon CD: B2 vs B3
data <- CD.NI.colon
data <- subset_samples(data,sample_data(data)$Behavior %in% c("B2","B3"))
data <- prune_taxa(taxa_sums(data) > 0, data)
start.time.b23 <- Sys.time()
divnet.CD.B23.colon <- divnet(data,
                              X= "Behavior",
                              ncores = 24)
end.time.b23 <- Sys.time()

#Colon CD: L2 vs L3
data <- CD.NI.colon
data <- subset_samples(data,sample_data(data)$Location %in% c("L2","L3"))
data <- prune_taxa(taxa_sums(data) > 0, data)
start.time.l23 <- Sys.time()
divnet.CD.L23.colon <- divnet(data,
                              X= "Location",
                              ncores = 24)
end.time.l23 <- Sys.time()

saveRDS(divnet.CD.B12.colon,"CD_colon_B1B2.rds")
saveRDS(divnet.CD.B13.colon,"CD_colon_B1B3.rds")
saveRDS(divnet.CD.B23.colon,"CD_colon_B2B3.rds")
saveRDS(divnet.CD.L23.colon,"CD_colon_L2L3.rds")

#Ileum: B2 vs B3
data <- CD.NI.ileum
data <- subset_samples(data,sample_data(data)$Behavior %in% c("B2","B3"))
data <- prune_taxa(taxa_sums(data) > 0, data)
start.time.ileum.b23 <- Sys.time()
divnet.CD.B23.ileum <- divnet(data,
                              X= "Behavior",
                              ncores = 24)
end.time.ileum.b23 <- Sys.time()


#CD B3 only
data <- CD.NI.patients
data <- subset_samples(data,sample_data(data)$Behavior %in% c("B3"))
data <- prune_taxa(taxa_sums(data) > 0, data)
start.time.B3 <- Sys.time()
divnet.CD.B3 <- divnet(data,
                       X= "AnatomSite",
                       ncores = 24)
end.time.B3 <- Sys.time()

saveRDS(divnet.CD.B23.ileum,"CD_ileum_B2B3.rds")
saveRDS(divnet.CD.B3,"CD_B3_tissue.rds")

data <- CD.NI.patients
data <- subset_samples(data,sample_data(data)$Location %in% c("L3"))
data <- prune_taxa(taxa_sums(data) > 0, data)
start.time.L3 <- Sys.time()
divnet.CD.L3 <- divnet(data,
                       X= "AnatomSite",
                       ncores = 24)
end.time.L3 <- Sys.time()

saveRDS(divnet.CD.L3,"CD_L3_tissue.rds")


##### remove wrong L3 ###
data <- CD.NI.patients
data <- subset_samples(data,sample_data(data)$Location %in% c("L3"))
data <- prune_taxa(taxa_sums(data) > 0, data)

data <- subset_samples(data,!(sample_data(data)$Behavior %in% c("B1")))
data <- prune_taxa(taxa_sums(data) > 0, data)

data.concordant <- subset_samples(data,!(sample_data(data)$PatientNo %in% c("92","96")))
data.concordant <- prune_taxa(taxa_sums(data.concordant) > 0, data.concordant)
data <- data.concordant

data.concordant <- subset_samples(data,sample_names(data) != "51 CD NI ileum UNC")
data.concordant <- prune_taxa(taxa_sums(data.concordant) > 0, data.concordant)

data <- data.concordant
start.time.L3 <- Sys.time()
divnet.CD.L3.no.disconcordant <- divnet(data,
                                        X= "AnatomSite",
                                        ncores = 24)
end.time.L3 <- Sys.time()

saveRDS(divnet.CD.L3.no.disconcordant,"CD_L3_tissue_remove_wrong_stricture.rds")

#CD L3, remove B1s
data <- CD.NI.patients
data <- subset_samples(data,sample_data(data)$Location %in% c("L3"))
data <- prune_taxa(taxa_sums(data) > 0, data)

data <- subset_samples(data,!(sample_data(data)$Behavior %in% c("B1")))
data <- prune_taxa(taxa_sums(data) >0, data)
start.time.L3 <- Sys.time()
divnet.CD.L3.n.b1 <- divnet(data,
                            X= "AnatomSite",
                            ncores = 24)
end.time.L3 <- Sys.time()

saveRDS(divnet.CD.L3.n.b1,"CD_L3_tissue_no_B1.rds")

#CD B2,B3 - Tissue
data <- CD.NI.patients
data <- subset_samples(data,sample_data(data)$Behavior %in% c("B2","B3"))
data <- prune_taxa(taxa_sums(data) > 0, data)
start.time.B23 <- Sys.time()
divnet.CD.B23 <- divnet(data,
                        X= "AnatomSite",
                        ncores = 24)
end.time.B23 <- Sys.time()
saveRDS(divnet.CD.B23,"CD_B2B3_tissue.rds")

# #CD colon: B1,B2 - Tissue (remove strictures)
# data <- CD.NI.colon
# data.concordant <- subset_samples(data,!(sample_data(data)$PatientNo %in% c("92","96")))
# data.concordant <- prune_taxa(taxa_sums(data.concordant) > 0, data.concordant)
# data <- data.concordant
# 
# data <- subset_samples(data,sample_data(data)$Behavior %in% c("B1","B2"))
# data <- prune_taxa(taxa_sums(data) > 0, data)
# start.time.B12 <- Sys.time()
# divnet.CD.colon.B12.concordant <- divnet(data,
#                                          X= "Behavior",
#                                          ncores = 24)
# end.time.B23 <- Sys.time()
# saveRDS(divnet.CD.colon.B12.concordant,"CD_colon_B1B2_match_stricturing_tissue.rds")
# 
# #CD colon: B3,B2 - Tissue (remove strictures)
# data <- CD.NI.colon
# data.concordant <- subset_samples(data,!(sample_data(data)$PatientNo %in% c("92","96")))
# data.concordant <- prune_taxa(taxa_sums(data.concordant) > 0, data.concordant)
# data <- data.concordant
# 
# data <- subset_samples(data,sample_data(data)$Behavior %in% c("B3","B2"))
# data <- prune_taxa(taxa_sums(data) > 0, data)
# start.time.B12 <- Sys.time()
# divnet.CD.colon.B23.concordant <- divnet(data,
#                                          X= "Behavior",
#                                          ncores = 24)
# end.time.B23 <- Sys.time()
# saveRDS(divnet.CD.colon.B23.concordant,"CD_colon_B2B3_match_stricturing_tissue.rds")

# #CD B2,B3 - Tissue (remove strictures)
# data <- CD.NI.patients
# data.concordant <- subset_samples(data,!(sample_data(data)$PatientNo %in% c("92","96")))
# data.concordant <- prune_taxa(taxa_sums(data.concordant) > 0, data.concordant)
# data <- data.concordant
# 
# data.concordant <- subset_samples(data,sample_names(data) != "51 CD NI ileum UNC")
# data.concordant <- prune_taxa(taxa_sums(data.concordant) > 0, data.concordant)
# 
# data <- data.concordant
# 
# data <- subset_samples(data,sample_data(data)$Behavior %in% c("B2","B3"))
# data <- prune_taxa(taxa_sums(data) > 0, data)
# start.time.B23 <- Sys.time()
# divnet.CD.B23.concordant <- divnet(data,
#                                    X= "AnatomSite",
#                                    ncores = 24)
# end.time.B23 <- Sys.time()
# saveRDS(divnet.CD.B23.concordant,"CD_B2B3_match_stricturing_tissue.rds")

# CD Stricturing 
data <- CD.NI.patients
data.B2 <- subset_samples(data,sample_data(data)$Behavior == "B2")
data.B2 <- prune_taxa(taxa_sums(data.B2) > 0, data.B2)

sample_data(data.B2)$Stricture <- ifelse(sample_data(data.B2)$PatientNo %in% c("96","92","57","76","93","138","149","41400","42600","43600","41","45","137","145","42100","42200","44000"),"ileal","colon")
table(sample_data(data.B2)$Stricture)

# #Comparing ileal strictures across tissues (2 vs 15....)
# data.B2.ileal <- subset_samples(data.B2,sample_data(data.B2)$Stricture == "ileal")
# data.B2.ileal <- prune_taxa(taxa_sums(data.B2.ileal)>0,data.B2.ileal)
# 
# data <- data.B2.ileal
# 
# start.time.CD.B2.ileal <- Sys.time()
# divnet.CD.B2.ileal.strict <- divnet(data,
#                                     X= "AnatomSite",
#                                     ncores = 24)
# end.time.CD.B2.ileal <- Sys.time()
# saveRDS(divnet.CD.B2.ileal.strict,"CD_B2_ileal_stricture_tissue.rds")
# 
# #Comparing strictures within the colon (2 ileal vs 5 colonic)
# data.B2.colon.tissue <- subset_samples(data.B2,sample_data(data.B2)$AnatomSite == "colon")
# data.B2.colon.tissue <- prune_taxa(taxa_sums(data.B2.colon.tissue)>0,data.B2.colon.tissue)
# data <- data.B2.colon.tissue
# 
# start.time.CD.B2.colon <- Sys.time()
# divnet.CD.B2.colon <- divnet(data,
#                              X= "Stricture",
#                              ncores = 24)
# end.time.CD.B2.colon <- Sys.time()
# 
# saveRDS(divnet.CD.B2.colon,"CD_B2_colon_strictures.rds")

#Comparing colonic strictures in colon vs ileal strictures in ileum
sample_data(data.B2)$TypeStricture <- paste(sample_data(data.B2)$Stricture,sample_data(data.B2)$AnatomSite,sep= " ")

data.B2.type.tissue <- subset_samples(data.B2,sample_data(data.B2)$TypeStricture %in% c("colon colon","ileal ileum"))
data.B2.type.tissue <- prune_taxa(taxa_sums(data.B2.type.tissue)>0,data.B2.type.tissue)

data <- data.B2.type.tissue

start.time.CD.B2.type.tissue <- Sys.time()
divnet.CD.B2.type.tissue <- divnet(data,
                                   X= "AnatomSite",
                                   ncores = 24)
end.time.CD.B2.type.tissue <- Sys.time()

saveRDS(divnet.CD.B2.type.tissue,"CD_B2_colonic_colon_vs_ileal_ileum_strictures.rds")

####################################### nonIBD and CD behaviors #############################
sample_data(nonIBD.CD.NI.colon)$Cat <- paste(sample_data(nonIBD.CD.NI.colon)$DiseaseStatus,sample_data(nonIBD.CD.NI.colon)$Behavior)
sample_data(nonIBD.CD.NI.ileum)$Cat <- paste(sample_data(nonIBD.CD.NI.ileum)$DiseaseStatus,sample_data(nonIBD.CD.NI.ileum)$Behavior)

#nonIBD vs CD B1 - colon
data <- nonIBD.CD.NI.colon
data <- subset_samples(data,sample_data(data)$Cat %in% c("nonIBD ", "CD B1"))
data <- prune_taxa(taxa_sums(data) > 0, data)

divnet.nonIBD.CD.B1.colon <- divnet(data, 
                                    X = "DiseaseStatus",
                                    ncores = 24)
#nonIBD vs CD B2 - colon
data <- nonIBD.CD.NI.colon
data <- subset_samples(data,sample_data(data)$Cat %in% c("nonIBD ", "CD B2"))
data <- prune_taxa(taxa_sums(data) > 0, data)

divnet.nonIBD.CD.B2.colon <- divnet(data, 
                                    X = "DiseaseStatus",
                                    ncores = 24)

# #nonIBD vs CD B2 (no ileal strictures) - colon
# data <- nonIBD.CD.NI.colon
# 
# data.concordant <- subset_samples(data,!(sample_data(data)$PatientNo %in% c("92","96")))
# data.concordant <- prune_taxa(taxa_sums(data.concordant) > 0, data.concordant)
# data <- data.concordant
# 
# data <- subset_samples(data,sample_data(data)$Cat %in% c("nonIBD ", "CD B2"))
# data <- prune_taxa(taxa_sums(data) > 0, data)
# 
# divnet.nonIBD.CD.B2.concordant.colon <- divnet(data, 
#                                                X = "DiseaseStatus",
#                                                ncores = 24)
# saveRDS(divnet.nonIBD.CD.B2.concordant.colon,"colon_nonIBD_CDB2_match_stricture.rds")

#nonIBD vs CD B3 - colon
data <- nonIBD.CD.NI.colon
data <- subset_samples(data,sample_data(data)$Cat %in% c("nonIBD ", "CD B3"))
data <- prune_taxa(taxa_sums(data) > 0, data)

divnet.nonIBD.CD.B3.colon <- divnet(data, 
                                    X = "DiseaseStatus",
                                    ncores = 24)

saveRDS(divnet.nonIBD.CD.B1.colon,"colon_nonIBD_CD_B1.rds")
saveRDS(divnet.nonIBD.CD.B2.colon,"colon_nonIBD_CD_B2.rds")
saveRDS(divnet.nonIBD.CD.B3.colon,"colon_nonIBD_CD_B3.rds")

#nonIBD vs CD B2 - ileum
data <- nonIBD.CD.NI.ileum
data <- subset_samples(data,sample_data(data)$Cat %in% c("nonIBD ", "CD B2"))
data <- prune_taxa(taxa_sums(data) > 0, data)

divnet.nonIBD.CD.B2.ileum <- divnet(data, 
                                    X = "DiseaseStatus",
                                    ncores = 24)

#nonIBD vs CD B3 - ileum
data <- nonIBD.CD.NI.ileum
data <- subset_samples(data,sample_data(data)$Cat %in% c("nonIBD ", "CD B3"))
data <- prune_taxa(taxa_sums(data) > 0, data)

divnet.nonIBD.CD.B3.ileum <- divnet(data, 
                                    X = "DiseaseStatus",
                                    ncores = 24)

saveRDS(divnet.nonIBD.CD.B2.ileum,"ileum_nonIBD_CD_B2.rds")
saveRDS(divnet.nonIBD.CD.B3.ileum,"ileum_nonIBD_CD_B3.rds")

####################################### nonIBD and CD location #############################
sample_data(nonIBD.CD.NI.colon)$Cat <- paste(sample_data(nonIBD.CD.NI.colon)$DiseaseStatus,sample_data(nonIBD.CD.NI.colon)$Location)

#nonIBD vs CD L2 - colon
data <- nonIBD.CD.NI.colon
data <- subset_samples(data,sample_data(data)$Cat %in% c("nonIBD ", "CD L2"))
data <- prune_taxa(taxa_sums(data) > 0, data)

divnet.nonIBD.CD.L2.colon <- divnet(data, 
                                    X = "DiseaseStatus",
                                    ncores = 24)

#nonIBD vs CD L3 - colon
data <- nonIBD.CD.NI.colon
data <- subset_samples(data,sample_data(data)$Cat %in% c("nonIBD ", "CD L3"))
data <- prune_taxa(taxa_sums(data) > 0, data)

divnet.nonIBD.CD.L3.colon <- divnet(data, 
                                    X = "DiseaseStatus",
                                    ncores = 24)

saveRDS(divnet.nonIBD.CD.L2.colon,"colon_nonIBD_CD_L2.rds")
saveRDS(divnet.nonIBD.CD.L3.colon,"colon_nonIBD_CD_L3.rds")

# #nonIBD vs CD L3 - colon stricture only
# data <- nonIBD.CD.NI.colon
# data.concordant <- subset_samples(data,!(sample_data(data)$PatientNo %in% c("92","96")))
# data.concordant <- prune_taxa(taxa_sums(data.concordant) > 0, data.concordant)
# data <- data.concordant
# 
# 
# data <- subset_samples(data,sample_data(data)$Cat %in% c("nonIBD ", "CD L3"))
# data <- prune_taxa(taxa_sums(data) > 0, data)
# 
# divnet.nonIBD.CD.L3.colon.concordant <- divnet(data, 
#                                                X = "DiseaseStatus",
#                                                ncores = 24)
# saveRDS(divnet.nonIBD.CD.L3.colon.concordant,"colon_nonIBD_CD_L3_concordant.rds")

#### Colon - Medication - Immunomodulators ####
data <- CD.NI.colon
data <- subset_samples(data,sample_data(data)$Immunomodulators %in% c("0","1"))
data <- prune_taxa(taxa_sums(data) > 0, data)

start.time.inf <- Sys.time()
divnet.CD.imm.colon <- divnet(data,
                              X= "Immunomodulators",
                              ncores = 24)
end.time.inf <- Sys.time()

#### Colon - Medication - Steroids ####
data <- CD.NI.colon
data <- subset_samples(data,sample_data(data)$Steroids %in% c("0","1"))
data <- prune_taxa(taxa_sums(data) > 0, data)

start.time.inf <- Sys.time()
divnet.CD.str.colon <- divnet(data,
                              X= "Steroids",
                              ncores = 24)
end.time.inf <- Sys.time()

#### Colon - Medication - TNF ####
data <- CD.NI.colon
data <- subset_samples(data,sample_data(data)$TNF %in% c("0","1"))
data <- prune_taxa(taxa_sums(data) > 0, data)

start.time.inf <- Sys.time()
divnet.CD.tnf.colon <- divnet(data,
                              X= "TNF",
                              ncores = 24)
end.time.inf <- Sys.time()

saveRDS(divnet.CD.imm.colon,"CD_colon_IMM_usage.rds")
saveRDS(divnet.CD.str.colon,"CD_colon_STR_usage.rds")
saveRDS(divnet.CD.tnf.colon,"CD_colon_TNF_usage.rds")


#### Ileum - Medication - Immunomodulators ####
data <- CD.NI.ileum
data <- subset_samples(data,sample_data(data)$Immunomodulators %in% c("0","1"))
data <- prune_taxa(taxa_sums(data) > 0, data)

start.time.inf <- Sys.time()
divnet.CD.imm.ileum <- divnet(data,
                              X= "Immunomodulators",
                              ncores = 24)
end.time.inf <- Sys.time()

#### ileum - Medication - 5-ASA ####
data <- CD.NI.ileum
data <- subset_samples(data,sample_data(data)$X5_ASA %in% c("0","1"))
data <- prune_taxa(taxa_sums(data) > 0, data)

start.time.inf <- Sys.time()
divnet.CD.asa.ileum <- divnet(data,
                              X= "X5_ASA",
                              ncores = 24)
end.time.inf <- Sys.time()

#### ileum - Medication - Steroids ####
data <- CD.NI.ileum
data <- subset_samples(data,sample_data(data)$Steroids %in% c("0","1"))
data <- prune_taxa(taxa_sums(data) > 0, data)

start.time.inf <- Sys.time()
divnet.CD.str.ileum <- divnet(data,
                              X= "Steroids",
                              ncores = 24)
end.time.inf <- Sys.time()

#### ileum - Medication - TNF ####
data <- CD.NI.ileum
data <- subset_samples(data,sample_data(data)$TNF %in% c("0","1"))
data <- prune_taxa(taxa_sums(data) > 0, data)

start.time.inf <- Sys.time()
divnet.CD.tnf.ileum <- divnet(data,
                              X= "TNF",
                              ncores = 24)
end.time.inf <- Sys.time()

saveRDS(divnet.CD.imm.ileum,"CD_ileum_IMM_usage.rds")
saveRDS(divnet.CD.str.ileum,"CD_ileum_STR_usage.rds")
saveRDS(divnet.CD.tnf.ileum,"CD_ileum_TNF_usage.rds")
saveRDS(divnet.CD.asa.ileum,"CD_ileum_5ASA_usage.rds")



# #### ileum - Medication - TNF no immunomodulators ####
# data <- CD.NI.ileum
# data <- subset_samples(data,sample_data(data)$Immunomodulators == "0")
# data <- prune_taxa(taxa_sums(data) > 0, data)
# data <- subset_samples(data,sample_data(data)$TNF %in% c("0","1"))
# data <- prune_taxa(taxa_sums(data) > 0, data)
# 
# 
# start.time.inf <- Sys.time()
# divnet.CD.tnf.ileum <- divnet(data,
#                               X= "TNF",
#                               ncores = 24)
# end.time.inf <- Sys.time()
# 
# saveRDS(divnet.CD.tnf.ileum,"CD_ileum_TNF_usage_no_IMM_WashU.rds")

####Inflammation #####
data <- CD.colon.m
divnet.CD.colon.m <- divnet(data,
                            X= "Pathology",
                            ncores = 24)
data <- CD.ileum.m
divnet.CD.ileum.m <- divnet(data,
                            X= "Pathology",
                            ncores = 24)

saveRDS(divnet.CD.colon.m,"colon_match_inflammation.rds")
saveRDS(divnet.CD.ileum.m,"ileum_match_inflammation.rds")

# #Ileum 
# data <- CD.ileum.m
# estimates.t <-divnet.CD.ileum.m$shannon %>% summary %$% estimate
# ses.t <- sqrt(divnet.CD.ileum.m$`shannon-variance`)
# covar_matrix <- model.matrix(lmer(estimates.t ~  (1|sample_data(data)$PatientNo) + sample_data(data)$Pathology))
# results.i <- betta(estimates.t,
#                    ses.t,
#                    covar_matrix)

# 
 
# sample_data(CD.NI.ileum)$Recurrence2 <- ifelse(sample_data(CD.NI.ileum)$Rutgeerts_Score %in% c("i0","i1","i0-i1","i1-i2"),"no recurrence",
#                                                ifelse(sample_data(CD.NI.ileum)$Rutgeerts_Score %in%c("i2","i2+","i2-i3","i3","i4"),"recurrence"," "))
# 
# data <- CD.NI.colon
# estimates.t <-divnet.CD.colon.m$shannon %>% summary %$% estimate
# ses.t <- sqrt(divnet.CD.colon.m$`shannon-variance`)
# covar_matrix <- model.matrix(lmer(estimates.t ~ (1|sample_data(data)$PatientNo) + sample_data(data)$Pathology))
# results.c <- betta(estimates.t,
#                    ses.t,
#                    covar_matrix)
# results$table

#### Recurrence with i2 ####
sample_data(CD.NI.ileum)$Recurrence2 <- ifelse(sample_data(CD.NI.ileum)$Rutgeerts_Score %in% c("i0","i1","i0-i1","i1-i2"),"no recurrence",
                                               ifelse(sample_data(CD.NI.ileum)$Rutgeerts_Score %in%c("i2","i2+","i2-i3","i3","i4"),"recurrence"," "))

data <- CD.NI.ileum

data <- subset_samples(data,sample_data(data)$Recurrence2 %in% c("no recurrence","recurrence"))
data <- prune_taxa(taxa_sums(data) > 0, data)

start.time.inf <- Sys.time()
divnet.CD.ileum.i2 <- divnet(data,
                             X= "Recurrence2",
                             ncores = 24)
end.time.inf <- Sys.time()
saveRDS(divnet.CD.ileum.i2,"UNC_CD_ileum_recurrence.rds")


sample_data(CD.NI.colon)$Recurrence2 <- ifelse(sample_data(CD.NI.colon)$Rutgeerts_Score %in% c("i0","i1","i0-i1","i1-i2"),"no recurrence",
                                               ifelse(sample_data(CD.NI.colon)$Rutgeerts_Score %in%c("i2","i2+","i2-i3","i3","i4"),"recurrence"," "))

data <- CD.NI.colon

data <- subset_samples(data,sample_data(data)$Recurrence2 %in% c("no recurrence","recurrence"))
data <- prune_taxa(taxa_sums(data) > 0, data)

start.time.inf <- Sys.time()
divnet.CD.colon.i2 <- divnet(data,
                             X= "Recurrence2",
                             ncores = 24)
end.time.inf <- Sys.time()
saveRDS(divnet.CD.colon.i2,"UNC_CD_colon_recurrence.rds")