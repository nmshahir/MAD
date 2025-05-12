#### Import required Libraries####
library("vegan")
library("MASS")
library("phyloseq")
library("ggplot2")
library("dada2")
library("ggridges")
library("ggpubr")
library("gridExtra")
####Pre-processing ####
# Import data files
#If on Windows
setwd("C:/Users/nmshahir/Dropbox/IBD_microbiota/Analyses/Metadata")

#If on MacOS
setwd("~/Dropbox/IBD_microbiota/Analyses/Metadata")
st.all <- readRDS("RDS/SHARE_no_chimera_20180925.rds")
taxa.all <- readRDS("RDS/tax_rdp_20180925.rds")
#guessmap = import_qiime_sample_data("Microbiome_Map_20180926.txt")
guessmap = import_qiime_sample_data("Microbiome_Map_20190713.txt")
test.tree <- read_tree("alignment_20181102.tre")

#Form phyloseq object
ps.all<- phyloseq(tax_table(taxa.all), guessmap,otu_table(st.all, taxa_are_rows = FALSE),phy_tree(test.tree))

#Remove all samples that do not have a tissue assigned
ps.tissue <- subset_samples(ps.all, sample_data(ps.all)$AnatomSite %in% c("colon","ileum"))
ps.tissue <- prune_taxa(taxa_sums(ps.tissue) > 0, ps.tissue)

ps.neg <- subset_samples(ps.all,
                         sample_data(ps.all)$SampleID %in% c("NC_221","UNCNEG1_221","UNCNEG2_221",
                                                             "UNCNEG3_221","NEG2_118", "NEG3_118",
                                                             "NEG4_118","NEG1_120","NEG2_120",
                                                             "NEG3_120","NEG4_120","NEG5_120",
                                                             "NEG6_120","NEG7_120","NEG8_120",
                                                             "NEG9_120","NEG2_106","NEG1_117"))
ps.neg <- prune_taxa(taxa_sums(ps.neg)>0, ps.neg)

ps.neg.221 <- subset_samples(ps.neg,sample_data(ps.neg)$SampleID %in% c("NC_221","UNCNEG1_221","UNCNEG2_221","UNCNEG3_221"))
ps.neg.221 <- prune_taxa(taxa_sums(ps.neg.221)>0, ps.neg.221)

ps.pos <- subset_samples(ps.all,
                         sample_data(ps.all)$SampleID %in% c("PA_221","PB_221"))
ps.pos <- prune_taxa(taxa_sums(ps.pos)>0, ps.pos)

#Remove all taxa that do not have a phyla assigned
#28589 RSV before
#28049 RSVs after
#540 RSVs removed
ps0 <- subset_taxa(ps.tissue, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# Get idea of phyla breakdown in samples
table(tax_table(ps0)[, "Phylum"], exclude = NULL)

####Institution filtering####
ps.UNC <- subset_samples(ps0,sample_data(ps0)$Cohort == "UNC")
ps.UNC <- prune_taxa(taxa_sums(ps.UNC)>0,ps.UNC)

ps.Chicago <- subset_samples(ps0,sample_data(ps0)$Cohort == "Chicago")
ps.Chicago <- prune_taxa(taxa_sums(ps.Chicago)>0,ps.Chicago)

ps.WashU <- subset_samples(ps0,sample_data(ps0)$Cohort == "WashU")
ps.WashU <- prune_taxa(taxa_sums(ps.WashU)>0,ps.WashU)

ps.Cedar <- subset_samples(ps0,sample_data(ps0)$Cohort == "Cedar")
ps.Cedar <- prune_taxa(taxa_sums(ps.Cedar)>0,ps.Cedar)

#Get taxa at 15% prevalence in particular institution
taxa.15.UNC <- getThresholdTaxa(ps.UNC, 0.15)
taxa.15.Chicago <- getThresholdTaxa(ps.Chicago, 0.15)
taxa.15.WashU <- getThresholdTaxa(ps.WashU, 0.15)
taxa.15.Cedar <- getThresholdTaxa(ps.Cedar, 0.15)

# Get all taxa at a particular institution
taxa.UNC <- getThresholdTaxa(ps.UNC, 0)
taxa.Chicago <- getThresholdTaxa(ps.Chicago, 0)
taxa.WashU <- getThresholdTaxa(ps.WashU, 0)
taxa.Cedar <- getThresholdTaxa(ps.Cedar, 0)

#Get comparison sets
# A - UNC
# B - Chicago
# C - WashU
# D - Cedar

taxa.all <- Reduce(union,list(taxa.UNC,taxa.Chicago,taxa.WashU, taxa.Cedar))
taxa.ABC <- Reduce(union,list(taxa.UNC,taxa.Chicago,taxa.WashU)) 
taxa.ABD <- Reduce(union,list(taxa.UNC,taxa.Chicago,taxa.Cedar))
taxa.ACD <- Reduce(union,list(taxa.UNC,taxa.WashU,taxa.Cedar))
taxa.BCD <- Reduce(union,list(taxa.Chicago,taxa.WashU,taxa.Cedar))


#Get taxa at 15% level that are specific to a single Institution
# 2 RSV specific to UNC
# 0 RSV specific to Chicago
# 5 RSV specific to WashU
# 1 RSV specific to Cedar
taxa.15.specific.UNC <- setdiff(taxa.15.UNC, taxa.BCD)
taxa.15.specific.Chicago <- setdiff(taxa.15.Chicago, taxa.ACD)
taxa.15.specific.WashU <- setdiff(taxa.15.WashU, taxa.ABD)
taxa.15.specific.Cedar <- setdiff(taxa.15.Cedar,taxa.ABC)

#### Get combined list of institution specific "contaminants" that are being removed
taxa.15.specific.all.sites <- Reduce(union,list(taxa.15.specific.UNC,taxa.15.specific.Chicago,taxa.15.specific.WashU,taxa.15.specific.Cedar))
#### Get list of taxa that are not in that list
taxa.no.site.contaminants <- taxa.all[!(taxa.all %in% taxa.15.specific.all.sites)]

####Batch filtering ####
ps.062 <- subset_samples(ps0,sample_data(ps0)$Batch == "seq062")
ps.062 <- prune_taxa(taxa_sums(ps.062)>0,ps.062)

ps.099 <- subset_samples(ps0,sample_data(ps0)$Batch == "seq099")
ps.099 <- prune_taxa(taxa_sums(ps.099)>0,ps.099)

ps.105 <- subset_samples(ps0,sample_data(ps0)$Batch == "seq105")
ps.105 <- prune_taxa(taxa_sums(ps.105)>0,ps.105)

ps.106 <- subset_samples(ps0,sample_data(ps0)$Batch == "seq106")
ps.106 <- prune_taxa(taxa_sums(ps.106)>0,ps.106)

ps.117 <- subset_samples(ps0,sample_data(ps0)$Batch == "seq117")
ps.117 <- prune_taxa(taxa_sums(ps.117)>0,ps.117)

ps.118 <- subset_samples(ps0,sample_data(ps0)$Batch == "seq118")
ps.118 <- prune_taxa(taxa_sums(ps.118)>0,ps.118)

ps.120 <- subset_samples(ps0,sample_data(ps0)$Batch == "seq120")
ps.120 <- prune_taxa(taxa_sums(ps.120)>0,ps.120)

ps.221 <- subset_samples(ps0,sample_data(ps0)$Batch == "seq221")
ps.221 <- prune_taxa(taxa_sums(ps.221)>0,ps.221)

####Extra Cleaning on Batch 221

#Remove contam RSVs
list.01.taxa.bad <- c(taxa_names(ps.221)[737],taxa_names(ps.221)[378],taxa_names(ps.221)[3476],taxa_names(ps.221)[727],taxa_names(ps.221)[3630],taxa_names(ps.221)[3483],taxa_names(ps.221)[3472],taxa_names(ps.221)[391],taxa_names(ps.221)[746],taxa_names(ps.221)[3712],taxa_names(ps.221)[3651],taxa_names(ps.221)[3749],taxa_names(ps.221)[714],taxa_names(ps.221)[405],taxa_names(ps.221)[3529],taxa_names(ps.221)[423],taxa_names(ps.221)[716],taxa_names(ps.221)[3890],taxa_names(ps.221)[3631],taxa_names(ps.221)[3489],taxa_names(ps.221)[3939],taxa_names(ps.221)[726],taxa_names(ps.221)[3659],taxa_names(ps.221)[3658],taxa_names(ps.221)[3434],taxa_names(ps.221)[2921],taxa_names(ps.221)[721],taxa_names(ps.221)[732],taxa_names(ps.221)[713],taxa_names(ps.221)[2841],taxa_names(ps.221)[699],taxa_names(ps.221)[3892],taxa_names(ps.221)[3936],taxa_names(ps.221)[743],taxa_names(ps.221)[3633],taxa_names(ps.221)[3189],taxa_names(ps.221)[3407],taxa_names(ps.221)[730],taxa_names(ps.221)[3883],taxa_names(ps.221)[5865],taxa_names(ps.221)[3672],taxa_names(ps.221)[720],taxa_names(ps.221)[3402],taxa_names(ps.221)[5836])
taxa.221 <- getThresholdTaxa(ps.221,0)
taxa.221.clean <- setdiff(taxa.221,list.01.taxa.bad)
ps.221.clean <- prune_taxa(taxa.221.clean,ps.221)

ps.221 <- ps.221.clean

#Remove contam samples
#### Removed identified "contaminant" samples
contam.samples.list <- c("17601_221",
                         "24001_221",
                         "20902_221",
                         "4404_221",
                         "5021_221",
                         "31601_221",
                         "19901_221",
                         "27601_221",
                         "21101_221",
                         "20711_221",
                         "22111_221",
                         "22101_221",
                         "20703_221",
                         "30001_221",
                         "15901_221",
                         "21108_221")
ps.221.clean.samples <- subset_samples(ps.221,!(sample_data(ps.221)$SampleID %in% contam.samples.list))
ps.221.clean.samples <- prune_taxa(taxa_sums(ps.221.clean.samples) > 0, ps.221.clean.samples)

ps.221 <- ps.221.clean.samples
#Get taxa at 15% prevalence in particular batch
taxa.15.062 <- getThresholdTaxa(ps.062, 0.15)
taxa.15.099 <- getThresholdTaxa(ps.099, 0.15)
taxa.15.105 <- getThresholdTaxa(ps.105, 0.15)
taxa.15.106 <- getThresholdTaxa(ps.106, 0.15)
taxa.15.117 <- getThresholdTaxa(ps.117, 0.15)
taxa.15.118 <- getThresholdTaxa(ps.118, 0.15)
taxa.15.120 <- getThresholdTaxa(ps.120, 0.15)
taxa.15.221 <- getThresholdTaxa(ps.221, 0.15)

# Get all taxa in particular batch
taxa.062 <- getThresholdTaxa(ps.062, 0)
taxa.099 <- getThresholdTaxa(ps.099, 0)
taxa.105 <- getThresholdTaxa(ps.105, 0)
taxa.106 <- getThresholdTaxa(ps.106, 0)
taxa.117 <- getThresholdTaxa(ps.117, 0)
taxa.118 <- getThresholdTaxa(ps.118, 0)
taxa.120 <- getThresholdTaxa(ps.120, 0)
taxa.221 <- getThresholdTaxa(ps.221, 0)

#Get comparison sets
# A - 062
# B - 099
# C - 105
# D - 106
# E - 117
# F - 118
# G - 120
# H - 221
taxa.all <- Reduce(union,list(taxa.062,taxa.099,taxa.105,taxa.106,taxa.117,taxa.118, taxa.120, taxa.221))
taxa.ABCDEFG <- Reduce(union,list(taxa.062,taxa.099,taxa.105,taxa.106,taxa.117,taxa.118,taxa.120))
taxa.ABCDEFH <- Reduce(union,list(taxa.062,taxa.099,taxa.105,taxa.106,taxa.117,taxa.118,taxa.221))
taxa.ABCDEGH <- Reduce(union,list(taxa.062,taxa.099,taxa.105,taxa.106,taxa.117,taxa.120,taxa.221))
taxa.ABCDFGH <- Reduce(union,list(taxa.062,taxa.099,taxa.105,taxa.106,taxa.120,taxa.118,taxa.221))
taxa.ABCEFGH <- Reduce(union,list(taxa.062,taxa.099,taxa.105,taxa.120,taxa.117,taxa.118,taxa.221))
taxa.ABDEFGH <- Reduce(union,list(taxa.062,taxa.099,taxa.120,taxa.106,taxa.117,taxa.118,taxa.221))
taxa.ACDEFGH <- Reduce(union,list(taxa.062,taxa.120,taxa.105,taxa.106,taxa.117,taxa.118,taxa.221))
taxa.BCDEFGH <- Reduce(union,list(taxa.120,taxa.099,taxa.105,taxa.106,taxa.117,taxa.118,taxa.221))

#Get taxa at 15% level that are specific to a single batch
taxa.15.specific.062 <- setdiff(taxa.15.062, taxa.BCDEFGH)
taxa.15.specific.099 <- setdiff(taxa.15.099, taxa.ACDEFGH)
taxa.15.specific.105 <- setdiff(taxa.15.105, taxa.ABDEFGH)
taxa.15.specific.106 <- setdiff(taxa.15.106, taxa.ABCEFGH)
taxa.15.specific.117 <- setdiff(taxa.15.117, taxa.ABCDFGH)
taxa.15.specific.118 <- setdiff(taxa.15.118, taxa.ABCDEGH)
taxa.15.specific.120 <- setdiff(taxa.15.120, taxa.ABCDEFH)
taxa.15.specific.221 <- setdiff(taxa.15.221, taxa.ABCDEFG)

#### Get combined list of batch specific "contaminants" that are being removed
taxa.15.specific.all.batch <- Reduce(union,list(taxa.15.specific.062,taxa.15.specific.099,taxa.15.specific.105,taxa.15.specific.106,taxa.15.specific.117,taxa.15.specific.118,taxa.15.specific.120,taxa.15.specific.221))

#### Get list of taxa that are not in that list
taxa.no.batch.contaminants <- taxa.all[!(taxa.all %in% taxa.15.specific.all.batch)]

#### All contaminants
taxa.specific.15.all.contaminants <- union(taxa.15.specific.all.batch,taxa.15.specific.all.sites)

#### Remove all contaminants
taxa.non.contaminants <- taxa.all[!(taxa.all %in% taxa.specific.15.all.contaminants)]

ps1 <- prune_taxa(taxa.non.contaminants,ps0)
ps1 <- prune_samples(sample_sums(ps1)>0,ps1)

ps1.no.contam.samples <- subset_samples(ps1,!(sample_data(ps1)$SampleID %in% contam.samples.list))
ps1.no.contam.samples <- prune_taxa(taxa_sums(ps1.no.contam.samples) > 0, ps1.no.contam.samples)
taxa.no.contam.samples <- getThresholdTaxa(ps1.no.contam.samples,0)
taxa.no.contam.rsv <- setdiff(taxa.no.contam.samples,list.01.taxa.bad)

ps1.no.contam.samples.no.contam.rsv <- prune_taxa(taxa.no.contam.rsv,ps1.no.contam.samples)
ps1 <- ps1.no.contam.samples.no.contam.rsv

saveRDS(ps1,"cleaned_data_no_merge_20190713.rds")

##### Step 2 ####
sample_data(ps1)$ID <- paste(sample_data(ps1)$PatientNo,sample_data(ps1)$DiseaseStatus,sample_data(ps1)$Pathology,sample_data(ps1)$AnatomSite,sample_data(ps1)$Cohort,sep=" ")

#Merge Samples if any
ps.patients <- merge_samples(ps1,"ID")
sample_data(ps.patients)$DiseaseStatus <- as.integer(sample_data(ps.patients)$DiseaseStatus)
sample_data(ps.patients)$DiseaseStatus <- as.factor(sample_data(ps.patients)$DiseaseStatus)
levels(sample_data(ps.patients)$DiseaseStatus) <- levels(sample_data(ps1)$DiseaseStatus)

sample_data(ps.patients)$AnatomSite <- as.integer(sample_data(ps.patients)$AnatomSite)
sample_data(ps.patients)$AnatomSite <- as.factor(sample_data(ps.patients)$AnatomSite)
levels(sample_data(ps.patients)$AnatomSite) <- levels(sample_data(ps1)$AnatomSite)

sample_data(ps.patients)$PatientNo <- as.integer(sample_data(ps.patients)$PatientNo)
sample_data(ps.patients)$PatientNo <- as.factor(sample_data(ps.patients)$PatientNo)
levels(sample_data(ps.patients)$PatientNo) <- levels(sample_data(ps1)$PatientNo)

sample_data(ps.patients)$Pathology <- as.integer(sample_data(ps.patients)$Pathology )
sample_data(ps.patients)$Pathology  <- as.factor(sample_data(ps.patients)$Pathology )
levels(sample_data(ps.patients)$Pathology ) <- levels(sample_data(ps1)$Pathology )

sample_data(ps.patients)$Behavior <- as.integer(sample_data(ps.patients)$Behavior)
sample_data(ps.patients)$Behavior   <- as.factor(sample_data(ps.patients)$Behavior)
levels(sample_data(ps.patients)$Behavior  ) <- levels(sample_data(ps1)$Behavior  )

sample_data(ps.patients)$Location <- as.integer(sample_data(ps.patients)$Location)
sample_data(ps.patients)$Location   <- as.factor(sample_data(ps.patients)$Location)
levels(sample_data(ps.patients)$Location  ) <- levels(sample_data(ps1)$Location  )

sample_data(ps.patients)$Sex <- as.integer(sample_data(ps.patients)$Sex)
sample_data(ps.patients)$Sex   <- as.factor(sample_data(ps.patients)$Sex)
levels(sample_data(ps.patients)$Sex) <- levels(sample_data(ps1)$Sex)

sample_data(ps.patients)$Race <- as.integer(sample_data(ps.patients)$Race)
sample_data(ps.patients)$Race   <- as.factor(sample_data(ps.patients)$Race)
levels(sample_data(ps.patients)$Race) <- levels(sample_data(ps1)$Race)

sample_data(ps.patients)$SmokingStatus <- as.integer(sample_data(ps.patients)$SmokingStatus)
sample_data(ps.patients)$SmokingStatus   <- as.factor(sample_data(ps.patients)$SmokingStatus)
levels(sample_data(ps.patients)$SmokingStatus) <- levels(as.factor(sample_data(ps1)$SmokingStatus))

sample_data(ps.patients)$Antibiotics <- as.integer(sample_data(ps.patients)$Antibiotics)
sample_data(ps.patients)$Antibiotics   <- as.factor(sample_data(ps.patients)$Antibiotics)
levels(sample_data(ps.patients)$Antibiotics) <- levels(as.factor(sample_data(ps1)$Antibiotics))

sample_data(ps.patients)$X5_ASA <- as.integer(sample_data(ps.patients)$X5_ASA)
sample_data(ps.patients)$X5_ASA   <- as.factor(sample_data(ps.patients)$X5_ASA)
levels(sample_data(ps.patients)$X5_ASA) <- levels(as.factor(sample_data(ps1)$X5_ASA))

sample_data(ps.patients)$Steroids <- as.integer(sample_data(ps.patients)$Steroids)
sample_data(ps.patients)$Steroids   <- as.factor(sample_data(ps.patients)$Steroids)
levels(sample_data(ps.patients)$Steroids) <- levels(as.factor(sample_data(ps1)$Steroids))

sample_data(ps.patients)$Immunomodulators <- as.integer(sample_data(ps.patients)$Immunomodulators)
sample_data(ps.patients)$Immunomodulators   <- as.factor(sample_data(ps.patients)$Immunomodulators)
levels(sample_data(ps.patients)$Immunomodulators) <- levels(as.factor(sample_data(ps1)$Immunomodulators))

sample_data(ps.patients)$TNF <- as.integer(sample_data(ps.patients)$TNF)
sample_data(ps.patients)$TNF  <- as.factor(sample_data(ps.patients)$TNF)
levels(sample_data(ps.patients)$TNF) <- levels(as.factor(sample_data(ps1)$TNF))

sample_data(ps.patients)$Cdiff <- as.integer(sample_data(ps.patients)$Cdiff)
sample_data(ps.patients)$Cdiff  <- as.factor(sample_data(ps.patients)$Cdiff)
levels(sample_data(ps.patients)$Cdiff) <- levels(as.factor(sample_data(ps1)$Cdiff))

sample_data(ps.patients)$Ethnicity <- as.integer(sample_data(ps.patients)$Ethnicity)
sample_data(ps.patients)$Ethnicity  <- as.factor(sample_data(ps.patients)$Ethnicity)
levels(sample_data(ps.patients)$Ethnicity) <- levels(as.factor(sample_data(ps1)$Ethnicity))

sample_data(ps.patients)$Cohort <- as.integer(sample_data(ps.patients)$Cohort)
sample_data(ps.patients)$Cohort<- as.factor(sample_data(ps.patients)$Cohort)
levels(sample_data(ps.patients)$Cohort) <- levels(as.factor(sample_data(ps1)$Cohort))

sample_data(ps.patients)$ClassificationUC <- as.integer(sample_data(ps.patients)$ClassificationUC)
sample_data(ps.patients)$ClassificationUC <- as.factor(sample_data(ps.patients)$ClassificationUC)
levels(sample_data(ps.patients)$ClassificationUC) <- levels(as.factor(sample_data(ps1)$ClassificationUC))

sample_data(ps.patients)$B3_v_nonB3 <- as.integer(sample_data(ps.patients)$B3_v_nonB3)
sample_data(ps.patients)$B3_v_nonB3 <- as.factor(sample_data(ps.patients)$B3_v_nonB3)
levels(sample_data(ps.patients)$B3_v_nonB3) <- levels(as.factor(sample_data(ps1)$B3_v_nonB3))

sample_data(ps.patients)$Perianal <- as.integer(sample_data(ps.patients)$Perianal)
sample_data(ps.patients)$Perianal <- as.factor(sample_data(ps.patients)$Perianal)
levels(sample_data(ps.patients)$Perianal) <- levels(as.factor(sample_data(ps1)$Perianal))

sample_data(ps.patients)$Subtype <- as.integer(sample_data(ps.patients)$Subtype)
sample_data(ps.patients)$Subtype <- as.factor(sample_data(ps.patients)$Subtype)
levels(sample_data(ps.patients)$Subtype) <- levels(as.factor(sample_data(ps1)$Subtype))

sample_data(ps.patients)$Age <- as.integer(sample_data(ps.patients)$Age)
sample_data(ps.patients)$Age <- as.factor(sample_data(ps.patients)$Age)
levels(sample_data(ps.patients)$Age) <- levels(as.factor(sample_data(ps1)$Age))

sample_data(ps.patients)$Antibiotics.Fine <- as.integer(sample_data(ps.patients)$Antibiotics.Fine)
sample_data(ps.patients)$Antibiotics.Fine <- as.factor(sample_data(ps.patients)$Antibiotics.Fine)
levels(sample_data(ps.patients)$Antibiotics.Fine) <- levels(as.factor(sample_data(ps1)$Antibiotics.Fine))

sample_data(ps.patients)$Rutgeerts_Score <- as.integer(sample_data(ps.patients)$Rutgeerts_Score )
sample_data(ps.patients)$Rutgeerts_Score <- as.factor(sample_data(ps.patients)$Rutgeerts_Score )
levels(sample_data(ps.patients)$Rutgeerts_Score) <- levels(as.factor(sample_data(ps1)$Rutgeerts_Score))

sample_data(ps.patients)$Region <- as.integer(sample_data(ps.patients)$Region )
sample_data(ps.patients)$Region <- as.factor(sample_data(ps.patients)$Region)
levels(sample_data(ps.patients)$Region) <- levels(as.factor(sample_data(ps1)$Region))

sample_data(ps.patients)$Depth <- sample_sums(ps.patients)
sample_data(ps.patients)$Recurrence <- ifelse(sample_data(ps.patients)$Rutgeerts_Score %in% c("i0","i1","i0-i1","i1-i2"),"no recurrence",
                                              ifelse(sample_data(ps.patients)$Rutgeerts_Score %in%c("i2-i3","i3","i4"),"recurrence"," ")
)

####STEP 3 ####
#adjusting the taxa file
tax.clean <- data.frame(tax_table(ps.patients))

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}

for (i in 1:nrow(tax.clean)){
  if (is.na(tax.clean[i,7])){
    tax.clean$Species[i] <- paste("Species", i, sep = "_")
  }
}
#All unknown genera are relabeled
for (i in 1:nrow(tax.clean)){
  if (is.na(tax.clean[i,6])){
    tax.clean$Genus[i] <- paste("Genus", i, sep = "_")
  }
}

tax_table(ps.patients) <- as.matrix(tax.clean)

ps.temp <- subset_samples(ps.patients,sample_data(ps.patients)$PatientNo != "ENCODE")
ps.temp <- prune_taxa(taxa_sums(ps.temp)>0,ps.temp)

ps.patients <- ps.temp

#This takes 11.8 hours
start_time.s <- Sys.time()
ps.patients.s <- tax_glom(ps.patients,"Species",NArm=FALSE)
end_time.s <- Sys.time()

#This takes 1.6 hours
start_time <- Sys.time()
ps.patients.g <- tax_glom(ps.patients,"Genus",NArm=FALSE)
end_time <- Sys.time()
ps.patients.old <- ps.patients


saveRDS(ps.patients.s,"ps_patients_species_20190713.rds")
saveRDS(ps.patients.g,"ps_patients_genus_20190713.rds")
saveRDS(ps.patients,"ps_patients_RSV_20190713.rds")
