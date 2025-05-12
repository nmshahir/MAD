library(phyloseq)
library(magrittr)
library(DESeq2)
##### Make New File ####
setwd("C:/Users/nmshahir/Dropbox/IBD_microbiota/Analyses/Downstream_Analyses/Differential_Abundance/DiffAbund_4/reformatted_DESeq2/all_sites/Edited")

#New Taxonomy Matrix
tax_all_sites_format <- read.csv("C:/Users/nmshahir/Dropbox/IBD_microbiota/Analyses/Downstream_Analyses/Differential_Abundance/DiffAbund_4/reformatted_DESeq2/all_sites/Edited/all_sites_all_taxa_no_duplicates_20191007.csv",row.names = 1)
tax_all_sites_format <- as.matrix(tax_all_sites_format)
tax_all_sites_final = tax_table(tax_all_sites_format)

#New OTU Matrix
otu_all_sites_format <- read.csv("C:/Users/nmshahir/Dropbox/IBD_microbiota/Analyses/Downstream_Analyses/Differential_Abundance/DiffAbund_4/reformatted_DESeq2/all_sites/Edited/all_sites_all_otu_counts_no_duplicates_20191007.csv",row.names = 1)
colnames(otu_all_sites_format) <- gsub("X","",colnames(otu_all_sites_format))
colnames(otu_all_sites_format) <- gsub("\\."," ",colnames(otu_all_sites_format))

otu_all_sites_final = otu_table(otu_all_sites_format,taxa_are_rows = TRUE)

#New Meta Data 
all_sites_meta_test <- read.csv("C:/Users/nmshahir/Dropbox/IBD_microbiota/Analyses/Downstream_Analyses/Differential_Abundance/DiffAbund_4/reformatted_DESeq2/all_sites/Original/UNC_WashU_meta_20200310.csv", row.names = 1)
all_sites_meta_patients <- sample_data(all_sites_meta_test)

#make certain this is in data.frame format
#Merge new data
new_patients <- merge_phyloseq(otu_all_sites_final,tax_all_sites_final, all_sites_meta_patients)

ps.patients.g <- new_patients

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

#ileal disease
ileum.patients <- subset_samples(ps.UNC.patients,sample_data(ps.UNC.patients)$AnatomSite == "ileum")

nonIBD.CD.ileum <- subset_samples(ileum.patients, sample_data(ileum.patients)$DiseaseStatus %in% c("nonIBD","CD"))
nonIBD.CD.NI.ileum <- subset_samples(nonIBD.CD.ileum,sample_data(nonIBD.CD.ileum)$Pathology == "NI")
nonIBD.CD.NI.ileum <- prune_taxa(taxa_sums(nonIBD.CD.NI.ileum)>0, nonIBD.CD.NI.ileum)

nonIBD.UC.ileum <- subset_samples(ileum.patients,sample_data(ileum.patients)$DiseaseStatus %in% c("nonIBD","UC"))
nonIBD.UC.NI.ileum <- subset_samples(nonIBD.UC.ileum,sample_data(nonIBD.UC.ileum)$Pathology == "NI")
nonIBD.UC.NI.ileum <- prune_taxa(taxa_sums(nonIBD.UC.NI.ileum)>0,nonIBD.UC.NI.ileum)

CD.UC.ileum <- subset_samples(ileum.patients,sample_data(ileum.patients)$DiseaseStatus %in% c("CD","UC"))
CD.UC.NI.ileum <- subset_samples(CD.UC.ileum, sample_data(CD.UC.ileum)$Pathology == "NI")
CD.UC.NI.ileum <- prune_taxa(taxa_sums(CD.UC.NI.ileum)>0, CD.UC.NI.ileum)

#colon disease
colon.patients <- subset_samples(ps.UNC.patients,sample_data(ps.UNC.patients)$AnatomSite == "colon")

nonIBD.CD.colon <- subset_samples(colon.patients, sample_data(colon.patients)$DiseaseStatus %in% c("nonIBD","CD"))
nonIBD.CD.NI.colon <- subset_samples(nonIBD.CD.colon,sample_data(nonIBD.CD.colon)$Pathology == "NI")
nonIBD.CD.NI.colon <- prune_taxa(taxa_sums(nonIBD.CD.NI.colon)>0, nonIBD.CD.NI.colon)

nonIBD.UC.colon <- subset_samples(colon.patients, sample_data(colon.patients)$DiseaseStatus %in% c("nonIBD","UC"))
nonIBD.UC.NI.colon <- subset_samples(nonIBD.UC.colon,sample_data(nonIBD.UC.colon)$Pathology == "NI")
nonIBD.UC.NI.colon <- prune_taxa(taxa_sums(nonIBD.UC.NI.colon)>0, nonIBD.UC.NI.colon)

CD.UC.colon <- subset_samples(colon.patients,sample_data(colon.patients)$DiseaseStatus %in% c("CD","UC"))
CD.UC.NI.colon <- subset_samples(CD.UC.colon, sample_data(CD.UC.colon)$Pathology == "NI")
CD.UC.NI.colon <- prune_taxa(taxa_sums(CD.UC.NI.colon)>0, CD.UC.NI.colon)

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
#### Differential Abundance Analyses ####
#### Ileum - UC vs nonIBD ####

data <- nonIBD.UC.NI.ileum

#run deseq2
infdds = phyloseq_to_deseq2(data, ~DiseaseStatus)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$DiseaseStatus <- relevel(infddsClean$DiseaseStatus, ref="nonIBD")

DESeq2::resultsNames(infddsClean)
res <- DESeq2::lfcShrink(infddsClean, coef=2)

#res = DESeq2::results(infdds)
res = res[order(res$padj, na.last=NA), ]

mcols(res, use.names=TRUE)

alpha = 0.05

sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))

res = cbind(as(res,"data.frame"),as(tax_table(data)[rownames(res), ], "matrix"))

#Save significant table
write.csv(sigtab, file = "UNC_UC_vs_nonIBD_BASELINE_ileum_signif.csv")

#Save entire table
write.csv(res, file = "UNC_UC_vs_nonIBD_BASELINE_ileum_all_taxa.csv")


####Colon - UC vs nonIBD  ####
data <- nonIBD.UC.NI.colon

#run deseq2
infdds = phyloseq_to_deseq2(data, ~DiseaseStatus)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$DiseaseStatus <- relevel(infddsClean$DiseaseStatus, ref="nonIBD")

DESeq2::resultsNames(infddsClean)
res <- DESeq2::lfcShrink(infddsClean, coef=2)

#res = DESeq2::results(infdds)
res = res[order(res$padj, na.last=NA), ]

mcols(res, use.names=TRUE)

alpha = 0.05

sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))

res = cbind(as(res,"data.frame"),as(tax_table(data)[rownames(res), ], "matrix"))

#Save significant table
write.csv(sigtab, file = "UNC_UC_vs_nonIBD_BASELINE_colon_signif.csv")

#Save entire table
write.csv(res, file = "UNC_UC_vs_nonIBD_BASELINE_colon_all_taxa.csv")

####Colon - UC vs CD ####
data <- CD.UC.NI.colon

#run deseq2
infdds = phyloseq_to_deseq2(data, ~DiseaseStatus)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$DiseaseStatus <- relevel(infddsClean$DiseaseStatus, ref="CD")

DESeq2::resultsNames(infddsClean)
res <- DESeq2::lfcShrink(infddsClean, coef=2)

#res = DESeq2::results(infdds)
res = res[order(res$padj, na.last=NA), ]

mcols(res, use.names=TRUE)

alpha = 0.05

sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))

res = cbind(as(res,"data.frame"),as(tax_table(data)[rownames(res), ], "matrix"))

#Save significant table
write.csv(sigtab, file = "UNC_UC_vs_CD_BASELINE_colon_signif.csv")

#Save entire table
write.csv(res, file = "UNC_UC_vs_CD_BASELINE_colon_all_taxa.csv")

####UC Colon - E2 vs E3 ####
data <- UC.NI.colon

#run deseq2
infdds = phyloseq_to_deseq2(data, ~ClassificationUC)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$ClassificationUC <- relevel(infddsClean$ClassificationUC, ref="E2")

DESeq2::resultsNames(infddsClean)
res <- DESeq2::lfcShrink(infddsClean, coef=2)

#res = DESeq2::results(infdds)
res = res[order(res$padj, na.last=NA), ]

mcols(res, use.names=TRUE)

alpha = 0.05

sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))

res = cbind(as(res,"data.frame"),as(tax_table(data)[rownames(res), ], "matrix"))

#Save significant table
write.csv(sigtab, file = "UNC_UC_E3_vs_E2_BASELINE_colon_signif.csv")

#Save entire table
write.csv(res, file = "UNC_UC_E3_vs_E2_BASELINE_colon_all_taxa.csv")

#### Aero Split - UC Colon ####
sample_data(UC.NI.colon)$AeroLvl <- ifelse(sample_data(UC.NI.colon)$PatientNo %in% c("127","41500","42700","43000","46","52"),"Low","Health")
data <- UC.NI.colon

#run deseq2
infdds = phyloseq_to_deseq2(data, ~AeroLvl)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$AeroLvl <- relevel(infddsClean$AeroLvl, ref="Health")

DESeq2::resultsNames(infddsClean)
res <- DESeq2::lfcShrink(infddsClean, coef=2)

#res = DESeq2::results(infdds)
res = res[order(res$padj, na.last=NA), ]

mcols(res, use.names=TRUE)

alpha = 0.05

sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))

res = cbind(as(res,"data.frame"),as(tax_table(data)[rownames(res), ], "matrix"))

#Save significant table
write.csv(sigtab, file = "UNC_UC_Low_vs_Health_aero_BASELINE_colon_signif.csv")

#Save entire table
write.csv(res, file = "UNC_UC_Low_vs_Health_aero_BASELINE_colon_all_taxa.csv")

#### MEDICATION - UC Colon #####

sample_data(UC.NI.colon)$Steroids <- as.factor(sample_data(UC.NI.colon)$Steroids)
sample_data(UC.NI.colon)$Immunomodulators <- as.factor(sample_data(UC.NI.colon)$Immunomodulators)
sample_data(UC.NI.colon)$TNF <- as.factor(sample_data(UC.NI.colon)$TNF)
sample_data(UC.NI.colon)$X5_ASA <- as.factor(sample_data(UC.NI.colon)$X5_ASA)
sample_data(UC.NI.colon)$Probiotic <- as.factor(sample_data(UC.NI.colon)$Probiotic)

#### TNF ####
data <- UC.NI.colon
data.TNF <- subset_samples(data,sample_data(data)$TNF %in% c("0","1"))
data.TNF <- prune_taxa(taxa_sums(data.TNF)>0,data.TNF)
data <- data.TNF

#run deseq2
infdds = phyloseq_to_deseq2(data, ~TNF)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$TNF <- relevel(infddsClean$TNF, ref="0")

DESeq2::resultsNames(infddsClean)
res <- DESeq2::lfcShrink(infddsClean, coef=2)

#res = DESeq2::results(infdds)
res = res[order(res$padj, na.last=NA), ]

mcols(res, use.names=TRUE)

alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))
res = cbind(as(res, "data.frame"), as(tax_table(data)[rownames(res), ], "matrix"))

unname(head(sigtab))
write.csv(sigtab, file = "UNC_UC_colon_TNF_usage_vs_no_usage_BASELINE_signif.csv")
write.csv(res, file = "UNC_UC_colon_TNF_usage_vs_no_usage_BASELINE_all_taxa.csv")

#### Immunomodulators ####
data <- UC.NI.colon
data.IMM <- subset_samples(data,sample_data(data)$Immunomodulators %in% c("0","1"))
data.IMM <- prune_taxa(taxa_sums(data.IMM)>0,data.IMM)
data <- data.IMM

#run deseq2
infdds = phyloseq_to_deseq2(data, ~Immunomodulators)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Immunomodulators <- relevel(infddsClean$Immunomodulators, ref="0")

DESeq2::resultsNames(infddsClean)
res <- DESeq2::lfcShrink(infddsClean, coef=2)

#res = DESeq2::results(infdds)
res = res[order(res$padj, na.last=NA), ]

mcols(res, use.names=TRUE)

alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))
res = cbind(as(res, "data.frame"), as(tax_table(data)[rownames(res), ], "matrix"))

unname(head(sigtab))
write.csv(sigtab, file = "UNC_UC_colon_IMM_usage_vs_no_usage_BASELINE_signif.csv")
write.csv(res, file = "UNC_UC_colon_IMM_usage_vs_no_usage_BASELINE_all_taxa.csv")

#### Steroids  ####
data <- UC.NI.colon
data.STR <- subset_samples(data,sample_data(data)$Steroids %in% c("0","1"))
data.STR <- prune_taxa(taxa_sums(data.STR)>0,data.STR)
data <- data.STR

#run deseq2
infdds = phyloseq_to_deseq2(data, ~Steroids)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Steroids <- relevel(infddsClean$Steroids, ref="0")

DESeq2::resultsNames(infddsClean)
res <- DESeq2::lfcShrink(infddsClean, coef=2)

#res = DESeq2::results(infdds)
res = res[order(res$padj, na.last=NA), ]

mcols(res, use.names=TRUE)

alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))
res = cbind(as(res, "data.frame"), as(tax_table(data)[rownames(res), ], "matrix"))

unname(head(sigtab))
write.csv(sigtab, file = "UNC_UC_colon_STR_usage_vs_no_usage_BASELINE_signif.csv")
write.csv(res, file = "UNC_UC_colon_STR_usage_vs_no_usage_BASELINE_all_taxa.csv")


#### 5-ASA  ####
data <- UC.NI.colon
data.ASA <- subset_samples(data,sample_data(data)$X5_ASA %in% c("0","1"))
data.ASA <- prune_taxa(taxa_sums(data.ASA)>0,data.ASA)
data <- data.ASA

#run deseq2
infdds = phyloseq_to_deseq2(data, ~X5_ASA)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$X5_ASA <- relevel(infddsClean$X5_ASA, ref="0")

DESeq2::resultsNames(infddsClean)
res <- DESeq2::lfcShrink(infddsClean, coef=2)

#res = DESeq2::results(infdds)
res = res[order(res$padj, na.last=NA), ]

mcols(res, use.names=TRUE)

alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))
res = cbind(as(res, "data.frame"), as(tax_table(data)[rownames(res), ], "matrix"))

unname(head(sigtab))
write.csv(sigtab, file = "UNC_UC_colon_5ASA_usage_vs_no_usage_BASELINE_signif.csv")
write.csv(res, file = "UNC_UC_colon_5ASA_usage_vs_no_usage_BASELINE_all_taxa.csv")


#### Probiotics  ####
data <- UC.NI.colon
data.PRO <- subset_samples(data,sample_data(data)$Probiotic %in% c("0","1"))
data.PRO <- prune_taxa(taxa_sums(data.PRO)>0,data.PRO)
data <- data.PRO

#run deseq2
infdds = phyloseq_to_deseq2(data, ~Probiotic)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Probiotic <- relevel(infddsClean$Probiotic, ref="0")

DESeq2::resultsNames(infddsClean)
res <- DESeq2::lfcShrink(infddsClean, coef=2)

#res = DESeq2::results(infdds)
res = res[order(res$padj, na.last=NA), ]

mcols(res, use.names=TRUE)

alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))
res = cbind(as(res, "data.frame"), as(tax_table(data)[rownames(res), ], "matrix"))

unname(head(sigtab))
write.csv(sigtab, file = "UNC_UC_colon_PRO_usage_vs_no_usage_BASELINE_signif.csv")
write.csv(res, file = "UNC_UC_colon_PRO_usage_vs_no_usage_BASELINE_all_taxa.csv")



