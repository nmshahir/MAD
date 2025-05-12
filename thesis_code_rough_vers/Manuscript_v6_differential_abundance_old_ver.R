####Manuscript - Differential Abundance Analysis ####
ps.patients.g <- readRDS("ps_patients_species_20190713.rds")
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

CD.ileum.m <- subset_samples(CD.ileum,sample_data(CD.ileum)$PatientNo %in% c("1004098","129","137","138","144","15","19","40100","41","41700","42100","42400","42600","44","44000","47","51","56","57","59","60","70","77"))
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

#### Differential Abundance Analyses ####
####UNC Ileum ####

nonIBD.NI.ileum.taxa.25 <- getThresholdTaxa(nonIBD.NI.ileum,0.25)
CD.NI.ileum.taxa.25 <- getThresholdTaxa(CD.NI.ileum,0.25)

nonIBD.CD.NI.ileum.taxa.25 <- union(nonIBD.NI.ileum.taxa.25,CD.NI.ileum.taxa.25)

nonIBD.CD.NI.ileum.25 <- prune_taxa(nonIBD.CD.NI.ileum.taxa.25,nonIBD.CD.NI.ileum)

# data <- subset_samples(nonIBD.CD.NI.ileum.25 , sample_data(nonIBD.CD.NI.ileum.25)$Age != "None")
# data  <- prune_taxa(taxa_sums(data) > 0, data )

data <- nonIBD.CD.NI.ileum.25
#data <- nonIBD.CD.NI.ileum
#lose 21 samples 

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
write.csv(sigtab, file = "UNC_nonIBD_vs_CD_BASELINE_ileum_25_stratified_species.csv")

#Save entire table
write.csv(res, file = "UNC_nonIBD_vs_CD_BASELINE_ileum_25_stratified_species_all_taxa.csv")


####UNC Colon ####
nonIBD.NI.colon <- subset_samples(nonIBD.CD.NI.colon,sample_data(nonIBD.CD.NI.colon)$DiseaseStatus=="nonIBD")
nonIBD.NI.colon <- prune_taxa(taxa_sums(nonIBD.NI.colon)>0,nonIBD.NI.colon)

CD.NI.colon <- subset_samples(nonIBD.CD.NI.colon,sample_data(nonIBD.CD.NI.colon)$DiseaseStatus=="CD")
CD.NI.colon <- prune_taxa(taxa_sums(CD.NI.colon)>0,CD.NI.colon)

nonIBD.NI.colon.taxa.25 <- getThresholdTaxa(nonIBD.NI.colon,0.25)
CD.NI.colon.taxa.25 <- getThresholdTaxa(CD.NI.colon,0.25)

nonIBD.CD.NI.colon.taxa.25 <- union(nonIBD.NI.colon.taxa.25,CD.NI.colon.taxa.25)

nonIBD.CD.NI.colon.25 <- prune_taxa(nonIBD.CD.NI.colon.taxa.25,nonIBD.CD.NI.colon)

data <- nonIBD.CD.NI.colon.25

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
res = cbind(as(res, "data.frame"), as(tax_table(data)[rownames(res), ], "matrix"))

unname(head(sigtab))
write.csv(sigtab, file = "UNC_nonIBD_vs_CD_BASELINE_colon_25_stratified_species.csv")
write.csv(res, file = "UNC_nonIBD_vs_CD_BASELINE_colon_25_stratified_species_all_taxa.csv")

#### nonIBD tissue ####
nonIBD.NI.tissue.taxa.25 <- union(nonIBD.NI.colon.taxa.25,nonIBD.NI.ileum.taxa.25)

nonIBD.NI.tissue.25 <- prune_taxa(nonIBD.NI.tissue.taxa.25,nonIBD.NI.patients)

data <- nonIBD.NI.tissue.25
#run deseq2
infdds = phyloseq_to_deseq2(data, ~AnatomSite)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$AnatomSite<- relevel(infddsClean$AnatomSite, ref="colon")

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
write.csv(sigtab, file = "UNC_nonIBD_tissue_colon_BASELINE_25.csv")
write.csv(res, file = "UNC_nonIBD_tissue_colon_BASELINE_25_all_taxa.csv")

#### CD tissue ####

CD.NI.tissue.taxa.25 <- union(CD.NI.colon.taxa.25,CD.NI.ileum.taxa.25)

CD.NI.tissue.25 <- prune_taxa(CD.NI.tissue.taxa.25,CD.NI.patients)

data <- CD.NI.tissue.25


#lose 21 samples 

#run deseq2
infdds = phyloseq_to_deseq2(data, ~AnatomSite)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$AnatomSite<- relevel(infddsClean$AnatomSite, ref="colon")

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
write.csv(sigtab, file = "UNC_CD_tissue_colon_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_tissue_colon_BASELINE_25_all_taxa.csv")

#### CD B2 and B3 tissue ####

CD.NI.tissue.taxa.25 <- union(CD.NI.colon.taxa.25,CD.NI.ileum.taxa.25)

CD.NI.tissue.25 <- prune_taxa(CD.NI.tissue.taxa.25,CD.NI.patients)

data <- CD.NI.tissue.25

data.b23 <- subset_samples(data,sample_data(data)$Behavior %in% c("B2","B3"))
data.b23 <- prune_taxa(taxa_sums(data.b23)>0,data.b23)

data <- data.b23

#lose 21 samples 

#run deseq2
infdds = phyloseq_to_deseq2(data, ~AnatomSite)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$AnatomSite<- relevel(infddsClean$AnatomSite, ref="colon")

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
write.csv(sigtab, file = "UNC_CD_B2_B3_tissue_colon_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_B2_B3_tissue_colon_BASELINE_25_all_taxa.csv")

#### CD B2 - tissue all ####
data <- CD.NI.tissue.25
data.b2 <- subset_samples(data,sample_data(data)$Behavior == "B2")
data.b2 <- prune_taxa(taxa_sums(data.b2)>0,data.b2)

data <- data.b2

#run deseq2
infdds = phyloseq_to_deseq2(data, ~AnatomSite)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$AnatomSite<- relevel(infddsClean$AnatomSite, ref="colon")

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
write.csv(sigtab, file = "UNC_CD_B2_tissue_colon_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_B2_tissue_colon_BASELINE_25_all_taxa.csv")

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus","Species")]


theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Species))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$GS <- paste(sigtabgen$Genus,sigtabgen$Species,sep=" ")
ggplot(sigtabgen, aes(y=GS, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "black", size = 0.2) +
  geom_point(size=8) + 
  scale_color_manual(values = microbiota.palette) +
  ggtitle("UNC CD B2 (p = 0.05)") +
  theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5)) +
  xlab(expression(paste(log[2],"Fold Change")))

#### CD B2 - tissue - concordant ####
data <- CD.NI.tissue.25
data.b <- subset_samples(data,sample_data(data)$Behavior == "B2")
data.b <- prune_taxa(taxa_sums(data.b)>0,data.b)

data <- data.b

data.concordant <- subset_samples(data,!(sample_data(data)$PatientNo %in% c("92","96")))
data.concordant <- prune_taxa(taxa_sums(data.concordant) > 0, data.concordant)
data <- data.concordant

data.concordant <- subset_samples(data,sample_names(data) != "51 CD NI ileum UNC")
data.concordant <- prune_taxa(taxa_sums(data.concordant) > 0, data.concordant)

data <- data.concordant

#run deseq2
infdds = phyloseq_to_deseq2(data, ~AnatomSite)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$AnatomSite<- relevel(infddsClean$AnatomSite, ref="colon")

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
write.csv(sigtab, file = "UNC_CD_B2_concordant_colon_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_B2_concordcant_colon_BASELINE_25_all_taxa.csv")

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus","Species")]


theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Species))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$GS <- paste(sigtabgen$Genus,sigtabgen$Species,sep=" ")
ggplot(sigtabgen, aes(y=GS, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "black", size = 0.2) +
  geom_point(size=8) + 
  scale_color_manual(values = microbiota.palette) +
  ggtitle("UNC CD B2 - Concordant Strictures(p = 0.05)") +
  theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5)) +
  xlab(expression(paste(log[2],"Fold Change")))

#### CD B3 - tissue all ####
data <- CD.NI.tissue.25
data.b3 <- subset_samples(data,sample_data(data)$Behavior == "B3")
data.b3 <- prune_taxa(taxa_sums(data.b3)>0,data.b3)

data <- data.b3

#run deseq2
infdds = phyloseq_to_deseq2(data, ~AnatomSite)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$AnatomSite<- relevel(infddsClean$AnatomSite, ref="colon")

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
write.csv(sigtab, file = "UNC_CD_B3_tissue_colon_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_B3_tissue_colon_BASELINE_25_all_taxa.csv")

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus","Species")]


theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Species))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$GS <- paste(sigtabgen$Genus,sigtabgen$Species,sep=" ")
ggplot(sigtabgen, aes(y=GS, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "black", size = 0.2) +
  geom_point(size=8) + 
  scale_color_manual(values = microbiota.palette) +
  ggtitle("UNC CD B3 (p = 0.05)") +
  theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5)) +
  xlab(expression(paste(log[2],"Fold Change")))
#### CD Colon - Inflammation | NOTHING ####

CD.colon.m.NI <- subset_samples(CD.colon.m,sample_data(CD.colon.m)$Pathology == "NI")
CD.colon.m.NI <- prune_taxa(taxa_sums(CD.colon.m.NI)>0, CD.colon.m.NI)

CD.colon.m.I <- subset_samples(CD.colon.m,sample_data(CD.colon.m)$Pathology == "I")
CD.colon.m.I <- prune_taxa(taxa_sums(CD.colon.m.I)>0, CD.colon.m.I)

CD.colon.m.NI.taxa.25 <- getThresholdTaxa(CD.colon.m.NI,0.25)
CD.colon.m.I.taxa.25 <- getThresholdTaxa(CD.colon.m.I,0.25)

CD.colon.m.taxa.25 <- union(CD.colon.m.NI.taxa.25,CD.colon.m.I.taxa.25)

CD.colon.m.25 <- prune_taxa(CD.colon.m.taxa.25,CD.colon.m)

data <- CD.colon.m.25
sample_data(data)$PatientNo <- as.factor(sample_data(data)$PatientNo)

#run deseq2
infdds = phyloseq_to_deseq2(data, ~PatientNo + Pathology)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Pathology <- relevel(infddsClean$Pathology, ref="I")

DESeq2::resultsNames(infddsClean)
res <- DESeq2::lfcShrink(infddsClean, coef=4)

#res = DESeq2::results(infdds)
res = res[order(res$padj, na.last=NA), ]

mcols(res, use.names=TRUE)

alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))

unname(head(sigtab))
write.csv(sigtab, file = "UNC_CD_colon_inflamed_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_colon_inflamed_BASELINE_25_all_taxa.csv")

#### CD Ileum - Inflammation  ####

CD.ileum.m.NI <- subset_samples(CD.ileum.m,sample_data(CD.ileum.m)$Pathology == "NI")
CD.ileum.m.NI <- prune_taxa(taxa_sums(CD.ileum.m.NI)>0, CD.ileum.m.NI)

CD.ileum.m.I <- subset_samples(CD.ileum.m,sample_data(CD.ileum.m)$Pathology == "I")
CD.ileum.m.I <- prune_taxa(taxa_sums(CD.ileum.m.I)>0, CD.ileum.m.I)

CD.ileum.m.NI.taxa.25 <- getThresholdTaxa(CD.ileum.m.NI,0.25)
CD.ileum.m.I.taxa.25 <- getThresholdTaxa(CD.ileum.m.I,0.25)

CD.ileum.m.taxa.25 <- union(CD.ileum.m.NI.taxa.25,CD.ileum.m.I.taxa.25)

CD.ileum.m.25 <- prune_taxa(CD.ileum.m.taxa.25,CD.ileum.m)

data <- CD.ileum.m.25
sample_data(data)$PatientNo <- as.factor(sample_data(data)$PatientNo)

#run deseq2
infdds = phyloseq_to_deseq2(data, ~PatientNo + Pathology)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Pathology <- relevel(infddsClean$Pathology, ref="I")

DESeq2::resultsNames(infddsClean)
res <- DESeq2::lfcShrink(infddsClean, coef=16)

#res = DESeq2::results(infdds)
res = res[order(res$padj, na.last=NA), ]

mcols(res, use.names=TRUE)

alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(data)[rownames(sigtab), ], "matrix"))

unname(head(sigtab))
write.csv(sigtab, file = "UNC_CD_ileum_inflamed_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_ileum_inflamed_BASELINE_25_all_taxa.csv")

#### CD Colon - Behavior####

CD.NI.colon.B1 <- subset_samples(CD.NI.colon,sample_data(CD.NI.colon)$Behavior == "B1")
CD.NI.colon.B1 <- prune_taxa(taxa_sums(CD.NI.colon.B1)>0, CD.NI.colon.B1)

CD.NI.colon.B2 <- subset_samples(CD.NI.colon,sample_data(CD.NI.colon)$Behavior == "B2")
CD.NI.colon.B2 <- prune_taxa(taxa_sums(CD.NI.colon.B2)>0, CD.NI.colon.B2)

CD.NI.colon.B3 <- subset_samples(CD.NI.colon,sample_data(CD.NI.colon)$Behavior == "B3")
CD.NI.colon.B3 <- prune_taxa(taxa_sums(CD.NI.colon.B3)>0, CD.NI.colon.B3)

CD.NI.colon.B1.taxa.25 <- getThresholdTaxa(CD.NI.colon.B1,0.25)
CD.NI.colon.B2.taxa.25 <- getThresholdTaxa(CD.NI.colon.B2,0.25)
CD.NI.colon.B3.taxa.25 <- getThresholdTaxa(CD.NI.colon.B3,0.25)

CD.colon.B.taxa.25 <- union(union(CD.NI.colon.B1.taxa.25,CD.NI.colon.B2.taxa.25),CD.NI.colon.B3.taxa.25)

CD.colon.B.25 <- prune_taxa(CD.colon.B.taxa.25,CD.NI.colon)


#### CD Colon - B2 (all) vs B1

data <- CD.colon.B.25

data.B12 <- subset_samples(data,sample_data(data)$Behavior %in% c("B1","B2"))
data.B12 <- prune_taxa(taxa_sums(data.B12)>0,data.B12)

data <- data.B12
#run deseq2
infdds = phyloseq_to_deseq2(data, ~Behavior)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Behavior <- relevel(infddsClean$Behavior, ref="B1")

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
write.csv(sigtab, file = "UNC_CD_colon_B2_vs_B1_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_colon_B2_vs_B1_BASELINE_25_all_taxa.csv")

#### CD Colon - B3 vs B1

data <- CD.colon.B.25

data.B13 <- subset_samples(data,sample_data(data)$Behavior %in% c("B1","B3"))
data.B13 <- prune_taxa(taxa_sums(data.B13)>0,data.B13)

data <- data.B13
#run deseq2
infdds = phyloseq_to_deseq2(data, ~Behavior)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Behavior <- relevel(infddsClean$Behavior, ref="B1")

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
write.csv(sigtab, file = "UNC_CD_colon_B3_vs_B1_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_colon_B3_vs_B1_BASELINE_25_all_taxa.csv")


#### CD Colon - B3 vs B2

data <- CD.colon.B.25

data.B23 <- subset_samples(data,sample_data(data)$Behavior %in% c("B2","B3"))
data.B23 <- prune_taxa(taxa_sums(data.B23)>0,data.B23)

data <- data.B23
#run deseq2
infdds = phyloseq_to_deseq2(data, ~Behavior)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Behavior <- relevel(infddsClean$Behavior, ref="B2")

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
write.csv(sigtab, file = "UNC_CD_colon_B3_vs_B2_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_colon_B3_vs_B2_BASELINE_25_all_taxa.csv")


#### COLON - nonIBD vs B1
data <- nonIBD.CD.NI.colon
sample_data(data)$Cat <- paste(sample_data(data)$DiseaseStatus, sample_data(data)$Behavior, sep= " ")

data.sub <- subset_samples(data,sample_data(data)$Cat %in% c("nonIBD ","CD B1"))
data.sub.25.taxa <- union(nonIBD.NI.colon.taxa.25,CD.colon.B.taxa.25)
data.sub <- prune_taxa(data.sub.25.taxa,data.sub)

data <- data.sub

#run deseq2
infdds = phyloseq_to_deseq2(data, ~Cat)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Cat <- relevel(infddsClean$Cat, ref="CD B1")

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
write.csv(sigtab, file = "UNC_colon_nonIBD_vs_CD_B1_BASELINE_25.csv")
write.csv(res, file = "UNC_colon_nonIBD_vs_CD_B1_BASELINE_25_all_taxa.csv")


#### COLON - nonIBD vs B2

data <- nonIBD.CD.NI.colon
sample_data(data)$Cat <- paste(sample_data(data)$DiseaseStatus, sample_data(data)$Behavior, sep= " ")

data.sub <- subset_samples(data,sample_data(data)$Cat %in% c("nonIBD ","CD B2"))
data.sub.25.taxa <- union(nonIBD.NI.colon.taxa.25,CD.colon.B.taxa.25)
data.sub <- prune_taxa(data.sub.25.taxa,data.sub)
data.sub2 <- data.sub
data <- data.sub

#run deseq2
infdds = phyloseq_to_deseq2(data, ~Cat)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Cat <- relevel(infddsClean$Cat, ref="CD B2")

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
write.csv(sigtab, file = "UNC_colon_nonIBD_vs_CD_B2_BASELINE_25.csv")
write.csv(res, file = "UNC_colon_nonIBD_vs_CD_B2_BASELINE_25_all_taxa.csv")

#### COLON -nonIBD vs B3
data <- nonIBD.CD.NI.colon
sample_data(data)$Cat <- paste(sample_data(data)$DiseaseStatus, sample_data(data)$Behavior, sep= " ")

data.sub <- subset_samples(data,sample_data(data)$Cat %in% c("nonIBD ","CD B3"))
data.sub.25.taxa <- union(nonIBD.NI.colon.taxa.25,CD.colon.B.taxa.25)
data.sub <- prune_taxa(data.sub.25.taxa,data.sub)

data <- data.sub

#run deseq2
infdds = phyloseq_to_deseq2(data, ~Cat)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Cat <- relevel(infddsClean$Cat, ref="CD B3")

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
write.csv(sigtab, file = "UNC_colon_nonIBD_vs_CD_B3_BASELINE_25.csv")
write.csv(res, file = "UNC_colon_nonIBD_vs_CD_B3_BASELINE_25_all_taxa.csv")

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus","Species")]


theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Species))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$GS <- paste(sigtabgen$Genus,sigtabgen$Species,sep=" ")
ggplot(sigtabgen, aes(y=GS, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "black", size = 0.2) +
  geom_point(size=8) + 
  scale_color_manual(values = microbiota.palette) +
  ggtitle("UNC Colon: nonIBD vs CD B3 (p = 0.05)") +
  theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5)) +
  xlab(expression(paste(log[2],"Fold Change")))



#### Ileum disease Behavior####

CD.NI.ileum.B2 <- subset_samples(CD.NI.ileum,sample_data(CD.NI.ileum)$Behavior == "B2")
CD.NI.ileum.B2 <- prune_taxa(taxa_sums(CD.NI.ileum.B2)>0, CD.NI.ileum.B2)

CD.NI.ileum.B3 <- subset_samples(CD.NI.ileum,sample_data(CD.NI.ileum)$Behavior == "B3")
CD.NI.ileum.B3 <- prune_taxa(taxa_sums(CD.NI.ileum.B3)>0, CD.NI.ileum.B3)

CD.NI.ileum.B2.taxa.25 <- getThresholdTaxa(CD.NI.ileum.B2,0.25)
CD.NI.ileum.B3.taxa.25 <- getThresholdTaxa(CD.NI.ileum.B3,0.25)

CD.ileum.B.taxa.25 <- union(CD.NI.ileum.B3.taxa.25,CD.NI.ileum.B2.taxa.25)

CD.ileum.B.25 <- prune_taxa(CD.ileum.B.taxa.25,CD.NI.ileum)

#### Ileum - nonIBD vs B2 ####
data <- nonIBD.CD.NI.ileum
sample_data(data)$Cat <- paste(sample_data(data)$DiseaseStatus, sample_data(data)$Behavior, sep= " ")

data.sub <- subset_samples(data,sample_data(data)$Cat %in% c("nonIBD ","CD B2"))
data.sub.25.taxa <- union(nonIBD.NI.ileum.taxa.25,CD.ileum.B.taxa.25)
data.sub <- prune_taxa(data.sub.25.taxa,data.sub)

data <- data.sub

#run deseq2
infdds = phyloseq_to_deseq2(data, ~Cat)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Cat <- relevel(infddsClean$Cat, ref="CD B2")

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
write.csv(sigtab, file = "UNC_ileum_nonIBD_vs_CD_B2_BASELINE_25.csv")
write.csv(res, file = "UNC_ileum_nonIBD_vs_CD_B2_BASELINE_25_all_taxa.csv")

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus","Species")]


theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Species))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$GS <- paste(sigtabgen$Genus,sigtabgen$Species,sep=" ")
ggplot(sigtabgen, aes(y=GS, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "black", size = 0.2) +
  geom_point(size=8) + 
  scale_color_manual(values = microbiota.palette) +
  ggtitle("UNC Ileum: nonIBD vs CD B2 (p = 0.05)") +
  theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5)) +
  xlab(expression(paste(log[2],"Fold Change")))

#### Ileum - nonIBD vs B3
data <- nonIBD.CD.NI.ileum
sample_data(data)$Cat <- paste(sample_data(data)$DiseaseStatus, sample_data(data)$Behavior, sep= " ")

data.sub <- subset_samples(data,sample_data(data)$Cat %in% c("nonIBD ","CD B3"))
data.sub.25.taxa <- union(nonIBD.NI.ileum.taxa.25,CD.ileum.B.taxa.25)
data.sub <- prune_taxa(data.sub.25.taxa,data.sub)

data <- data.sub

#run deseq2
infdds = phyloseq_to_deseq2(data, ~Cat)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Cat <- relevel(infddsClean$Cat, ref="CD B3")

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
write.csv(sigtab, file = "UNC_ileum_nonIBD_vs_CD_B3_BASELINE_25.csv")
write.csv(res, file = "UNC_ileum_nonIBD_vs_CD_B3_BASELINE_25_all_taxa.csv")

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus","Species")]


theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Species))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$GS <- paste(sigtabgen$Genus,sigtabgen$Species,sep=" ")
ggplot(sigtabgen, aes(y=GS, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "black", size = 0.2) +
  geom_point(size=8) + 
  scale_color_manual(values = microbiota.palette) +
  ggtitle("UNC Ileum: nonIBD vs CD B3 (p = 0.05)") +
  theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5)) +
  xlab(expression(paste(log[2],"Fold Change")))

#### CD Ileum - B3 vs B2
data <- CD.ileum.B.25

#run deseq2
infdds = phyloseq_to_deseq2(data, ~Behavior)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Behavior <- relevel(infddsClean$Behavior, ref="B2")

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
write.csv(sigtab, file = "UNC_CD_ileum_B3_vs_CD_B2_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_ileum_B3_vs_CD_B2_BASELINE_25_all_taxa.csv")

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus","Species")]


theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Species))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$GS <- paste(sigtabgen$Genus,sigtabgen$Species,sep=" ")
ggplot(sigtabgen, aes(y=GS, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "black", size = 0.2) +
  geom_point(size=8) + 
  scale_color_manual(values = microbiota.palette) +
  ggtitle("UNC CD Ileum: B3 vs B2 (p = 0.05)") +
  theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5)) +
  xlab(expression(paste(log[2],"Fold Change")))


#### LOCATION ####
#### Colon - Disease Location ####
CD.NI.colon.L2 <- subset_samples(CD.NI.colon,sample_data(CD.NI.colon)$Location == "L2")
CD.NI.colon.L2 <- prune_taxa(taxa_sums(CD.NI.colon.L2)>0, CD.NI.colon.L2)

CD.NI.colon.L3 <- subset_samples(CD.NI.colon,sample_data(CD.NI.colon)$Location == "L3")
CD.NI.colon.L3 <- prune_taxa(taxa_sums(CD.NI.colon.L3)>0, CD.NI.colon.L3)

CD.NI.colon.L2.taxa.25 <- getThresholdTaxa(CD.NI.colon.L2,0.25)
CD.NI.colon.L3.taxa.25 <- getThresholdTaxa(CD.NI.colon.L3,0.25)

CD.colon.L.taxa.25 <- union(CD.NI.colon.L3.taxa.25,CD.NI.colon.L2.taxa.25)

CD.colon.L.25 <- prune_taxa(CD.colon.L.taxa.25,CD.NI.colon)

CD.NI.ileum.L3 <- subset_samples(CD.NI.ileum,sample_data(CD.NI.ileum)$Location == "L3")
CD.NI.ileum.L3 <- prune_taxa(taxa_sums(CD.NI.ileum.L3)>0, CD.NI.ileum.L3)

CD.NI.ileum.L3.taxa.25 <- getThresholdTaxa(CD.NI.ileum.L3,0.25)

#### CD Colon - L3 vs L2

data <- CD.colon.L.25

data.L23 <- subset_samples(data,sample_data(data)$Location %in% c("L2","L3"))
data.L23 <- prune_taxa(taxa_sums(data.L23)>0,data.L23)

data <- data.L23
#run deseq2
infdds = phyloseq_to_deseq2(data, ~Location)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Location <- relevel(infddsClean$Location, ref="L2")

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
write.csv(sigtab, file = "UNC_CD_colon_L3_vs_L2_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_colon_L3_vs_L2_BASELINE_25_all_taxa.csv")


#### Colon - nonIBD vs CD L2

data <- nonIBD.CD.NI.colon.25

sample_data(data)$Cat <- paste(sample_data(data)$DiseaseStatus, sample_data(data)$Location, sep= " ")

data.sub <- subset_samples(data,sample_data(data)$Cat %in% c("nonIBD ","CD L2"))
data.sub.25.taxa <- union(nonIBD.NI.colon.taxa.25,CD.colon.L.taxa.25)
data.sub <- prune_taxa(data.sub.25.taxa,data.sub)

data <- data.sub

#run deseq2
infdds = phyloseq_to_deseq2(data, ~Cat)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Cat <- relevel(infddsClean$Cat, ref="CD L2")

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
write.csv(sigtab, file = "UNC_colon_nonIBD_vs_CD_L2_BASELINE_25.csv")
write.csv(res, file = "UNC_colon_nonIBD_vs_CD_L2_BASELINE_25_all_taxa.csv")


#### Colon - nonIBD vs CD L3

data <- nonIBD.CD.NI.colon.25

sample_data(data)$Cat <- paste(sample_data(data)$DiseaseStatus, sample_data(data)$Location, sep= " ")

data.sub <- subset_samples(data,sample_data(data)$Cat %in% c("nonIBD ","CD L3"))
data.sub.25.taxa <- union(nonIBD.NI.colon.taxa.25,CD.colon.L.taxa.25)
data.sub <- prune_taxa(data.sub.25.taxa,data.sub)

data <- data.sub

#run deseq2
infdds = phyloseq_to_deseq2(data, ~Cat)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Cat <- relevel(infddsClean$Cat, ref="CD L3")

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
write.csv(sigtab, file = "UNC_colon_nonIBD_vs_CD_L3_BASELINE_25.csv")
write.csv(res, file = "UNC_colon_nonIBD_vs_CD_L3_BASELINE_25_all.csv")


#### CD L3 - Tissue

data <- CD.NI.patients
data.L3 <- subset_samples(data,sample_data(data)$Location == "L3")

data.L3.taxa.25 <- union(CD.NI.colon.L3.taxa.25,CD.NI.ileum.L3.taxa.25)

data.L3.25 <- prune_taxa(data.L3.taxa.25,data.L3)

data <- data.L3.25

#run deseq2
infdds = phyloseq_to_deseq2(data, ~AnatomSite)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$AnatomSite <- relevel(infddsClean$AnatomSite, ref="colon")

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
write.csv(sigtab, file = "UNC_CD_L3_ileum_vs_colon_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_L3_ileum_vs_colon_BASELINE_25_all_taxa.csv")

#### CD L3 (no B1) - Tissue

data <- CD.NI.patients
data.L3 <- subset_samples(data,sample_data(data)$Location == "L3")

data.L3.taxa.25 <- union(CD.NI.colon.L3.taxa.25,CD.NI.ileum.L3.taxa.25)

data.L3.25 <- prune_taxa(data.L3.taxa.25,data.L3)

data <- data.L3.25

data.L3.n.B1 <- subset_samples(data,sample_data(data)$Behavior != "B1")
data.L3.n.B1 <- prune_taxa(taxa_sums(data.L3.n.B1)>0,data.L3.n.B1)
data <- data.L3.n.B1
#run deseq2
infdds = phyloseq_to_deseq2(data, ~AnatomSite)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$AnatomSite <- relevel(infddsClean$AnatomSite, ref="colon")

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
write.csv(sigtab, file = "UNC_CD_L3_no_B1_ileum_vs_colon_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_L3_no_B1_ileum_vs_colon_BASELINE_25_all_taxa.csv")



#### MEDICATION - COLON ####
CD.NI.colon.TNF.y <- subset_samples(CD.NI.colon,sample_data(CD.NI.colon)$TNF == "1")
CD.NI.colon.TNF.y <- prune_taxa(taxa_sums(CD.NI.colon.TNF.y)>0, CD.NI.colon.TNF.y)

CD.NI.colon.STR.y <- subset_samples(CD.NI.colon,sample_data(CD.NI.colon)$Steroids == "1")
CD.NI.colon.STR.y <- prune_taxa(taxa_sums(CD.NI.colon.STR.y)>0, CD.NI.colon.STR.y)

CD.NI.colon.IMM.y <- subset_samples(CD.NI.colon,sample_data(CD.NI.colon)$Immunomodulators == "1")
CD.NI.colon.IMM.y <- prune_taxa(taxa_sums(CD.NI.colon.IMM.y)>0, CD.NI.colon.IMM.y)

CD.NI.colon.TNF.y.taxa.25 <- getThresholdTaxa(CD.NI.colon.TNF.y,0.25)
CD.NI.colon.STR.y.taxa.25 <- getThresholdTaxa(CD.NI.colon.STR.y,0.25)
CD.NI.colon.IMM.y.taxa.25 <- getThresholdTaxa(CD.NI.colon.IMM.y,0.25)

CD.NI.colon.TNF.n <- subset_samples(CD.NI.colon,sample_data(CD.NI.colon)$TNF == "0")
CD.NI.colon.TNF.n <- prune_taxa(taxa_sums(CD.NI.colon.TNF.n)>0, CD.NI.colon.TNF.n)

CD.NI.colon.STR.n <- subset_samples(CD.NI.colon,sample_data(CD.NI.colon)$Steroids == "0")
CD.NI.colon.STR.n <- prune_taxa(taxa_sums(CD.NI.colon.STR.n)>0, CD.NI.colon.STR.n)

CD.NI.colon.IMM.n <- subset_samples(CD.NI.colon,sample_data(CD.NI.colon)$Immunomodulators == "0")
CD.NI.colon.IMM.n <- prune_taxa(taxa_sums(CD.NI.colon.IMM.n)>0, CD.NI.colon.IMM.n)

CD.NI.colon.TNF.n.taxa.25 <- getThresholdTaxa(CD.NI.colon.TNF.n,0.25)
CD.NI.colon.STR.n.taxa.25 <- getThresholdTaxa(CD.NI.colon.STR.n,0.25)
CD.NI.colon.IMM.n.taxa.25 <- getThresholdTaxa(CD.NI.colon.IMM.n,0.25)

CD.colon.TNF.taxa.25 <- union(CD.NI.colon.TNF.y.taxa.25,CD.NI.colon.TNF.n.taxa.25)
CD.colon.STR.taxa.25 <- union(CD.NI.colon.STR.y.taxa.25,CD.NI.colon.STR.n.taxa.25)
CD.colon.IMM.taxa.25 <- union(CD.NI.colon.IMM.y.taxa.25,CD.NI.colon.IMM.n.taxa.25)

CD.colon.TNF.25 <- prune_taxa(CD.colon.TNF.taxa.25,CD.NI.colon)
CD.colon.STR.25 <- prune_taxa(CD.colon.STR.taxa.25,CD.NI.colon)
CD.colon.IMM.25 <- prune_taxa(CD.colon.IMM.taxa.25,CD.NI.colon)
#### CD Colon - TNF

data <- CD.colon.TNF.25
data.TNF <- subset_samples(CD.colon.TNF.25,sample_data(CD.colon.TNF.25)$TNF %in% c("0","1"))
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
write.csv(sigtab, file = "UNC_CD_colon_TNF_usage_vs_no_usage_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_colon_TNF_usage_vs_no_usage_BASELINE_25_all_taxa.csv")

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus","Species")]


theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Species))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$GS <- paste(sigtabgen$Genus,sigtabgen$Species,sep=" ")
ggplot(sigtabgen, aes(y=GS, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "black", size = 0.2) +
  geom_point(size=8) + 
  scale_color_manual(values = microbiota.palette) +
  ggtitle("UNC CD Colon: anti-TNF Usage vs No anti-TNF Usage (p = 0.05)") +
  theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5)) +
  xlab(expression(paste(log[2],"Fold Change")))

#### CD Colon - Steroids

data <- CD.colon.STR.25
data.str <- subset_samples(CD.colon.STR.25,sample_data(CD.colon.STR.25)$Steroids %in% c("0","1"))
data <- data.str

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
write.csv(sigtab, file = "UNC_CD_colon_steroid_usage_vs_no_usage_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_colon_STR_usage_vs_no_usage_BASELINE_25_all_taxa.csv")

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus","Species")]


theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Species))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$GS <- paste(sigtabgen$Genus,sigtabgen$Species,sep=" ")
ggplot(sigtabgen, aes(y=GS, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "black", size = 0.2) +
  geom_point(size=8) + 
  scale_color_manual(values = microbiota.palette) +
  ggtitle("UNC CD Colon: Steroid Usage vs No Steroid Usage (p = 0.05)") +
  theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5)) +
  xlab(expression(paste(log[2],"Fold Change")))

#### CD Colon - Immunomodulators

data <- CD.colon.IMM.25
data.imm <- subset_samples(CD.colon.IMM.25,sample_data(CD.colon.IMM.25)$Immunomodulators %in% c("0","1"))
data <- data.imm

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
write.csv(sigtab, file = "UNC_CD_colon_Immunomodulators_usage_vs_no_usage_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_colon_IMM_usage_vs_no_usage_BASELINE_25_all_taxa.csv")

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus","Species")]


theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Species))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$GS <- paste(sigtabgen$Genus,sigtabgen$Species,sep=" ")
ggplot(sigtabgen, aes(y=GS, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "black", size = 0.2) +
  geom_point(size=8) + 
  scale_color_manual(values = microbiota.palette) +
  ggtitle("UNC CD Colon: Immunomodulators Usage vs No Immunomodulators Usage (p = 0.05)") +
  theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5)) +
  xlab(expression(paste(log[2],"Fold Change")))



#### MEDICATION - ILEUM ####
CD.NI.ileum.TNF.y <- subset_samples(CD.NI.ileum,sample_data(CD.NI.ileum)$TNF == "1")
CD.NI.ileum.TNF.y <- prune_taxa(taxa_sums(CD.NI.ileum.TNF.y)>0, CD.NI.ileum.TNF.y)

CD.NI.ileum.STR.y <- subset_samples(CD.NI.ileum,sample_data(CD.NI.ileum)$Steroids == "1")
CD.NI.ileum.STR.y <- prune_taxa(taxa_sums(CD.NI.ileum.STR.y)>0, CD.NI.ileum.STR.y)

CD.NI.ileum.IMM.y <- subset_samples(CD.NI.ileum,sample_data(CD.NI.ileum)$Immunomodulators == "1")
CD.NI.ileum.IMM.y <- prune_taxa(taxa_sums(CD.NI.ileum.IMM.y)>0, CD.NI.ileum.IMM.y)

CD.NI.ileum.ASA.y <- subset_samples(CD.NI.ileum,sample_data(CD.NI.ileum)$X5_ASA == "1")
CD.NI.ileum.ASA.y <- prune_taxa(taxa_sums(CD.NI.ileum.ASA.y)>0, CD.NI.ileum.ASA.y)

CD.NI.ileum.TNF.y.taxa.25 <- getThresholdTaxa(CD.NI.ileum.TNF.y,0.25)
CD.NI.ileum.STR.y.taxa.25 <- getThresholdTaxa(CD.NI.ileum.STR.y,0.25)
CD.NI.ileum.IMM.y.taxa.25 <- getThresholdTaxa(CD.NI.ileum.IMM.y,0.25)
CD.NI.ileum.ASA.y.taxa.25 <- getThresholdTaxa(CD.NI.ileum.ASA.y,0.25)

CD.NI.ileum.TNF.n <- subset_samples(CD.NI.ileum,sample_data(CD.NI.ileum)$TNF == "0")
CD.NI.ileum.TNF.n <- prune_taxa(taxa_sums(CD.NI.ileum.TNF.n)>0, CD.NI.ileum.TNF.n)

CD.NI.ileum.STR.n <- subset_samples(CD.NI.ileum,sample_data(CD.NI.ileum)$Steroids == "0")
CD.NI.ileum.STR.n <- prune_taxa(taxa_sums(CD.NI.ileum.STR.n)>0, CD.NI.ileum.STR.n)

CD.NI.ileum.IMM.n <- subset_samples(CD.NI.ileum,sample_data(CD.NI.ileum)$Immunomodulators == "0")
CD.NI.ileum.IMM.n <- prune_taxa(taxa_sums(CD.NI.ileum.IMM.n)>0, CD.NI.ileum.IMM.n)

CD.NI.ileum.ASA.n <- subset_samples(CD.NI.ileum,sample_data(CD.NI.ileum)$X5_ASA == "0")
CD.NI.ileum.ASA.n <- prune_taxa(taxa_sums(CD.NI.ileum.ASA.n)>0, CD.NI.ileum.ASA.n)

CD.NI.ileum.TNF.n.taxa.25 <- getThresholdTaxa(CD.NI.ileum.TNF.n,0.25)
CD.NI.ileum.STR.n.taxa.25 <- getThresholdTaxa(CD.NI.ileum.STR.n,0.25)
CD.NI.ileum.IMM.n.taxa.25 <- getThresholdTaxa(CD.NI.ileum.IMM.n,0.25)
CD.NI.ileum.ASA.n.taxa.25 <- getThresholdTaxa(CD.NI.ileum.ASA.n,0.25)

CD.ileum.TNF.taxa.25 <- union(CD.NI.ileum.TNF.y.taxa.25,CD.NI.ileum.TNF.n.taxa.25)
CD.ileum.STR.taxa.25 <- union(CD.NI.ileum.STR.y.taxa.25,CD.NI.ileum.STR.n.taxa.25)
CD.ileum.IMM.taxa.25 <- union(CD.NI.ileum.IMM.y.taxa.25,CD.NI.ileum.IMM.n.taxa.25)
CD.ileum.ASA.taxa.25 <- union(CD.NI.ileum.ASA.y.taxa.25,CD.NI.ileum.ASA.n.taxa.25)

CD.ileum.TNF.25 <- prune_taxa(CD.ileum.TNF.taxa.25,CD.NI.ileum)
CD.ileum.STR.25 <- prune_taxa(CD.ileum.STR.taxa.25,CD.NI.ileum)
CD.ileum.IMM.25 <- prune_taxa(CD.ileum.IMM.taxa.25,CD.NI.ileum)
CD.ileum.ASA.25 <- prune_taxa(CD.ileum.ASA.taxa.25,CD.NI.ileum)


#### CD Ileum - TNF

data <- CD.ileum.TNF.25
data.TNF <- subset_samples(CD.ileum.TNF.25,sample_data(CD.ileum.TNF.25)$TNF %in% c("0","1"))
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
write.csv(sigtab, file = "UNC_CD_ileum_TNF_usage_vs_no_usage_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_ileum_TNF_usage_vs_no_usage_BASELINE_25_all_taxa.csv")

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus","Species")]


theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Species))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$GS <- paste(sigtabgen$Genus,sigtabgen$Species,sep=" ")
ggplot(sigtabgen, aes(y=GS, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "black", size = 0.2) +
  geom_point(size=8) + 
  scale_color_manual(values = microbiota.palette) +
  ggtitle("UNC CD Ileum: anti-TNF Usage vs No anti-TNF Usage (p = 0.05)") +
  theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5)) +
  xlab(expression(paste(log[2],"Fold Change")))

#### CD Ileum - Steroids

data <- CD.ileum.STR.25
data.str <- subset_samples(CD.ileum.STR.25,sample_data(CD.ileum.STR.25)$Steroids %in% c("0","1"))
data <- data.str

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
write.csv(sigtab, file = "UNC_CD_ileum_steroid_usage_vs_no_usage_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_ileum_STR_usage_vs_no_usage_BASELINE_25_all_taxa.csv")

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus","Species")]


theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Species))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$GS <- paste(sigtabgen$Genus,sigtabgen$Species,sep=" ")
ggplot(sigtabgen, aes(y=GS, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "black", size = 0.2) +
  geom_point(size=8) + 
  scale_color_manual(values = microbiota.palette) +
  ggtitle("UNC CD Ileum: Steroid Usage vs No Steroid Usage (p = 0.05)") +
  theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5)) +
  xlab(expression(paste(log[2],"Fold Change")))

#### CD Ileum - Immunomodulators

data <- CD.ileum.IMM.25
data.imm <- subset_samples(CD.ileum.IMM.25,sample_data(CD.ileum.IMM.25)$Immunomodulators %in% c("0","1"))
data <- data.imm

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
write.csv(sigtab, file = "UNC_CD_ileum_Immunomodulators_usage_vs_no_usage_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_ileum_IMM_usage_vs_no_usage_BASELINE_25_all_taxa.csv")

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus","Species")]


theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Species))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$GS <- paste(sigtabgen$Genus,sigtabgen$Species,sep=" ")
ggplot(sigtabgen, aes(y=GS, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "black", size = 0.2) +
  geom_point(size=8) + 
  scale_color_manual(values = microbiota.palette) +
  ggtitle("UNC CD Ileum: Immunomodulators Usage vs No Immunomodulators Usage (p = 0.05)") +
  theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5)) +
  xlab(expression(paste(log[2],"Fold Change")))

#### CD Ileum - 5-ASA

data <- CD.ileum.ASA.25
data.asa <- subset_samples(CD.ileum.ASA.25,sample_data(CD.ileum.ASA.25)$X5_ASA %in% c("0","1"))
data <- data.asa

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
write.csv(sigtab, file = "UNC_CD_ileum_5ASA_usage_vs_no_usage_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_ileum_5ASA_usage_vs_no_usage_BASELINE_25_all_taxa.csv")

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus","Species")]


theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Species))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))
sigtabgen$GS <- paste(sigtabgen$Genus,sigtabgen$Species,sep=" ")
ggplot(sigtabgen, aes(y=GS, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "black", size = 0.2) +
  geom_point(size=8) + 
  scale_color_manual(values = microbiota.palette) +
  ggtitle("UNC CD Ileum: 5-ASA Usage vs No 5-ASA Usage (p = 0.05)") +
  theme_classic() +
  theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5)) +
  xlab(expression(paste(log[2],"Fold Change")))

#### UNC - Recurrence - ILEUM ####
CD.NI.ileum.rec.y <- subset_samples(CD.NI.ileum,sample_data(CD.NI.ileum)$Recurrence == "recurrence")
CD.NI.ileum.rec.y <- prune_taxa(taxa_sums(CD.NI.ileum.rec.y)>0, CD.NI.ileum.rec.y)

CD.NI.ileum.rec.n <- subset_samples(CD.NI.ileum,sample_data(CD.NI.ileum)$Recurrence == "no recurrence")
CD.NI.ileum.rec.n <- prune_taxa(taxa_sums(CD.NI.ileum.rec.n)>0, CD.NI.ileum.rec.n)

CD.NI.ileum.rec.y.taxa.25 <- getThresholdTaxa(CD.NI.ileum.rec.y,0.25)
CD.NI.ileum.rec.n.taxa.25 <- getThresholdTaxa(CD.NI.ileum.rec.n,0.25)

CD.ileum.rec.taxa.25 <- union(CD.NI.ileum.rec.y.taxa.25,CD.NI.ileum.rec.n.taxa.25)

CD.ileum.rec.25 <- prune_taxa(CD.ileum.rec.taxa.25,CD.NI.ileum)


data <- CD.ileum.rec.25
data.rec <- subset_samples(data,sample_data(data)$Recurrence %in% c("no recurrence","recurrence"))
data <- data.rec

#run deseq2
infdds = phyloseq_to_deseq2(data, ~Recurrence)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Recurrence <- relevel(infddsClean$Recurrence, ref="no recurrence")

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
write.csv(sigtab, file = "UNC_CD_ileum_rec_vs_no_rec_BASELINE_25.csv")
write.csv(res, file = "UNC_CD_ileum_rec_vs_no_rec_BASELINE_25_all_taxa.csv")

#### WashU - Recurrence - ILEUM ####
data.r
WashU.ileum.rec.y <- subset_samples(data.r,sample_data(data.r)$Recurrence == "recurrence")
WashU.ileum.rec.y <- prune_taxa(taxa_sums(WashU.ileum.rec.y)>0, WashU.ileum.rec.y)

WashU.ileum.rec.n <- subset_samples(data.r,sample_data(data.r)$Recurrence == "no recurrence")
WashU.ileum.rec.n <- prune_taxa(taxa_sums(WashU.ileum.rec.n)>0, WashU.ileum.rec.n)

WashU.ileum.rec.y.taxa.25 <- getThresholdTaxa(WashU.ileum.rec.y,0.25)
WashU.ileum.rec.n.taxa.25 <- getThresholdTaxa(WashU.ileum.rec.n,0.25)

WashU.ileum.rec.taxa.25 <- union(WashU.ileum.rec.y.taxa.25,WashU.ileum.rec.n.taxa.25)

WashU.ileum.rec.25 <- prune_taxa(WashU.ileum.rec.taxa.25,data.r)


data <- WashU.ileum.rec.25
data.rec <- subset_samples(data,sample_data(data)$Recurrence %in% c("no recurrence","recurrence"))
data <- data.rec

#run deseq2
infdds = phyloseq_to_deseq2(data, ~Recurrence)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Recurrence <- relevel(infddsClean$Recurrence, ref="no recurrence")

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
write.csv(sigtab, file = "WashU_CD_ileum_rec_vs_no_rec_BASELINE_25.csv")
write.csv(res, file = "WashU_CD_ileum_rec_vs_no_rec_BASELINE_25_all_taxa.csv")


#### WashU - Recurrence - ILEUM ####
data.r <- ps.WashU.patients

WashU.B2 <- subset_samples(data.r,sample_data(data.r)$Behavior == "B2")
WashU.B2 <- prune_taxa(taxa_sums(WashU.B2)>0, WashU.B2)

WashU.B3 <- subset_samples(data.r,sample_data(data.r)$Behavior == "B3")
WashU.B3 <- prune_taxa(taxa_sums(WashU.B3)>0, WashU.B3)

WashU.ileum.B2.taxa.25 <- getThresholdTaxa(WashU.B2,0.25)
WashU.ileum.B3.taxa.25 <- getThresholdTaxa(WashU.B3,0.25)

WashU.ileum.B23.taxa.25 <- union(WashU.ileum.B2.taxa.25,WashU.ileum.B3.taxa.255)

WashU.ileum.B23.25 <- prune_taxa(WashU.ileum.B23.taxa.25,data.r)


data <- WashU.ileum.B23.25
data.b23.w <- subset_samples(data,sample_data(data)$Behavior %in% c("B2","B3"))
data <- data.b23.w

#run deseq2
infdds = phyloseq_to_deseq2(data, ~Behavior)
geoMeans = apply(DESeq2::counts(infdds), 1, gm_mean)
infdds = DESeq2::estimateSizeFactors(infdds, geoMeans = geoMeans)
infdds = DESeq2::DESeq(infdds, fitType="local")

infddsClean <- infdds[which(mcols(infdds)$betaConv),]
infddsClean$Behavior <- relevel(infddsClean$Behavior, ref="B2")

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
write.csv(sigtab, file = "WashU_CD_ileum_B3_vs_B2_BASELINE_25.csv")
write.csv(res, file = "WashU_CD_ileum_B3_vs_B2_BASELINE_25_all_taxa.csv")

#### LEFSE FORMATTING ####
colon_disease <- MaaslinFormating2(nonIBD.CD.NI.colon.25)
write.csv(colon_disease, "colon_disease_maaslin.csv")

ileum_disease <- MaaslinFormating2(nonIBD.CD.NI.ileum.25)
write.csv(ileum_disease, "ileum_disease_maaslin.csv")

nonIBD_tissue <- MaaslinFormating2(nonIBD.NI.tissue.25)
write.csv(nonIBD_tissue,"nonIBD_tissue_maaslin.csv")

CD_tissue <- MaaslinFormating2(CD.NI.tissue.25)
write.csv(CD_tissue,"CD_tissue_maaslin.csv")

CD_B23_tissue <- MaaslinFormating2(data.b23)
write.csv(CD_B23_tissue,"CD_B2B3_tissue_maaslin.csv")

CD_ileum_inf <- MaaslinFormating2(CD.ileum.m.25)
write.csv(CD_ileum_inf,"CD_ileum_inflammation_maaslin.csv")

#### TO DO ###
CD_colon_B1B2 <- MaaslinFormating2(data.B12)
CD_colon_B1B3 <- MaaslinFormating2(data.B13)
CD_colon_B2B3 <- MaaslinFormating2(data.B23)

write.csv(CD_colon_B1B2,"CD_colon_B12_maaslin.csv")
write.csv(CD_colon_B1B3,"CD_colon_B13_maaslin.csv")
write.csv(CD_colon_B2B3,"CD_colon_B23_maaslin.csv")

colon_nonIBD_CDB1 <- MaaslinFormating2(data.sub)
colon_nonIBD_CDB2 <- MaaslinFormating2(data.sub2)
colon_nonIBD_CDB3 <- MaaslinFormating2(data.sub3)

write.csv(colon_nonIBD_CDB1,"colon_nonIBD_CDB1_maaslin.csv")
write.csv(colon_nonIBD_CDB2,"colon_nonIBD_CDB2_maaslin.csv")
write.csv(colon_nonIBD_CDB3,"colon_nonIBD_CDB3_maaslin.csv")


CD_ileum_B2B3 <- MaaslinFormating2(CD.ileum.B.25)
ileum_nonIBD_CDB2 <- MaaslinFormating2(data.sub2)
ileum_nonIBD_CDB3 <- MaaslinFormating2(data.sub3)

write.csv(CD_ileum_B2B3,"CD_ileum_B23_maaslin.csv")
write.csv(ileum_nonIBD_CDB2,"ileum_nonIBD_CDB2_maaslin.csv")
write.csv(ileum_nonIBD_CDB3,"ileum_nonIBD_CDB3_maaslin.csv")

CD_B2_tissue <- MaaslinFormating2(data.b2)
CD_B2_concordant <- MaaslinFormating2(data.concordant)
CD_B3_tissue <- MaaslinFormating2(data.b3)

write.csv(CD_B2_tissue, "CD_B2_tissue_maaslin.csv")
write.csv(CD_B2_concordant, "CD_B2_concordant_maaslin.csv")
write.csv(CD_B3_tissue,"CD_B3_tissue_maaslin.csv")

CD_colon_L2L3 <- MaaslinFormating2(data.L23)
colon_nonIBD_CDL2 <- MaaslinFormating2(data.sub2)
colon_nonIBD_CDL3 <- MaaslinFormating2(data.su3)

write.csv(CD_colon_L2L3,"CD_colon_L23_maaslin.csv")
write.csv(colon_nonIBD_CDL2,"colon_nonIBD_CDL2_maaslin.csv")
write.csv(colon_nonIBD_CDL3,"colon_nonIBD_CDL3_maaslin.csv")

CD_L3_tissue <- MaaslinFormating2(data.L3.25)
write.csv(CD_L3_tissue,"CD_L3_tissue_maaslin.csv")

CD_colon_TNF_usage <- MaaslinFormating2(data.TNF)
CD_colon_imm_usage <- MaaslinFormating2(data.imm)
CD_colon_str_usage <- MaaslinFormating2(data.str)
CD_colon_asa_usage <- MaaslinFormating2(data.asa)

write.csv(CD_colon_TNF_usage,"CD_colon_TNF_usage.csv")
write.csv(CD_colon_imm_usage,"CD_colon_IMM_usage.csv")
write.csv(CD_colon_str_usage,"CD_colon_STR_usage.csv")
write.csv(CD_colon_asa_usage,"CD_colon_5ASA_usage.csv")

CD_ileum_TNF_usage <- MaaslinFormating2(data.TNF)
CD_ileum_imm_usage <- MaaslinFormating2(data.imm)
CD_ileum_str_usage <- MaaslinFormating2(data.str)
CD_ileum_asa_usage <- MaaslinFormating2(data.asa)

write.csv(CD_ileum_TNF_usage,"CD_ileum_TNF_usage_maaslin.csv")
write.csv(CD_ileum_imm_usage,"CD_ileum_IMM_usage_maaslin.csv")
write.csv(CD_ileum_str_usage,"CD_ileum_STR_usage_maaslin.csv")
write.csv(CD_ileum_asa_usage,"CD_ileum_5ASA_usage_maaslin.csv")

ileum_recur <- MaaslinFormating2(data.rec)
write.csv(ileum_recur,"ileal_recurrence.csv")

washu_ileum_recur <- MaaslinFormating2(data.rec)
write.csv(washu_ileum_recur,"ileal_recurrence_WASHU.csv")
