library(tidyr)
library(dplyr)
library(tibble)
#palettes
microbiota.palette <- c("Acidobacteria"="#be5400","Actinobacteria"="#d92728","Armatimonadetes"="#d36ff8","Bacteroidetes"="#0955da","Candidatus_Saccharibacteria"="#88374d","Chloroflexi"="#99c618","Cyanobacteria/Chloroplast"="#6adab2","Deferribacteres"= "#93aa00","Deinococcus-Thermus"="#00984c","Firmicutes"="#3b8b00","Fusobacteria"="#ff5fa3","Gemmatimonadetes"="#63b6ff","Lentisphaerae"="#ff8f54","Parcubacteria"="#006045","Proteobacteria"="#72349e","Spirochaetes"="#fbaee3","SR1"="#50570e","Synergistetes"="#fabb5e","Tenericutes"="#7b4601","Verrucomicrobia"="#01cacc")

#function for DESeq2
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
#### TEST FUNCTION ####
getThresholdGenus = function(physeq,thresholdValue){
  prevdf = apply(X = otu_table(physeq), 
                 MARGIN = ifelse(taxa_are_rows(physeq),yes=1, no = 2),
                 FUN = function(x){sum(x>0)})
  prevdf = data.frame(Prevalence = prevdf,TotalAbundance = taxa_sums(physeq),tax_table(physeq))
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq,"Phylum"))
  
  prevalenceThreshold = thresholdValue * nsamples(physeq)
  keepTaxaGenus <- prevdf1$Genus[(prevdf1$Prevalence >= prevalenceThreshold)]
  unique.Genera <- unique(keepTaxaGenus)
  return(unique.Genera)
}

getThresholdPhyseq = function(physeq,thresholdValue){
  prevdf = apply(X = otu_table(physeq), 
                 MARGIN = ifelse(taxa_are_rows(physeq),yes=1, no = 2),
                 FUN = function(x){sum(x>0)})
  prevdf = data.frame(Prevalence = prevdf,TotalAbundance = taxa_sums(physeq),tax_table(physeq))
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq,"Phylum"))
  
  prevalenceThreshold = thresholdValue * nsamples(physeq)
  keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
  ps2 = prune_taxa(keepTaxa, physeq)
  #unique.Taxa <- unique(keepTaxa)
  #return(unique.Taxa)
  return(ps2)
}

getThresholdPhyseq2 = function(physeq,thresholdValue){
  prevdf = apply(X = otu_table(physeq), 
                 MARGIN = ifelse(taxa_are_rows(physeq),yes=1, no = 2),
                 FUN = function(x){sum(x>0)})
  prevdf = data.frame(Prevalence = prevdf,TotalAbundance = taxa_sums(physeq),tax_table(physeq))
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq,"Phylum"))
  keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= thresholdValue)]
  ps2 = prune_taxa(keepTaxa, physeq)
  #unique.Taxa <- unique(keepTaxa)
  #return(unique.Taxa)
  return(ps2)
}

# getThresholdTaxa = function(physeq,thresholdValue){
#   prevdf = apply(X = otu_table(physeq), 
#                  MARGIN = ifelse(taxa_are_rows(physeq),yes=1, no = 2),
#                  FUN = function(x){sum(x>0)})
#   prevdf = data.frame(Prevalence = prevdf,TotalAbundance = taxa_sums(physeq),tax_table(physeq))
#   prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq,"Phylum"))
#   
#   prevalenceThreshold = thresholdValue * nsamples(physeq)
#   keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
#   unique.Taxa <- unique(keepTaxa)
#   return(unique.Taxa)
# }

getThresholdTaxa_old = function(physeq,thresholdValue){
  prevdf = apply(X = otu_table(physeq), 
                 MARGIN = ifelse(taxa_are_rows(physeq),yes=1, no = 2),
                 FUN = function(x){sum(x>0)})
  prevdf = data.frame(Prevalence = prevdf,TotalAbundance = taxa_sums(physeq),tax_table(physeq))
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq,"Phylum"))
  
  prevalenceThreshold = thresholdValue * nsamples(physeq)
  if (prevalenceThreshold > 1)
    keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence) >= prevalenceThreshold]
  else
    keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence) >= 2]
  unique.Taxa <- unique(keepTaxa)
  return(unique.Taxa)
}

#THIS IS THE CORRECT VERSION  -2018
getThresholdTaxa = function(physeq,thresholdValue){
  prevdf = apply(X = otu_table(physeq), 
                 MARGIN = ifelse(taxa_are_rows(physeq),yes=1, no = 2),
                 FUN = function(x){sum(x>0)})
  prevdf = data.frame(Prevalence = prevdf,TotalAbundance = taxa_sums(physeq),tax_table(physeq))
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq,"Phylum"))
  
  prevalenceThreshold = thresholdValue * nsamples(physeq)
  if(thresholdValue == 0)
    keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence) >= prevalenceThreshold]
  
  else if (prevalenceThreshold > 1)
    keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence) >= prevalenceThreshold]
  
  else
    keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence) >= 2]
  unique.Taxa <- unique(keepTaxa)
  return(unique.Taxa)
}

#GENERALIZED VERISION of getThresholdTaxaGen - 2019.08.16
getThresholdTaxaGen = function(physeq,thresholdValue){
  prevdf = apply(X = otu_table(physeq), 
                 MARGIN = ifelse(taxa_are_rows(physeq),yes=1, no = 2),
                 FUN = function(x){sum(x>0)})
  taxa.table = tax_table(physeq)
  prevdf = data.frame(Prevalence = prevdf,TotalAbundance = taxa_sums(physeq),taxa.table@.Data)
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq,"Phylum"))
  
  prevalenceThreshold = thresholdValue * nsamples(physeq)
  if(thresholdValue == 0)
    keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence) >= prevalenceThreshold]
  
  else if (prevalenceThreshold > 1)
    keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence) >= prevalenceThreshold]
  
  else
    keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence) >= 2]
  unique.Taxa <- unique(keepTaxa)
  return(unique.Taxa)
}


getThresholdTaxa2 = function(physeq,prevalenceThreshold){
  prevdf = apply(X = otu_table(physeq), 
                 MARGIN = ifelse(taxa_are_rows(physeq),yes=1, no = 2),
                 FUN = function(x){sum(x>0)})
  prevdf = data.frame(Prevalence = prevdf,TotalAbundance = taxa_sums(physeq),tax_table(physeq))
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq,"Phylum"))
  keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
  unique.Taxa <- unique(keepTaxa)
  return(unique.Taxa)
}

getThresholdTaxa3 = function(physeq,thresholdValue){
  prevdf = apply(X = otu_table(physeq), 
                 MARGIN = ifelse(taxa_are_rows(physeq),yes=1, no = 2),
                 FUN = function(x){sum(x>0)})
  prevdf = data.frame(Prevalence = prevdf,TotalAbundance = taxa_sums(physeq),tax_table(physeq))
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq,"Phylum"))
  
  prevalenceThreshold = thresholdValue * nsamples(physeq)
  keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
  return(keepTaxa)
}
plotPrevalence = function(physeq) {
  prevdf = apply(X = otu_table(physeq),
                 MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  # # Add taxonomy and total read counts to this data.frame
  prevdf = data.frame(Prevalence = prevdf,TotalAbundance = taxa_sums(physeq),tax_table(physeq))
  plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
  
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq, "Phylum"))
  ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(physeq),color=Phylum)) +
    # Include a guess for parameter
    scale_color_manual(values = microbiota.palette)+
    geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
    scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
    facet_wrap(~Phylum) + theme(legend.position="none")
}

plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Fusobacteria"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "DiseaseStatus",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}

MaaslinFormating <- function(data1,data2) {
  data1$Taxa <- paste(data1$Kingdom,data1$Phylum,data1$Class,data1$Order,data1$Family,data1$Genus,data1$Species,sep="|")
  data1 <- data1 %>% select (-c(rowname,Kingdom,Phylum,Class,Order,Family,Genus,Species))
  data1 <- data1 %>% select(Taxa, everything())
  data1 <- tibble::column_to_rownames(data1, var= "Taxa")
  
  data1 <- tibble::rownames_to_column(data1)
  data2 <- tibble::rownames_to_column(data2)
  
  data.maaslin.format <- dplyr::bind_rows(mutate_all(data2, as.character), mutate_all(data1, as.character))
  return(data.maaslin.format)
  
}

mergeTaxaRawOtuCounts <- function(physeq){
  taxa_of_physeq <- tax_table(physeq)
  otu_of_physeq <- t(otu_table(physeq))
  data_tbl_taxa <- as.data.frame(taxa_of_physeq)
  data_tbl_taxa <- tibble::rownames_to_column(data_tbl_taxa)
  data_tbl_otu <- as.data.frame(otu_of_physeq)
  data_tbl_otu <- tibble::rownames_to_column(data_tbl_otu)
  data_tbl_taxa_otu <- dplyr::full_join(data_tbl_taxa, data_tbl_otu, by = "rowname")
  return(data_tbl_taxa_otu)
}

mergeTaxaNormOtuCounts <- function(physeq){
  physeq <- transform_sample_counts(physeq,function(x) x/sum(x))
  taxa_of_physeq <- tax_table(physeq)
  otu_of_physeq <- t(otu_table(physeq))
  data_tbl_taxa <- as.data.frame(taxa_of_physeq)
  data_tbl_taxa <- tibble::rownames_to_column(data_tbl_taxa)
  data_tbl_otu <- as.data.frame(otu_of_physeq)
  data_tbl_otu <- tibble::rownames_to_column(data_tbl_otu)
  data_tbl_taxa_otu <- dplyr::full_join(data_tbl_taxa, data_tbl_otu, by = "rowname")
  return(data_tbl_taxa_otu)
}

MaaslinFormating <- function(data1,data2) {
  data1$Taxa <- paste(data1$Kingdom,data1$Phylum,data1$Class,data1$Order,data1$Family,data1$Genus,data1$Species,sep="|")
  data1 <- data1 %>% dplyr::select (-c(rowname,Kingdom,Phylum,Class,Order,Family,Genus,Species))
  data1 <- data1 %>% dplyr::select(Taxa, everything())
  data1 <- tibble::column_to_rownames(data1, var= "Taxa")
  
  data1 <- tibble::rownames_to_column(data1)
  data2 <- tibble::rownames_to_column(data2)
  
  data.maaslin.format <- dplyr::bind_rows(mutate_all(data2, as.character), mutate_all(data1, as.character))
  return(data.maaslin.format)
  
}

MaaslinFormatingGenus <- function(data1,data2) {
  data1$Taxa <- paste(data1$Kingdom,data1$Phylum,data1$Class,data1$Order,data1$Family,data1$Genus,sep="|")
  data1 <- data1 %>% select (-c(rowname,Kingdom,Phylum,Class,Order,Family,Genus))
  data1 <- data1 %>% select(Taxa, everything())
  data1 <- tibble::column_to_rownames(data1, var= "Taxa")
  
  data1 <- tibble::rownames_to_column(data1)
  data2 <- tibble::rownames_to_column(data2)
  
  data.maaslin.format <- dplyr::bind_rows(mutate_all(data2, as.character), mutate_all(data1, as.character))
  return(data.maaslin.format)
  
}

MaaslinFormating2 <- function(data) {
  data.norm.otu <- mergeTaxaNormOtuCounts(data)
  data.norm.otu <- as.data.frame(data.norm.otu)
  
  data.metadata <- sample_data(data)
  data.metadata.t <- as.data.frame(t(data.metadata))
  
  data.maaslin.format <- MaaslinFormating(data.norm.otu,data.metadata.t)
  return(data.maaslin.format)
}

MaaslinFormatingGenus2 <- function(data) {
  data.norm.otu <- mergeTaxaNormOtuCounts(data)
  data.norm.otu <- as.data.frame(data.norm.otu)
  
  data.metadata <- sample_data(data)
  data.metadata.t <- as.data.frame(t(data.metadata))
  
  data.maaslin.format <- MaaslinFormatingGenus(data.norm.otu,data.metadata.t)
  return(data.maaslin.format)
}

# MaaslinFormating3 <- function(data) {
#   data.norm.otu <- mergeTaxaRawOtuCounts(data)
#   data.norm.otu <- as.data.frame(data.norm.otu)
#   
#   data.metadata <- sample_data(data)
#   data.metadata.t <- as.data.frame(t(data.metadata))
#   
#   data.maaslin.format <- MaaslinFormating(data.norm.otu,data.metadata.t)
#   return(data.maaslin.format)
# }


test.plot <- plotPrevalence(ps.core.062)
genus.10 = getThresholdGenus(test.phy,0.10)
test.phy10= subset_taxa(test.phy, Genus %in% genus.10)
test.phy10  = prune_samples(sample_sums(test.phy10 )>0, test.phy10 )

#ADDED 2022.01.08
MaaslinFormating2 <- function(data) {
  data.norm.otu <- mergeTaxaNormOtuCounts(data)
  data.norm.otu <- as.data.frame(data.norm.otu)
  
  data.metadata <- sample_data(data)
  data.metadata.t <- as.data.frame(t(data.metadata))
  
  data.maaslin.format <- MaaslinFormating(data.norm.otu,data.metadata.t)
  return(data.maaslin.format)
}

LefseFormat <- function(data,data_key){
  data_t <- t(data)
  data_t <- janitor::row_to_names(data_t,row_number = 1)
  data_t_key <- merge(data_key,data_t,by.x="pathway",by.y=0, all.y=TRUE)
  return(data_t_key)
}

#fiddling with Seurat Object

# SeuratObject@assay$RNA ---> original RNA counts
# SeuratObject@assay$RNA[m_gene:n_gene,x_cell:y_cell] --> Gives Slice of count matrix
# colSums(SeuratObject@assay$RNA) ---> Gives number of counts

# SeuratObject.transcripts <- rownames(SeuratObject@assay$RNA) ==> gives names of transcripts
#
# viral.sums <- colSums(SeuratObject@assays$RNA[grepl("^NC-",SeuratObject.transcripts),]) ==> Gives sum of viral transcripts across cells
# SeuratObject.viral.matrix <- SeuratObject@assays$RNA[grepl("^NC-",SeuratObject.transcripts),] => get only viral counts
# SeuratObject.viral.matrix[SeuratObject.viral.matrix != 0] <- 1 => converts to every non-zero value to 1
# (summing the above gives number of viral FEATURES)
data[!grepl("^ca", samps)]

S1.transcripts <- rownames(S1@assays$RNA)
#S1.transcripts[grepl("^NC-",S1.transcripts)]
#S1.vir.only <- S1.transcripts[grepl("^NC-",S1.transcripts)]

S1.viral.matrix <- S1@assays$RNA[grepl("^NC-",S1.transcripts),]
viral.sums <- colSums(S1@assays$RNA[grepl("^NC-",S1.transcripts),])
S1.viral.matrix[S1.viral.matrix != 0] <- 1
viral.features <- colSums(S1.viral.matrix)

S1$nCount_Viral <- viral.sums
S1$nFeature_Viral <- viral.features

S1$nFeature_RNA_human <- S1$nFeature_RNA - S1$nFeature_Viral
S1$nCount_RNA_human <- S1$nCount_RNA - S1$nCount_Viral

S1.mito.matrix <- S1@assays$RNA[grepl("^MT-",S1.transcripts),]
mito.sums <- colSums(S1.mito.matrix)
S1.mito.matrix[S1.mito.matrix != 0] <- 1
mito.features <- colSums(S1.mito.matrix)

S1$nCount_mito <- mito.sums
S1$nFeature_mito <- mito.features
S1$percent.mt.adj <- (S1$nCount_mito/S1$nCount_RNA_human)*100
#random_group_labels <- sample(x = c("g1", "g2"), size = ncol(x = S1), replace = TRUE)


#my_list <- list(l1 = c(1, 3, 5, 7),                
l2 = c(1, 2, 3),                    
l3 = c(1, 1, 10, 5, 8, 65, 90)) 
#S2
S2.transcripts <- rownames(S2@assays$RNA)
#S2.transcripts[grepl("^NC-",S2.transcripts)]
#S2.vir.only <- S2.transcripts[grepl("^NC-",S2.transcripts)]

S2.viral.matrix <- S2@assays$RNA[grepl("^NC-",S2.transcripts),]
viral.sums <- colSums(S2@assays$RNA[grepl("^NC-",S2.transcripts),])
S2.viral.matrix[S2.viral.matrix != 0] <- 1
viral.features <- colSums(S2.viral.matrix)

S2$nCount_Viral <- viral.sums
S2$nFeature_Viral <- viral.features

S2$nFeature_RNA_human <- S2$nFeature_RNA - S2$nFeature_Viral
S2$nCount_RNA_human <- S2$nCount_RNA - S2$nCount_Viral

S2.mito.matrix <- S2@assays$RNA[grepl("^MT-",S2.transcripts),]
mito.sums <- colSums(S2.mito.matrix)
S2.mito.matrix[S2.mito.matrix != 0] <- 1
mito.features <- colSums(S2.mito.matrix)

S2$nCount_mito <- mito.sums
S2$nFeature_mito <- mito.features
S2$percent.mt.adj <- (S2$nCount_mito/S2$nCount_RNA_human)*100

#S2
S3.transcripts <- rownames(S3@assays$RNA)
#S3.transcripts[grepl("^NC-",S3.transcripts)]
#S3.vir.only <- S3.transcripts[grepl("^NC-",S3.transcripts)]

S3.viral.matrix <- S3@assays$RNA[grepl("^NC-",S3.transcripts),]
viral.sums <- colSums(S3@assays$RNA[grepl("^NC-",S3.transcripts),])
S3.viral.matrix[S3.viral.matrix != 0] <- 1
viral.features <- colSums(S3.viral.matrix)

S3$nCount_Viral <- viral.sums
S3$nFeature_Viral <- viral.features

S3$nFeature_RNA_human <- S3$nFeature_RNA - S3$nFeature_Viral
S3$nCount_RNA_human <- S3$nCount_RNA - S3$nCount_Viral

S3.mito.matrix <- S3@assays$RNA[grepl("^MT-",S3.transcripts),]
mito.sums <- colSums(S3.mito.matrix)
S3.mito.matrix[S3.mito.matrix != 0] <- 1
mito.features <- colSums(S3.mito.matrix)

S3$nCount_mito <- mito.sums
S3$nFeature_mito <- mito.features
S3$percent.mt.adj <- (S3$nCount_mito/S3$nCount_RNA_human)*100

#S4
S4.transcripts <- rownames(S4@assays$RNA)
#S4.transcripts[grepl("^NC-",S4.transcripts)]
#S4.vir.only <- S4.transcripts[grepl("^NC-",S4.transcripts)]

S4.viral.matrix <- S4@assays$RNA[grepl("^NC-",S4.transcripts),]
viral.sums <- colSums(S4@assays$RNA[grepl("^NC-",S4.transcripts),])
S4.viral.matrix[S4.viral.matrix != 0] <- 1
viral.features <- colSums(S4.viral.matrix)

S4$nCount_Viral <- viral.sums
S4$nFeature_Viral <- viral.features

S4$nFeature_RNA_human <- S4$nFeature_RNA - S4$nFeature_Viral
S4$nCount_RNA_human <- S4$nCount_RNA - S4$nCount_Viral

S4.mito.matrix <- S4@assays$RNA[grepl("^MT-",S4.transcripts),]
mito.sums <- colSums(S4.mito.matrix)
S4.mito.matrix[S4.mito.matrix != 0] <- 1
mito.features <- colSums(S4.mito.matrix)

S4$nCount_mito <- mito.sums
S4$nFeature_mito <- mito.features
S4$percent.mt.adj <- (S4$nCount_mito/S4$nCount_RNA_human)*100

#S5
S5.transcripts <- rownames(S5@assays$RNA)
#S5.transcripts[grepl("^NC-",S5.transcripts)]
#S5.vir.only <- S5.transcripts[grepl("^NC-",S5.transcripts)]

S5.viral.matrix <- S5@assays$RNA[grepl("^NC-",S5.transcripts),]
viral.sums <- colSums(S5@assays$RNA[grepl("^NC-",S5.transcripts),])
S5.viral.matrix[S5.viral.matrix != 0] <- 1
viral.features <- colSums(S5.viral.matrix)

S5$nCount_Viral <- viral.sums
S5$nFeature_Viral <- viral.features

S5$nFeature_RNA_human <- S5$nFeature_RNA - S5$nFeature_Viral
S5$nCount_RNA_human <- S5$nCount_RNA - S5$nCount_Viral

S5.mito.matrix <- S5@assays$RNA[grepl("^MT-",S5.transcripts),]
mito.sums <- colSums(S5.mito.matrix)
S5.mito.matrix[S5.mito.matrix != 0] <- 1
mito.features <- colSums(S5.mito.matrix)

S5$nCount_mito <- mito.sums
S5$nFeature_mito <- mito.features
S5$percent.mt.adj <- (S5$nCount_mito/S5$nCount_RNA_human)*100

#S6
S6.transcripts <- rownames(S6@assays$RNA)
#S6.transcripts[grepl("^NC-",S6.transcripts)]
#S6.vir.only <- S6.transcripts[grepl("^NC-",S6.transcripts)]

S6.viral.matrix <- S6@assays$RNA[grepl("^NC-",S6.transcripts),]
viral.sums <- colSums(S6@assays$RNA[grepl("^NC-",S6.transcripts),])
S6.viral.matrix[S6.viral.matrix != 0] <- 1
viral.features <- colSums(S6.viral.matrix)

S6$nCount_Viral <- viral.sums
S6$nFeature_Viral <- viral.features

S6$nFeature_RNA_human <- S6$nFeature_RNA - S6$nFeature_Viral
S6$nCount_RNA_human <- S6$nCount_RNA - S6$nCount_Viral

S6.mito.matrix <- S6@assays$RNA[grepl("^MT-",S6.transcripts),]
mito.sums <- colSums(S6.mito.matrix)
S6.mito.matrix[S6.mito.matrix != 0] <- 1
mito.features <- colSums(S6.mito.matrix)

S6$nCount_mito <- mito.sums
S6$nFeature_mito <- mito.features
S6$percent.mt.adj <- (S6$nCount_mito/S6$nCount_RNA_human)*10
########################
tgen6 <- readRDS(file = "location")

#Get viral transcripts seen
#Step One: Get full list of transcripts seen 
tgen6.rows <- rownames(x=tgen6)

#Step Two: Get starting position of viral transcripts
head(grep(pattern="NC-",tgen6.rows))
#Step Three: Get end position of viral transcripts seen
tail(grep(pattern="NC-",tgen6.rows))

#Get # of viral transcripts seen
vir <- grep("NC-",tgen6.rows)
length(vir) 

virus_counter <- function(seurat_object) {
  seurat.rows <- rownames(x=seurat_object)
  vir.ids <- grep("^NC-",seurat.rows)
  vir.counts <- length(vir.ids)
  print(vir.counts)
}

PercentageMito <- function(seurat_object) {
  seurat_object.transcripts <- rownames(seurat_object@assays$RNA)
  seurat_object.mito.matrix <- seurat_object@assays$RNA[grepl("^MT-",S1.transcripts),]
  mito.sums <- colSums(seurat_object.mito.matrix)
  seurat_object.mito.matrix[seurat_object.mito.matrix != 0] <- 1
  mito.features <- colSums(seurat_object.mito.matrix)
  
  seurat_object$nCount_mito <- mito.sums
  seurat_object$nFeature_mito <- mito.features
}

#####
library(UpSetR)
library(ggplot2)

###MY DATA#####
COPD_1_mem <- read.csv("C:/Users/nmshahir/Dropbox/Davenport_Postdoc/scRNA_viral_reads/kumata_ref/COPD/IPF_TILD063_1_LF_Whole_C1_X5SCR_F01825_HN3KNDMXX_viral_map_counts.csv")
COPD_2_mem <- read.csv("C:/Users/nmshahir/Dropbox/Davenport_Postdoc/scRNA_viral_reads/kumata_ref/COPD/IPF_TILD063_1_MF_Whole_C2_X5SCR_F01826_HN3KNDMXX_viral_map_counts.csv")
Ctrl_1_mem <- read.csv("C:/Users/nmshahir/Dropbox/Davenport_Postdoc/scRNA_viral_reads/kumata_ref/Control/IPF_VUHD65_1_LU_Whole_C1_X5SCR_HD65_H5LLNDSXX_viral_map_counts.csv")
Ctrl_2_mem <- read.csv("C:/Users/nmshahir/Dropbox/Davenport_Postdoc/scRNA_viral_reads/kumata_ref/Control/IPF_VUHD70_1_LU_Whole_C1_X5SCR_HD70_HCKWNDSXX_viral_map_counts.csv")
IPF_1_mem <- read.csv("C:/Users/nmshahir/Dropbox/Davenport_Postdoc/scRNA_viral_reads/kumata_ref/IPF/IPF_VUILD60_1_LU_Whole_C1_X5SCR_ILD60-1_H5LM5DSXX_viral_map_counts.csv")
IPF_2_mem <- read.csv("C:/Users/nmshahir/Dropbox/Davenport_Postdoc/scRNA_viral_reads/kumata_ref/IPF/IPF_VUILD60_1_LU_Whole_C1_X5SCR_ILD60-2_H5LM5DSXX_viral_map_counts.csv")


Ctrl.1.kumata <- merge(Ctrl_1_mem,kumata_ref_list,by="Refseq",all.x = TRUE)
Ctrl.2.kumata <- merge(Ctrl_2_mem,kumata_ref_list,by="Refseq",all.x = TRUE)

COPD.1.kumata <- merge(COPD_1_mem,kumata_ref_list,by="Refseq",all.x=TRUE)
COPD.2.kumata <- merge(COPD_2_mem,kumata_ref_list,by="Refseq",all.x=TRUE)

IPF.1.kumata <- merge(IPF_1_mem,kumata_ref_list,by="Refseq",all.x = TRUE)
IPF.2.kumata <- merge(IPF_2_mem,kumata_ref_list,by="Refseq",all.x = TRUE)

####ALL SAMPLES####

kumata_virus <- list(COPD.1=COPD.1.kumata$Refseq,COPD.2=COPD.2.kumata$Refseq,Ctrl.1=Ctrl.1.kumata$Refseq,Ctrl.2=Ctrl.2.kumata$Refseq,IPF.1=IPF.1.kumata$Refseq,IPF.2=IPF.2.kumata$Refseq)
upset(fromList(kumata_virus),nsets=6,order.by = "freq",
      number.angles = 30,
      point.size = 3.5,
      line.size = 2,
      mainbar.y.label = "Viruses Aligned from Kumata et al Reference",
      sets.x.label = "Total Number of Distinct Viruses Seen via BWA mem") 


#### Log the reads###
COPD_1_mem$lReads <- log10(COPD_1_mem$Reads)
COPD_2_mem$lReads <- log10(COPD_2_mem$Reads)
Ctrl_1_mem$lReads <- log10(Ctrl_1_mem$Reads)
Ctrl_2_mem$lReads <- log10(Ctrl_2_mem$Reads)
IPF_1_mem$lReads <- log10(IPF_1_mem$Reads)
IPF_2_mem$lReads <- log10(IPF_2_mem$Reads)

#### Control Historgrams ###
ggplot(Ctrl_1_mem, aes(x=lReads))+geom_histogram() +
  labs(x="log10 of Reads", y = "Frequency",title="Control 1")
ggplot(Ctrl_2_mem, aes(x=lReads))+geom_histogram() +
  labs(x="log10 of Reads", y = "Frequency",title="Control 2")

#### COPD Historgrams ###
ggplot(COPD_1_mem, aes(x=lReads))+geom_histogram() +
  labs(x="log10 of Reads", y = "Frequency",title="COPD Less Fibrosis")
ggplot(COPD_2_mem, aes(x=lReads))+geom_histogram() +
  labs(x="log10 of Reads", y = "Frequency",title="COPD More Fibrosis")

#### IPF Historgrams ###
ggplot(IPF_1_mem, aes(x=lReads))+geom_histogram() +
  labs(x="log10 of Reads", y = "Frequency",title="IPF Fibrosis 1")
ggplot(IPF_2_mem, aes(x=lReads))+geom_histogram() +
  labs(x="log10 of Reads", y = "Frequency",title="IPF Fibrosis 2")
#####
library(dplyr)
library(tibble)
library(Seurat)
library(patchwork)
library(ggplot2)

setwd("C:/Users/nmshahir/Dropbox/Davenport_Postdoc/Seuart_Object/filtered_feature_bc_matrix")

#Control
d10x <- Read10X("Ctrl1/filtered_feature_bc_matrix/")
colnames(d10x) <- paste("VUHD65", sep="_", colnames(d10x))



ctrl.1 <- CreateSeuratObject(d10x, min.features = 1, names.field = 1, names.delim = "_")
rm(d10x)

d10x <- Read10X("Ctrl2/filtered_feature_bc_matrix/")
colnames(d10x) <- paste("VUHD70", sep="_", colnames(d10x))



ctrl.2 <- CreateSeuratObject(d10x, min.features = 1, names.field = 1, names.delim = "_")
rm(d10x)

#QC
ctrl.1 <- PercentageFeatureSet(object = ctrl.1, pattern = "^MT-", col.name = "percent.mt")
ctrl.1 <- subset(ctrl.1, subset = nFeature_RNA > 1000 & percent.mt < 25)

ctrl.2 <- PercentageFeatureSet(object = ctrl.2, pattern = "^MT-", col.name = "percent.mt")
ctrl.2 <- subset(ctrl.2, subset = nFeature_RNA > 1000 & percent.mt < 25)
#Get Viral Counts
ctrl1.rows <- rownames(x=ctrl.1)
ctrl2.rows <- rownames(x=ctrl.2)

vir.ctrl1 <- grep("NC-",ctrl1.rows)
vir.ctrl2 <- grep("NC-",ctrl2.rows)

#COPD
d10x <- Read10X("COPD1/filtered_feature_bc_matrix/")
colnames(d10x) <- paste("F01825", sep="_", colnames(d10x))



copd.1 <- CreateSeuratObject(d10x, min.features = 1, names.field = 1, names.delim = "_")
rm(d10x)

d10x <- Read10X("COPD2/filtered_feature_bc_matrix/")
colnames(d10x) <- paste("F01826", sep="_", colnames(d10x))



copd.2 <- CreateSeuratObject(d10x, min.features = 1, names.field = 1, names.delim = "_")
rm(d10x)

#QC
copd.1 <- PercentageFeatureSet(object = copd.1, pattern = "^MT-", col.name = "percent.mt")
copd.1 <- subset(copd.1, subset = nFeature_RNA > 1000 & percent.mt < 25)

copd.2 <- PercentageFeatureSet(object = copd.2, pattern = "^MT-", col.name = "percent.mt")
copd.2 <- subset(copd.2, subset = nFeature_RNA > 1000 & percent.mt < 25)

copd1.rows <- row.names(x=copd.1)
copd2.rows <- row.names(x=copd.2)

vir.copd1 <- grep("NC-",copd1.rows)
vir.copd2 <- grep("NC-",copd2.rows)

#IPF
d10x <- Read10X("IPF1/filtered_feature_bc_matrix/")
colnames(d10x) <- paste("ILD60-1", sep="_", colnames(d10x))



ipf.1 <- CreateSeuratObject(d10x, min.features = 1, names.field = 1, names.delim = "_")
rm(d10x)

d10x <- Read10X("IPF2/filtered_feature_bc_matrix/")
colnames(d10x) <- paste("ILD60-2", sep="_", colnames(d10x))



ipf.2 <- CreateSeuratObject(d10x, min.features = 1, names.field = 1, names.delim = "_")
rm(d10x)

#QC
ipf.1 <- PercentageFeatureSet(object = ipf.1, pattern = "^MT-", col.name = "percent.mt")
ipf.1 <- subset(ipf.1, subset = nFeature_RNA > 1000 & percent.mt < 25)

ipf.2 <- PercentageFeatureSet(object = ipf.2, pattern = "^MT-", col.name = "percent.mt")
ipf.2 <- subset(ipf.2, subset = nFeature_RNA > 1000 & percent.mt < 25)

ipf1.rows <- row.names(x=ipf.1)
ipf2.rows <- row.names(x=ipf.2)

vir.ipf1 <- grep("NC-",ipf1.rows)
vir.ipf2 <- grep("NC-",ipf2.rows)

####################################
library(devtools)
install_github("ropensci/rentrez")
library(rentrez)
entrez_dbs()

#Get nuccore id
r_search <- entrez_search(db="nuccore", term="NC_000898.1")
#Get summary
taxize_summ <- entrez_summary(db="nuccore", id=r_search$ids)
#Get name
taxize_summ$title
#Get source format
taxize_summ$subtype
#Get actual source
taxize_summ$subname

#################################################################
#NCBI BWA MEM RUN
##COPD Patient 1 
COPD1.ncbi.mem.list <- c()
for(ref in 1:length(COPD_1_mem$Reference)) {
  COPD1.ncbi.mem.list <- append(COPD1.ncbi.mem.list,entrez_search(db="nuccore",term=COPD_1_mem$Reference[ref]))
}
#COPD_1_mem.search <- entrez_search(db="nuccore",term=COPD_1_mem$Reference)
entrez_summary(db="pubmed", id="NC_000898.1")

#################################################################
#### Import required Libraries####
library("vegan")
library("MASS")
library("phyloseq")
library("ggplot2")
library("dada2")
library("ggridges")
library("ggpubr")
library("gridExtra")
library("microbiome")
####Pre-processing ####
# Import data files

#If on Windows
setwd("C:/Dropbox/IBD_microbiota/Analyses/Metadata")


ps.patients.g <- readRDS("ps_patients_species_20190713.rds")
metadata.func <- as.data.frame(sample_data(ps.patients.g))
rownames(metadata.func) <- gsub(" ","_",rownames(metadata.func))
pred_metagenome_unstrat <- read.delim("C:/Users/nmshahir/Dropbox/IBD_microbiota/Analyses/Downstream_Analyses/picrust2_out_pipeline_take2/KO_metagenome_out/pred_metagenome_unstrat.tsv", header=TRUE)
colnames(pred_metagenome_unstrat) <- gsub("X","",colnames(pred_metagenome_unstrat))
names(pred_metagenome_unstrat)[1] <- "Taxa"

pred_metagenome_unstrat_descrip <- read.delim("C:/Users/nmshahir/Dropbox/IBD_microbiota/Analyses/Downstream_Analyses/picrust2_out_pipeline_take2/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv")
names(pred_metagenome_unstrat_descrip)[1] <- "Taxa"
pred_test.names <- pred_metagenome_unstrat$Taxa
pred_test <- as.data.frame(pred_metagenome_unstrat)
pred.test.test <- data.matrix(pred_test)
pred.test.test[,1] <- pred_test.names
rownames(pred.test.test) <- pred_test.names
#func_description <- pred_metagenome_unstrat_descrip[c("Taxa","description")]

#phyloseq object with functional data
ps.func.all<- phyloseq( metadata.func,otu_table(pred.test.test, taxa_are_rows =TRUE))

############################################################################################
ps_clr <- microbiome::transform(ps.func.all,"clr")
ps.clr.unc <- subset_samples(ps_clr,sample_data(ps_clr)$Cohort == "UNC")
ps.clr.unc <- prune_taxa(taxa_sums(ps.clr.unc)>0,ps.clr.unc)

ps.uc <- subset_samples(ps.clr.unc,sample_data(ps.clr.unc)$DiseaseStatus == "UC")
ps.uc <- prune_taxa(taxa_sums(ps.uc)>0, ps.uc)

ps.uc.colon <- subset_samples(ps.uc,sample_data(ps.uc)$AnatomSite == "colon")
ps.uc.n.colon <- subset_samples(ps.uc.colon,sample_data(ps.uc.colon)$Pathology == "NI")
ps.uc.n.colon <- prune_taxa(taxa_sums(ps.uc.n.colon)>0, ps.uc.n.colon)

sample_data(ps.uc.n.colon)$AeroLvl <- ifelse(sample_data(ps.uc.n.colon)$PatientNo %in% c("127","41500","42700","43000","46","52"),"Low","Health")
ps.uc.n.colon.clean <- subset_samples(ps.uc.n.colon,sample_data(ps.uc.n.colon)$Antibiotics != 1)
ps.uc.n.colon.clean <- prune_taxa(taxa_sums(ps.uc.n.colon.clean)>0,ps.uc.n.colon.clean)

ps.uc.n.colon.correct <- subset_samples(ps.uc.n.colon.clean,sample_data(ps.uc.n.colon.clean)$PatientNo != "18")
ps.uc.n.colon.correct <- prune_taxa(taxa_sums(ps.uc.n.colon.correct)>0,ps.uc.n.colon.correct)
UC.clr <- ps.uc.n.colon.correct
ord_clr <- phyloseq::ordinate(UC.clr, "RDA")
#Examine eigenvalues and % prop. variance explained
head(ord_clr$CA$eig)
sapply(ord_clr$CA$eig[1:8], function(x) x / sum(ord_clr$CA$eig))


ind.coord <- data.frame(ord_clr$CA$v)
sdev_ind <- apply(ind.coord, 1, sd)
ind_cont_PCA1 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC1^2 /sdev_ind))))
ind_cont_PCA1_top <- ind_cont_PCA1 %>% 
  rownames_to_column("Taxa") %>% 
  filter(PCA >= 0.001803355)%>% 
  column_to_rownames("Taxa")
ind_cont_PCA2 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC2^2 /sdev_ind))))
ind_cont_PCA2_top <- ind_cont_PCA2 %>% 
  rownames_to_column("otu") %>% 
  filter(PCA >= 0.002)%>% 
  column_to_rownames("otu")
ind_cont_PCA3 <- data.frame(PCA = (100*(1 / nrow(ind.coord)*(ind.coord$PC3^2 / sdev_ind))))
ind_cont_PCA3_top <- ind_cont_PCA3 %>% 
  rownames_to_column("otu") %>% 
  filter(PCA >= 0.002)%>% 
  column_to_rownames("otu")
ind_cont_PCA_top <- rbind(ind_cont_PCA1_top, ind_cont_PCA2_top, ind_cont_PCA3_top)

sum(ind_cont_PCA1_top$PCA) / sum(ind_cont_PCA1$PCA)
##PLOT PCA

p2 = plot_ordination(UC.clr, ord_clr, type="samples", color="AeroLvl") 
plot_ordination(ps.uc,n,colon, ord_clr, type="samples", color="Cohort", shape="AnatomSite",title="all")
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(UC.clr, ord_clr, type="samples", color="AeroLvl") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = AeroLvl), linetype = 2)

clr_dist_matrix <- phyloseq::distance(UC.clr, method = "euclidean") 
vegan::adonis(clr_dist_matrix ~ phyloseq::sample_data(UC.clr)$AeroLvl)
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(UC.clr)$AeroLvl)



################
ps.clr.unc <- subset_samples(ps.func.all,sample_data(ps.func.all)$Cohort == "UNC")
ps.clr.unc <- prune_taxa(taxa_sums(ps.clr.unc)>0,ps.clr.unc)

ps.uc <- subset_samples(ps.clr.unc,sample_data(ps.clr.unc)$DiseaseStatus == "UC")
ps.uc <- prune_taxa(taxa_sums(ps.uc)>0, ps.uc)

ps.uc.colon <- subset_samples(ps.uc,sample_data(ps.uc)$AnatomSite == "colon")
ps.uc.n.colon <- subset_samples(ps.uc.colon,sample_data(ps.uc.colon)$Pathology == "NI")
ps.uc.n.colon <- prune_taxa(taxa_sums(ps.uc.n.colon)>0, ps.uc.n.colon)

sample_data(ps.uc.n.colon)$AeroLvl <- ifelse(sample_data(ps.uc.n.colon)$PatientNo %in% c("127","41500","42700","43000","46","52"),"Low","Health")
ps.uc.n.colon.clean <- subset_samples(ps.uc.n.colon,sample_data(ps.uc.n.colon)$Antibiotics != 1)
ps.uc.n.colon.clean <- prune_taxa(taxa_sums(ps.uc.n.colon.clean)>0,ps.uc.n.colon.clean)

ps.uc.n.colon.correct <- subset_samples(ps.uc.n.colon.clean,sample_data(ps.uc.n.colon.clean)$PatientNo != "18")
ps.uc.n.colon.correct <- prune_taxa(taxa_sums(ps.uc.n.colon.correct)>0,ps.uc.n.colon.correct)
UC.clr <- ps.uc.n.colon.correct


aldex2_da <- ALDEx2::aldex(round(data.frame(phyloseq::otu_table(UC.clr))), phyloseq::sample_data(UC.clr)$AeroLvl, test="t", effect = TRUE, denom="iqlr")
sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = "KEGG") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(KEGG, diff.btw, diff.win, effect, wi.ep, wi.eBH)
sig_aldex2 <- left_join(sig_aldex2, taxa_info)
###################################################################################################################################################################
#EC data
#######
pred_metagenome_unstrat <- read.delim("C:/Users/nmshahir/Dropbox/IBD_microbiota/Analyses/Downstream_Analyses/picrust2_out_pipeline_take2/EC_metagenome_out/pred_metagenome_unstrat.tsv", header=TRUE)
colnames(pred_metagenome_unstrat) <- gsub("X","",colnames(pred_metagenome_unstrat))
names(pred_metagenome_unstrat)[1] <- "Taxa"

pred_test.names <- pred_metagenome_unstrat$Taxa
pred_test <- as.data.frame(pred_metagenome_unstrat)
pred.test.test <- data.matrix(pred_test)

rownames(pred.test.test) <- pred_test.names
#func_description <- pred_metagenome_unstrat_descrip[c("Taxa","description")]

#phyloseq object with functional data
ps.func.all<- phyloseq( metadata.func,otu_table(pred.test.test, taxa_are_rows =TRUE))

ps.clr.unc <- subset_samples(ps.func.all,sample_data(ps.func.all)$Cohort == "UNC")
ps.clr.unc <- prune_taxa(taxa_sums(ps.clr.unc)>0,ps.clr.unc)

ps.uc <- subset_samples(ps.clr.unc,sample_data(ps.clr.unc)$DiseaseStatus == "UC")
ps.uc <- prune_taxa(taxa_sums(ps.uc)>0, ps.uc)

ps.uc.colon <- subset_samples(ps.uc,sample_data(ps.uc)$AnatomSite == "colon")
ps.uc.n.colon <- subset_samples(ps.uc.colon,sample_data(ps.uc.colon)$Pathology == "NI")
ps.uc.n.colon <- prune_taxa(taxa_sums(ps.uc.n.colon)>0, ps.uc.n.colon)

sample_data(ps.uc.n.colon)$AeroLvl <- ifelse(sample_data(ps.uc.n.colon)$PatientNo %in% c("127","41500","42700","43000","46","52"),"Low","Health")
ps.uc.n.colon.clean <- subset_samples(ps.uc.n.colon,sample_data(ps.uc.n.colon)$Antibiotics != 1)
ps.uc.n.colon.clean <- prune_taxa(taxa_sums(ps.uc.n.colon.clean)>0,ps.uc.n.colon.clean)

ps.uc.n.colon.correct <- subset_samples(ps.uc.n.colon.clean,sample_data(ps.uc.n.colon.clean)$PatientNo != "18")
ps.uc.n.colon.correct <- prune_taxa(taxa_sums(ps.uc.n.colon.correct)>0,ps.uc.n.colon.correct)
UC.clr <- ps.uc.n.colon.correct


aldex2_da <- ALDEx2::aldex(round(data.frame(phyloseq::otu_table(UC.clr))), phyloseq::sample_data(UC.clr)$AeroLvl, test="t", effect = TRUE, denom="iqlr")
sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = "EC") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(EC, diff.btw, diff.win, effect, wi.ep, wi.eBH)
sig_aldex2 <- left_join(sig_aldex2, taxa_info)

###################################################################################################################################################################
#Metacyc data
#######
pred_metagenome_unstrat <- read.delim("C:/Users/nmshahir/Dropbox/IBD_microbiota/Analyses/Downstream_Analyses/picrust2_out_pipeline_take2/pathways_out/path_abun_unstrat.tsv", header=TRUE)
colnames(pred_metagenome_unstrat) <- gsub("X","",colnames(pred_metagenome_unstrat))
names(pred_metagenome_unstrat)[1] <- "Taxa"

pred_test.names <- pred_metagenome_unstrat$Taxa
pred_test <- as.data.frame(pred_metagenome_unstrat)
pred.test.test <- data.matrix(pred_test)

rownames(pred.test.test) <- pred_test.names
#func_description <- pred_metagenome_unstrat_descrip[c("Taxa","description")]

#phyloseq object with functional data
ps.func.all<- phyloseq( metadata.func,otu_table(pred.test.test, taxa_are_rows =TRUE))

ps.clr.unc <- subset_samples(ps.func.all,sample_data(ps.func.all)$Cohort == "UNC")
ps.clr.unc <- prune_taxa(taxa_sums(ps.clr.unc)>0,ps.clr.unc)

ps.uc <- subset_samples(ps.clr.unc,sample_data(ps.clr.unc)$DiseaseStatus == "UC")
ps.uc <- prune_taxa(taxa_sums(ps.uc)>0, ps.uc)

ps.uc.colon <- subset_samples(ps.uc,sample_data(ps.uc)$AnatomSite == "colon")
ps.uc.n.colon <- subset_samples(ps.uc.colon,sample_data(ps.uc.colon)$Pathology == "NI")
ps.uc.n.colon <- prune_taxa(taxa_sums(ps.uc.n.colon)>0, ps.uc.n.colon)

sample_data(ps.uc.n.colon)$AeroLvl <- ifelse(sample_data(ps.uc.n.colon)$PatientNo %in% c("127","41500","42700","43000","46","52"),"Low","Health")
ps.uc.n.colon.clean <- subset_samples(ps.uc.n.colon,sample_data(ps.uc.n.colon)$Antibiotics != 1)
ps.uc.n.colon.clean <- prune_taxa(taxa_sums(ps.uc.n.colon.clean)>0,ps.uc.n.colon.clean)

ps.uc.n.colon.correct <- subset_samples(ps.uc.n.colon.clean,sample_data(ps.uc.n.colon.clean)$PatientNo != "18")
ps.uc.n.colon.correct <- prune_taxa(taxa_sums(ps.uc.n.colon.correct)>0,ps.uc.n.colon.correct)
UC.clr <- ps.uc.n.colon.correct


aldex2_da <- ALDEx2::aldex(round(data.frame(phyloseq::otu_table(UC.clr))), phyloseq::sample_data(UC.clr)$AeroLvl, test="t", effect = TRUE, denom="iqlr")
sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = "EC") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(EC, diff.btw, diff.win, effect, wi.ep, wi.eBH)
sig_aldex2 <- left_join(sig_aldex2, taxa_info)

####################################################
##functions for processing compositional microbiome datasets
#all scripts written by BPB, 062519

#ERROR PROTECTION FUNCTIONS
format.ASV.tab <- function(ps){
  if(as.logical(class(phyloseq::otu_table(ps))[1] == "otu_table") && 
     as.logical(taxa_are_rows(phyloseq::otu_table(ps)) == TRUE)){
    asv.tab <- as.matrix(phyloseq::otu_table(ps))
  } else {
    asv.tab <- as.matrix(t(phyloseq::otu_table(ps)))
  }
  asv.tab
}

format.parameter.string <- function(string){
  string <- gsub(pattern = "c(|)", replacement = "", x = string)
  string[2] <- gsub(pattern = ":", replacement = ", ", x = string[2])
  string <- strsplit(string, split = ", ")
  string <- as.numeric(paste0(c(string[[2]], string[[3]])))
  string
}

format.long <- function(df){
  ctc <- rep(colnames(df)[2], nrow(df))
  ctm <- rep(colnames(df)[3], nrow(df))
  tvv <- rep(df[,1], 2)
  th <- c(df[,2], df[, 3])
  
  df2 <- cbind.data.frame(c(ctc, ctm), tvv, th)
  colnames(df2) <- c("cgroup", "threshold.value", "taxa.hits")
  df2
}

#PLOTTING FUNCTIONS
plot.threshold <- function(est.obj, y=NULL, x=NULL, PFT=NULL, RAT=NULL, CVT=NULL, taxrank=NULL){
  #never gray
  ggplot2::theme_set(theme_bw())
  
  #error checking
  if (!is.null(x) && is.null(y)){
    stop("Please provide variables for both axes")
  }
  
  if (is.null(x) && !is.null(y)){
    stop("Please provide variables for both axes")
  }
  
  #plot either WS, AS or taxonomic filter figs
  
  #begin WS threshold plots
  if (is.data.frame(est.obj)){
    df <- format.long(est.obj)
    
    plot1 <- ggplot2::ggplot(data = df, aes(x=threshold.value, y=taxa.hits, color=cgroup)) + 
      ggplot2::geom_line(size=2, alpha=0.7) + ggplot2::scale_color_manual(values = c("orange", "steelblue2"), labels = c("Total", "Matches")) +
      ggplot2::theme(axis.text.y = element_text(size = 15, colour = "black"),
                     axis.text.x = element_text(size = 15, colour = "black", angle = 315, vjust = 0.7),
                     axis.title.x = element_text(size = 15, colour = "black"),
                     axis.title.y = element_text(size = 15),
                     legend.title = element_text(size = 0),
                     legend.text = element_text(size = 15, colour = "black"),
                     legend.position = "top",
                     legend.key = element_rect(size = 5),
                     legend.key.size = unit(2.5, 'lines')) + 
      labs(x="Threshold value", y="Taxon count")
    
    plot2 <- ggplot2::ggplot(data = est.obj, aes(x = threshold.value, y = read.percent)) + ggplot2::geom_line(size=2, color="orangered1") + 
      ggplot2::theme(axis.text.y = element_text(size = 15, colour = "black"),
                     axis.text.x = element_text(size = 15, colour = "black", angle = 315, vjust = 0.7),
                     axis.title.x = element_text(size = 15, colour = "black"),
                     axis.title.y = element_text(size = 15, colour = "black")) +
      labs(x="Threshold value", y="Read percent")
    
    cowplot::plot_grid(plot1, plot2, labels = "AUTO")
    #begin ASV threshold plots
  } else if (!is.null(x) && !is.null(y)){
    
    #create df of ASV data and remove filtered taxa
    ASV.df <- est.obj$ASV.filtering.stats
    ASV.df.f <- ASV.df[complete.cases(ASV.df[,5]),] #remove filtered taxa
    
    #plots
    if (x == "P" && y =="RA"){
      #plot RA by prevalence
      asv.plot <- ggplot2::ggplot(data = ASV.df.f, aes(x = ASV.prevalence.percent, y = ASV.read.percent, color=Phylum)) + ggplot2::geom_point(size=3) + 
        ggplot2::facet_wrap(taxrank, scales = "fixed") + ggplot2::scale_y_log10() +
        ggplot2::theme(axis.text.x = element_text(size = 10, angle = 315, vjust = 0.2), 
                       legend.position = "right",
                       strip.text = element_text(size=11, color="white"),
                       strip.background = element_rect(fill = "black")) +
        ggplot2::geom_hline(yintercept = RAT, lty=2) +
        ggplot2::geom_vline(xintercept = PFT) +
        ggplot2::labs(y="ASV abundance (%)", x= "ASV prevalence (%)")
      
      return(asv.plot)
      
    } else if (x == "P" && y =="CV"){
      #plot CV by prevalence
      asv.plot <- ggplot2::ggplot(data = ASV.df.f, aes(x = ASV.prevalence.percent, y = ASV.CV, color=Phylum)) + ggplot2::geom_point(size=3) + 
        ggplot2::facet_wrap(taxrank, scales = "fixed") + 
        ggplot2::theme(axis.text.x = element_text(size = 10, angle = 315, vjust = 0.2), 
                       legend.position = "right",
                       strip.text = element_text(size=11, color="white"),
                       strip.background = element_rect(fill = "black")) +
        ggplot2::geom_hline(yintercept = CVT, lty=2) +
        ggplot2::geom_vline(xintercept = PFT) +
        ggplot2::labs(y="ASV CV", x= "ASV prevalence (%)")
      
      return(asv.plot)
      
    } else if (x == "RA" && y =="CV"){
      #plot CV by RA
      asv.plot <- ggplot2::ggplot(data = ASV.df.f, aes(x = ASV.read.percent, y = ASV.CV, color=Phylum)) + ggplot2::geom_point(size=3) + 
        ggplot2::facet_wrap(taxrank, scales = "fixed") + ggplot2::scale_x_log10() +
        ggplot2::theme(axis.text.x = element_text(size = 10, angle = 315, vjust = 0.2), 
                       legend.position = "right",
                       strip.text = element_text(size=11, color="white"),
                       strip.background = element_rect(fill = "black")) +
        ggplot2::geom_hline(yintercept = CVT, lty=2) +
        ggplot2::geom_vline(xintercept = RAT) +
        ggplot2::labs(y="ASV CV", x= "ASV abundance (%)")
      
      return(asv.plot)
      
    } else if (x == "RA" && y =="P"){
      #plot prevalence by RA
      asv.plot <- ggplot2::ggplot(data = ASV.df.f, aes(x = ASV.read.percent, y = ASV.prevalence.percent, color=Phylum)) + ggplot2::geom_point(size=3) + 
        ggplot2::facet_wrap(taxrank, scales = "fixed") + ggplot2::scale_x_log10() +
        ggplot2::theme(axis.text.x = element_text(size = 10, angle = 315, vjust = 0.2), 
                       legend.position = "right",
                       strip.text = element_text(size=11, color="white"),
                       strip.background = element_rect(fill = "black")) +
        ggplot2::geom_hline(yintercept = PFT, lty=2) +
        ggplot2::geom_vline(xintercept = RAT) +
        ggplot2::labs(y="ASV prevalence (%)", x= "ASV abundance (%)")
      
      return(asv.plot)
      
    } else if (x == "CV" && y =="RA"){
      #plot RA by CV
      asv.plot <- ggplot2::ggplot(data = ASV.df.f, aes(x = ASV.CV, y = ASV.read.percent, color=Phylum)) + ggplot2::geom_point(size=3) + 
        ggplot2::facet_wrap(taxrank, scales = "fixed") + ggplot2::scale_y_log10() +
        ggplot2::theme(axis.text.x = element_text(size = 10, angle = 315, vjust = 0.2), 
                       legend.position = "right",
                       strip.text = element_text(size=11, color="white"),
                       strip.background = element_rect(fill = "black")) +
        ggplot2::geom_hline(yintercept = RAT, lty=2) +
        ggplot2::geom_vline(xintercept = CVT) +
        ggplot2::labs(y="ASV abundance (%)", x= "ASV CV")
      
      return(asv.plot)
      
    } else if (x == "CV" && y =="P"){
      #plot prevalence by CV
      asv.plot <- ggplot2::ggplot(data = ASV.df.f, aes(x = ASV.CV, y = ASV.prevalence.percent, color=Phylum)) + ggplot2::geom_point(size=3) + 
        ggplot2::facet_wrap(taxrank, scales = "fixed") + 
        ggplot2::theme(axis.text.x = element_text(size = 10, angle = 315, vjust = 0.2), 
                       legend.position = "right",
                       strip.text = element_text(size=11, color="white"),
                       strip.background = element_rect(fill = "black")) +
        ggplot2::geom_hline(yintercept = PFT, lty=2) +
        ggplot2::geom_vline(xintercept = CVT) +
        ggplot2::labs(y="ASV prevalence (%)", x= "ASV CV")
      
      return(asv.plot)
      
    }
    #begin AS threshold plots
  } else {
    #RA
    if (!is.null(est.obj$relative.abundance.filtering.stats)){
      RA.stats <- est.obj$relative.abundance.filtering.stats
      plot3 <- ggplot2::ggplot(data = RA.stats, aes(x = relative.abundance.filter, y = ASV.count)) + ggplot2::geom_line(size=2, color="steelblue2") + 
        ggplot2::theme(axis.text.y = element_text(size = 15, colour = "black"),
                       axis.text.x = element_text(size = 15, colour = "black", angle = 315, vjust = 0.7),
                       axis.title.x = element_text(size = 15, colour = "black"),
                       axis.title.y = element_text(size = 15, colour = "black")) +
        labs(x="Relative abundance threshold", y="Taxon count")
    } else {
      RA.stats <- NULL
      plot3 <- NULL
    }
    #CV
    if (!is.null(est.obj$CV.filtering.stats)){
      CV.stats <- est.obj$CV.filtering.stats
      plot4 <- ggplot2::ggplot(data = CV.stats, aes(x = CV.filter, y = ASV.count)) + ggplot2::geom_line(size=2, color="orangered1") + 
        ggplot2::theme(axis.text.y = element_text(size = 15, colour = "black"),
                       axis.text.x = element_text(size = 15, colour = "black", angle = 315, vjust = 0.7),
                       axis.title.x = element_text(size = 15, colour = "black"),
                       axis.title.y = element_text(size = 15, colour = "black")) +
        labs(x="CV threshold", y="Taxon count")
    } else {
      CV.stats <- NULL
      plot4 <- NULL
    }
    #P
    if (!is.null(est.obj$prevalence.filtering.stats)){
      P.stats <- est.obj$prevalence.filtering.stats
      plot5 <- ggplot2::ggplot(data = P.stats, aes(x = prevalence.filter, y = ASV.count)) + ggplot2::geom_line(size=2, color="forestgreen") + 
        ggplot2::theme(axis.text.y = element_text(size = 15, colour = "black"),
                       axis.text.x = element_text(size = 15, colour = "black", angle = 315, vjust = 0.7),
                       axis.title.x = element_text(size = 15, colour = "black"),
                       axis.title.y = element_text(size = 15, colour = "black")) +
        labs(x="Prevalence threshold", y="Taxon count")
    } else {
      P.stats <- NULL
      plot5 <- NULL
    }
    
    plist <- list(plot3, plot4, plot5)
    cowplot::plot_grid(labels = "AUTO", plotlist = plist[which(!sapply(plist, is.null))])
    
  }
}

#PROCESSING FUNCTIONS
#standardization 
standardize.median <- function(ps){
  median.rc <- median(phyloseq::sample_sums(ps))
  ps.t <- phyloseq::transform_sample_counts(ps, fun = function(x) round(median.rc * (x/sum(x))))
  ps.t
}

#parameter estimation
getWS <- function(ps, WSrange, controlID, controlFASTA=NULL, useREFSEQ=FALSE){
  
  #Build param lists
  l.t <- seq(from = WSrange[1], to = WSrange[2], by = WSrange[3])
  nt <- length(l.t)
  tvec <- c()
  svec <- c()
  pvec <- c()
  
  #build files for sequence matching if specified
  if (!is.null(controlFASTA)){
    nvec <- list()
    mvec <- c()
    cfasta <- ShortRead::readFasta(controlFASTA)
  }
  
  
  for (i in 1:nt){
    tryCatch({
      #loop through values
      ps.ws <- suppressMessages(WSfilter(ps = ps, WST = l.t[i]))
      asv.tab <- format.ASV.tab(ps.ws)
      
      #FILTERING
      tvec[i] <- nrow(asv.tab[which(asv.tab[,match(controlID, colnames(asv.tab))] != 0),])
      svec[i] <- sum(phyloseq::sample_sums(ps.ws))
      pvec[i] <- sum(phyloseq::sample_sums(ps.ws))/sum(phyloseq::sample_sums(ps))*100
      
      if (!is.null(controlFASTA)){
        if (isTRUE(useREFSEQ)){
          #from phyloseq refseq slot
          control.taxanames <- rownames(asv.tab[which(asv.tab[,match(controlID, colnames(asv.tab))] != 0),])
          nvec[[i]] <- phyloseq::refseq(ps.ws)[control.taxanames]
        } else {
          nvec[[i]] <- rownames(asv.tab[which(asv.tab[,match(controlID, colnames(asv.tab))] != 0),])
        }
        #calculate matches to control FASTA
        mvec[i] <- sum(sapply(nvec[[i]], function(x) any(grepl(x, as.character(ShortRead::sread(cfasta))))))
      }
      
    },
    error=function(e){cat("Warning :",conditionMessage(e), "\n")})
  }
  names(tvec) <- c(l.t)
  names(svec) <- c(l.t)
  names(pvec) <- c(l.t)
  if (!is.null(controlFASTA)){
    names(mvec) <- c(l.t)
  } else {
    mvec <- rep(NA, length(tvec))
  }
  
  df <- as.data.frame(cbind(as.numeric(paste0(names(tvec))), tvec, mvec, svec, pvec))
  colnames(df) <- c("threshold.value", "control.taxa.count", "control.taxa.matches", "read.count", "read.percent")
  rownames(df) <- seq(1:length(l.t))
  df
}

getCV <- function(ps, WST=NULL, CVrange){
  
  #WS filtering
  if(is.null(WST)){
    message('Not applying WS filter')
    ps.ws <- ps
  } else {
    ps.ws <- WSfilter(ps = ps, WST = WST)
  }
  #standardize to median sample depth
  ps.wsm <- standardize.median(ps.ws)
  
  #build param vectors
  if(is.null(CVrange)){
    dfc <- NULL
  } else {
    l.c <- seq(from = CVrange[1], to = CVrange[2], by = CVrange[3])
    nc <- length(l.c)
    cvec <- c()
    
    for (i in 1:nc){
      tryCatch({
        #loop through values and filter
        ps.wsm <- phyloseq::filter_taxa(ps.wsm, function(x) sd(x)/mean(x) > l.c[i], TRUE)
        cvec[i] <- phyloseq::ntaxa(ps.wsm)
        
      },
      error=function(e){cat("Warning :c",conditionMessage(e), "\n")})
    }
    #create df
    #make both vectors same length and add NAs if CV filter zeroed out ASV table
    length(cvec) <- length(l.c)
    dfc <- as.data.frame(cbind(l.c, cvec))
    colnames(dfc) <- c("CV.filter", "ASV.count")
    rownames(dfc) <- seq(1:length(l.c))
  }
  
  dfc
}

getRA <- function(ps, WST=NULL, RArange){
  
  #WS filtering
  if(is.null(WST)){
    message('Not applying WS filter')
    ps.ws <- ps
  } else {
    ps.ws <- WSfilter(ps = ps, WST = WST)
  }
  
  #build param vectors
  if(is.null(RArange)){
    dfr <- NULL
  } else {
    l.r <- seq(from = RArange[1], to = RArange[2], by = RArange[3])
    nr <- length(l.r)
    rvec <- c()
    
    for (i in 1:nr){
      tryCatch({
        #loop through values and filter
        raf <- sum(phyloseq::taxa_sums(ps.ws)) * l.r[i]
        ps.ws <- phyloseq::prune_taxa(taxa_sums(ps.ws)>=raf, ps.ws)
        rvec[i] <- phyloseq::ntaxa(ps.ws)
        
      },
      error=function(e){cat("Warning :r",conditionMessage(e), "\n")})
    }
    #create df
    #make both vectors same length and add NAs if RF filter zeroed out ASV table
    length(rvec) <- length(l.r)
    dfr <- as.data.frame(cbind(l.r, rvec))
    colnames(dfr) <- c("relative.abundance.filter", "ASV.count")
    rownames(dfr) <- seq(1:length(l.r))
  }
  
  dfr
}

getPrev <- function(ps, WST=NULL, Prange){
  
  #WS filtering
  if(is.null(WST)){
    message('Not applying WS filter')
    ps.ws <- ps
  } else {
    ps.ws <- WSfilter(ps = ps, WST = WST)
  }
  
  asv.tab <- format.ASV.tab(ps.ws)
  
  #build param vectors
  l.p <- seq(from = Prange[1], to = Prange[2], by = Prange[3])
  np <- length(l.p)
  pvec <- c()
  taxa.cvec <- c()
  
  #get prevalence list for each taxon
  taxa.plist <- apply(X = asv.tab, MARGIN = 1, FUN = function(x){names(x)[which(x!=0)]})
  #populate vector of sample counts per ASV
  for(j in 1:length(taxa.plist)){
    taxa.cvec[j] <- length(taxa.plist[[j]])
  }
  
  #loop through params
  for (i in 1:np){
    tryCatch({
      #apply ASV names and filter ASVs below PF
      names(taxa.cvec) <- names(taxa.plist)
      prev.count <- phyloseq::nsamples(ps.ws)* l.p[i]
      taxa.cvec.f <- taxa.cvec[which(taxa.cvec > prev.count)]
      tn.cvec.f <- names(taxa.cvec.f)
      #filter ps
      ps.ws <- phyloseq::prune_taxa(tn.cvec.f, ps.ws)
      pvec[i] <- phyloseq::ntaxa(ps.ws)
      
    },
    error=function(e){cat("Warning :p",conditionMessage(e), "\n")})
  }
  
  #create df
  #make both vectors same length and add NAs if PF filter zeroed out ASV table
  length(pvec) <- length(l.p)
  dfp <- cbind.data.frame(l.p, pvec)
  colnames(dfp) <- c("prevalence.filter", "ASV.count")
  rownames(dfp) <- seq(1:length(l.p))
  #name taxa prevalence vector
  names(taxa.cvec) <- names(taxa.plist)
  
  # Build return list
  l.return = list()
  l.return[['prevalence.filtering.stats']] <- dfp
  l.return[['ASV.prevalence.count']] <- taxa.cvec
  
  return(l.return)
  
}

#filtering scripts
WSfilter <- function(ps, WST){
  
  #perform filter
  message('Applying WS filter threshold of ', WST)
  
  filterfx = function(x){
    x[(x / sum(x)) < WST] <- 0
    return(x)
  }
  
  ps <- phyloseq::transform_sample_counts(ps, fun = filterfx)
  ps
}

MDfilter <- function(ps, mdFACTOR, mdCAT, mdNEGATIVE=FALSE){
  #create sample df for subsetting
  sampledf <- suppressWarnings(as.matrix(phyloseq::sample_data(ps)))
  
  if (isTRUE(mdNEGATIVE)){
    filtered.names <- rownames(sampledf[which(sampledf[,match(mdCAT, colnames(sampledf))] == mdFACTOR),])
    message('Removing ',  (phyloseq::nsamples(ps) - length(filtered.names)), ' samples not matching metadata identifiers ', mdCAT, ":", mdFACTOR)
  } else {
    filtered.names <- rownames(sampledf[which(sampledf[,match(mdCAT, colnames(sampledf))] != mdFACTOR),])
    message('Removing ',  (phyloseq::nsamples(ps) - length(filtered.names)), ' samples matching metadata identifiers ', mdCAT, ":", mdFACTOR)
  }    
  
  #subset sampledf to include nonfiltered samples only
  sampledf.s <- as.data.frame(sampledf[filtered.names,])
  phyloseq::sample_data(ps) <- phyloseq::sample_data(sampledf.s)
  ps
}

CVfilter <- function(ps, WST=NULL, CVF){
  
  #WS filtering
  if(is.null(WST)){
    message('Not applying WS filter')
    ps.ws <- ps
  } else {
    ps.ws <- WSfilter(ps = ps, WST = WST)
  }
  
  #standardize to median sample depth
  ps.wsm <- standardize.median(ps.ws)
  
  #perform filter
  ps.wsm <- phyloseq::filter_taxa(ps.wsm, function(x) sd(x)/mean(x) > CVF, TRUE)
  #get taxa names to apply to original, unstandardized dataset
  filtered.taxa.names <- phyloseq::taxa_names(ps.wsm)
  #apply filter to unstandardized dataset
  ps.ws <- phyloseq::prune_taxa(taxa = filtered.taxa.names, x = ps.ws)
  ps.ws
  
}


RAfilter<- function(ps, WST=NULL, RAF){
  
  #WS filtering
  if(is.null(WST)){
    message('Not applying WS filter')
    ps.ws <- ps
  } else {
    ps.ws <- WSfilter(ps = ps, WST = WST)
  }
  
  #perform filter
  raf <- sum(phyloseq::taxa_sums(ps.ws)) * RAF
  ps.ws <- phyloseq::prune_taxa(taxa_sums(ps.ws)>=raf, ps.ws)
  ps.ws
  
}

Pfilter <- function(ps, WST=NULL, PF){
  
  #WS filtering
  if(is.null(WST)){
    message('Not applying WS filter')
    ps.ws <- ps
  } else {
    ps.ws <- WSfilter(ps = ps, WST = WST)
  }
  
  #format asv table
  asv.tab <- format.ASV.tab(ps.ws)
  
  #create sample count vector
  taxa.cvec <- c()
  #get prevalence list for each taxon
  taxa.plist <- apply(X = asv.tab, MARGIN = 1, FUN = function(x){names(x)[which(x!=0)]})
  #populate vector of sample counts per ASV
  for(j in 1:length(taxa.plist)){
    taxa.cvec[j] <- length(taxa.plist[[j]])
  }
  names(taxa.cvec) <- names(taxa.plist)
  prev.count <- phyloseq::nsamples(ps.ws)* PF
  taxa.cvec.f <- taxa.cvec[which(taxa.cvec > prev.count)]
  tn.cvec.f <- names(taxa.cvec.f)
  #perform filter
  ps.ws <- phyloseq::prune_taxa(tn.cvec.f, ps.ws)
  ps.ws
  
}

#WRAPPER FUNCTIONS
estimate.WSthreshold <- function(ps, WSrange, controlID, controlFASTA=NULL, useREFSEQ=FALSE) {
  
  #throw error if controlID doesn't match
  if(!(controlID %in% phyloseq::sample_names(ps))){
    stop("controlID provided is not a valid sample name")
  }
  
  #convert param string to numeric vector
  string.w <- substitute(WSrange)
  WST <- eval(expr = format.parameter.string(string = string.w), envir = parent.frame())
  message('Estimating filtering statistics from WS thresholds ', WST[1], ' to ', WST[2], ' by ', WST[3])
  gws <- getWS(ps = ps, WSrange = WST, controlID = controlID, controlFASTA = controlFASTA, useREFSEQ=useREFSEQ)
  gws
  
}



estimate.ASthreshold <- function(ps, WST=NULL, RAT=NULL, CVT=NULL, PFT=NULL, mdCAT=NULL, mdFACTOR=NULL, mdNEGATIVE=FALSE,
                                 minLIB=NULL, Prange=NULL, CVrange=NULL, RArange=NULL){
  
  #throw error if mdCAT doesn't match
  if(all(!is.null(mdCAT), !(mdCAT %in% colnames(phyloseq::sample_data(ps))))){
    stop("mdCAT provided is not a valid metadata category")
  }
  
  #remove samples < minlib
  if(is.null(minLIB)){
    ps = ps
  } else {
    pml.c <- nrow(phyloseq::sample_data(ps))
    ps = phyloseq::prune_samples(phyloseq::sample_sums(ps)>=minLIB, ps)
    message('Removing ',(pml.c -  phyloseq::nsamples(ps)), ' samples with read count < ', minLIB)
  }
  
  #WS filtering
  if (is.null(WST)){
    ps.ws <- ps
  } else {
    ps.ws <- WSfilter(ps = ps, WST = WST)
  }
  
  #save WS filtered object for reversion later
  ps.wso <- ps.ws
  
  #METADATA BASED SAMPLE FILTERING
  if (any(c(is.null(mdCAT), is.null(mdFACTOR)))){
    ps.ws <- ps.ws
  } else {
    ps.ws <- MDfilter(ps = ps.ws, mdFACTOR = mdFACTOR, mdCAT = mdCAT, mdNEGATIVE = mdNEGATIVE)
  }
  
  #INCORPORATE FIXED THRESHOLDS
  #RELATIVE ABUNDANCE
  if(!is.null(RAT)){
    ps.ws <- suppressMessages(RAfilter(ps = ps.ws, WST = NULL, RAF = RAT))
    message('Applying fixed relative abundance threshold of ', RAT)
  }
  #CV
  if(!is.null(CVT)){
    ps.ws <- suppressMessages(CVfilter(ps = ps.ws, WST = NULL, CVF = CVT))
    message('Applying fixed CV threshold of ', CVT)
  }
  #PREVALENCE
  if(!is.null(PFT)){
    ps.ws <- suppressMessages(Pfilter(ps = ps.ws, WST = NULL, PF = PFT))
    message('Applying fixed prevalence threshold of ', PFT)
  }
  
  #ESTIMATION
  #RELATIVE ABUNDANCE
  #build param lists
  if(is.null(RArange)){
    gr <- NULL
  } else {
    #convert param string to numeric vector
    string.r <- substitute(RArange)
    RAF <- eval(expr = format.parameter.string(string = string.r), envir = parent.frame())
    message('Estimating filtering statistics from relative abundance thresholds ', RAF[1], ' to ', RAF[2], ' by ', RAF[3])
    gr <- suppressMessages(getRA(ps = ps.ws, WST = NULL, RArange = RAF))
    
  }
  
  #CV
  #build param lists
  if(is.null(CVrange)){
    gc <- NULL
  } else {
    string.c <- substitute(CVrange)
    CVF <- eval(expr = format.parameter.string(string = string.c), envir = parent.frame())
    message('Estimating filtering statistics from CV thresholds ', CVF[1], ' to ', CVF[2], ' by ', CVF[3])
    gc <- suppressMessages(getCV(ps = ps.ws, WST = NULL, CVrange = CVF))
  }
  
  #PREVALENCE
  #Build param lists
  if(is.null(Prange)){
    gp <- NULL
  } else {
    string.p <- substitute(Prange)
    PF <- eval(expr = format.parameter.string(string = string.p), envir = parent.frame())
    message('Estimating filtering statistics from prevalence thresholds ', PF[1], ' to ', PF[2], ' by ', PF[3])
    gp <- suppressMessages(getPrev(ps = ps.ws, WST = NULL, Prange = PF))
  }
  
  #CREATE ASV DF
  #build df vectors
  ts <- taxa_sums(ps.wso)
  tsp <- taxa_sums(ps.wso)/sum(taxa_sums(ps.wso)) * 100
  namevec <- names(ts)
  #standardize to median sample depth for CV calculation
  ps.ws <- standardize.median(ps.wso)
  asv.tab <- format.ASV.tab(ps.ws)
  cv.asv <- apply(asv.tab[namevec,], MARGIN = 1, FUN = function(x) sd(x)/mean(x))
  tax.tab <- phyloseq::tax_table(ps.wso)[namevec,]
  
  #set prev vectors to null if no prev stats desired
  if (all(c(is.null(Prange), !is.null(PFT)))){
    gp.reload <- suppressMessages(getPrev(ps = ps.ws, WST = NULL, Prange = c(0.10,0.11,0.01)))
    taxa.cvec <- gp.reload$ASV.prevalence.count
    prev <- taxa.cvec[namevec]
    prevp <- prev/phyloseq::nsamples(ps.ws) * 100
  } else if (all(c(is.null(Prange), is.null(PFT)))){
    prev <- rep(NA, length(ts))
    prevp <- rep(NA, length(ts))
  } else {
    gp.reload <- suppressMessages(getPrev(ps = ps.ws, WST = NULL, Prange = PF))
    taxa.cvec <- gp.reload$ASV.prevalence.count
    prev <- taxa.cvec[namevec]
    prevp <- prev/phyloseq::nsamples(ps.ws) * 100
  }
  
  #build df and rename
  df.asv <- cbind.data.frame(ts, tsp, prev, prevp, cv.asv, tax.tab, rownames(tax.tab))
  colnames(df.asv)[1:5] <- c("ASV.read.count", "ASV.read.percent", "ASV.prevalence", "ASV.prevalence.percent", "ASV.CV")
  colnames(df.asv)[ncol(df.asv)] <- "ASV.ID"
  rownames(df.asv) <- seq(1:nrow(df.asv))
  
  # Build return list
  l.return = list()
  l.return[['relative.abundance.filtering.stats']] <- gr
  l.return[['CV.filtering.stats']] <- gc
  l.return[['prevalence.filtering.stats']] <- gp$prevalence.filtering.stats
  l.return[['ASV.filtering.stats']] <- df.asv
  
  return(l.return)
}

microfilter <- function(ps, controlID=NULL, mdCAT=NULL, mdFACTOR=NULL, mdNEGATIVE=FALSE, minLIB=NULL, WST=NULL, RAT=NULL, CVT=NULL, PFT=NULL, return.all=FALSE){
  
  #throw error if controlID doesn't match
  if(all(!is.null(controlID), !(controlID %in% phyloseq::sample_names(ps)))){
    stop("controlID provided is not a valid sample name")
  }
  
  #throw error if mdCAT doesn't match
  if(all(!is.null(mdCAT), !(mdCAT %in% colnames(phyloseq::sample_data(ps))))){
    stop("mdCAT provided is not a valid metadata category")
  }
  
  #remove samples < minlib
  if(is.null(minLIB)){
    ps = ps
  } else {
    pml.c <- nrow(phyloseq::sample_data(ps))
    ps = phyloseq::prune_samples(phyloseq::sample_sums(ps)>=minLIB, ps)
    message('Removing ',(pml.c -  phyloseq::nsamples(ps)), ' samples with read count < ', minLIB)
  }
  
  #create unfiltered sample sum vector
  ov <- phyloseq::sample_sums(ps)
  
  #WS filtering
  if (is.null(WST)){
    ps.ws <- ps
  } else {
    ps.ws <- WSfilter(ps = ps, WST = WST)
  }
  
  #create WS filtered sample sum vector
  ifv <- phyloseq::sample_sums(ps.ws)
  #calculate percent filtered, individual
  p.if <- phyloseq::sample_sums(ps.ws)/phyloseq::sample_sums(ps)*100
  
  asv.tab <- format.ASV.tab(ps.ws)
  
  #METADATA-BASED SAMPLE REMOVAL
  if(is.null(controlID)){
    npos <- NULL
    tax.tab.subset <- NULL
    ttsn <- NULL
  } else {
    #calculate control taxa count
    npos <- nrow(asv.tab[which(asv.tab[,match(controlID, colnames(asv.tab))] != 0),])
    
    #get taxonomy of taxa in positive control
    tax.tab <- phyloseq::tax_table(ps.ws)
    taxanames.control <- rownames(asv.tab[which(asv.tab[,match(controlID, colnames(asv.tab))] != 0),])
    tax.tab.subset <- tax.tab[taxanames.control] #taxonomy of taxa in positive control
    ttsn <- tax.tab.subset
    rownames(ttsn) <- NULL
  }
  
  #remove samples by metadata filters
  if (any(c(is.null(mdCAT), is.null(mdFACTOR)))){
    ps.ws <- ps.ws
  } else {
    ps.ws <- MDfilter(ps = ps.ws, mdFACTOR = mdFACTOR, mdCAT = mdCAT, mdNEGATIVE = mdNEGATIVE)
  }
  
  #AS filtering
  #relative abundance filter
  if(is.null(RAT)){
    ps.ws <- ps.ws
    raf <- NULL
  } else {
    message('Applying relative abundance threshold of ', RAT)
    ps.ws <- suppressMessages(RAfilter(ps = ps.ws, WST = NULL, RAF = RAT))
    raf <- RAT * sum(phyloseq::taxa_sums(ps.ws))
  }
  
  #CV filter
  if(is.null(CVT)){
    ps.ws <- ps.ws
  } else {
    message('Applying CV threshold of ', CVT)
    ps.ws <- suppressMessages(CVfilter(ps = ps.ws, WST = NULL, CVF = CVT))
  }
  
  #prevalence filter
  if(is.null(PFT)){
    ps.ws <- ps.ws
    prev.count <- NULL
  } else {
    message('Applying prevalence threshold of ', PFT)
    ps.ws <- suppressMessages(Pfilter(ps = ps.ws, WST = NULL, PF = PFT))
    prev.count <- phyloseq::nsamples(ps.ws) * PFT
  }
  
  #create AS filter sample sum vector
  pfv <- phyloseq::sample_sums(ps.ws)
  #calculate percent filtered, AS
  p.pf <- suppressWarnings(phyloseq::sample_sums(ps.ws)/phyloseq::sample_sums(ps)[names(phyloseq::sample_sums(ps.ws))]*100)
  
  #order vectors
  pfv <- pfv[names(p.if)]
  p.pf <- p.pf[names(p.if)]
  
  #cbind vectors into df
  sstab <- cbind(ov, ifv, p.if, pfv, p.pf)
  colnames(sstab) <- c("unfiltered.read.count", "WSfiltered.read.count", "WSfiltered.read.percent", "ASfiltered.read.count", "ASfiltered.read.percent")
  
  # Build return list
  l.return = list()
  if (return.all==FALSE){
    return(ps.ws)
  } else {
    l.return[['filtered.phyloseq']] <- ps.ws
    l.return[['ntaxa.in.control']] <- npos
    l.return[['control.taxa.sequences']] <- rownames(tax.tab.subset)
    l.return[['taxonomy.of.control.taxa']] <- ttsn
    l.return[['read.count.table']] <- sstab
    l.return[['relative.abundance.filter.read.count']] <- raf
    l.return[['prevalence.filter.sample.count']] <- prev.count
    
  }
  
  return(l.return)
}

write.dataset.biom <- function(ps, filePATH, filePREFIX, writeFASTA=TRUE, rename=FALSE, useREFSEQ=FALSE){
  
  #pull seqs from refseq slot or extract from ASV ID for fasta format
  if (isTRUE(useREFSEQ)){
    #from phyloseq refseq slot
    f.onames <- phyloseq::refseq(ps)
  } else {
    f.onames <- phyloseq::taxa_names(ps)
  }
  
  if (isTRUE(rename)){
    phyloseq::taxa_names(ps) <- paste("ASV", 1:length(phyloseq::taxa_names(ps)), sep = "")
    names(f.onames) <- paste0(">", phyloseq::taxa_names(ps))
  } else {
    names(f.onames) <- paste0(">", phyloseq::taxa_names(ps))
  }
  
  #generate biom file
  suppressWarnings(ps.b <- biomformat::make_biom(
    data = format.ASV.tab(ps),
    sample_metadata = as.data.frame(phyloseq::sample_data(ps)),
    observation_metadata = as.data.frame(phyloseq::tax_table(ps)), 
    matrix_element_type = "int"
  )
  )
  
  #create output string
  if (isTRUE(writeFASTA)){
    fa <- print(paste0(filePATH, filePREFIX, "_ASVs.fasta"))
  }
  bo <- print(paste0(filePATH, filePREFIX, "_ASV_table.biom"))
  
  #write output
  if (isTRUE(writeFASTA)){
    write.table(x = f.onames, file = fa, quote = FALSE, sep = "\n", col.names = FALSE)
  }
  #biom export
  biomformat::write_biom(x = ps.b, biom_file = bo)
  
  #return phyloseq object with taxa renamed to ASV1, etc., if desired
  if (isTRUE(rename)){
    return(ps)
  }
}

write.dataset <- function(ps, filePATH, filePREFIX, writeFASTA=TRUE, rename=FALSE, useREFSEQ=FALSE){
  
  #pull seqs from refseq slot or extract from ASV ID for fasta format
  if (isTRUE(useREFSEQ)){
    #from phyloseq refseq slot
    f.onames <- phyloseq::refseq(ps)
  } else {
    f.onames <- phyloseq::taxa_names(ps)
  }
  
  if (isTRUE(rename)){
    phyloseq::taxa_names(ps) <- paste("ASV", 1:length(phyloseq::taxa_names(ps)), sep = "")
    names(f.onames) <- paste0(">", phyloseq::taxa_names(ps))
  } else {
    names(f.onames) <- paste0(">", phyloseq::taxa_names(ps))
  }
  
  
  #generate asv table formatted for biom generation
  asv.tab <- format.ASV.tab(ps)
  suppressWarnings(asv.tab <- as.matrix(asv.tab))
  cb <- as.matrix(cbind(rownames(asv.tab), asv.tab))
  rcb <- as.matrix(rbind(colnames(cb), cb))
  rcb[1,1] <- "#ASVID"
  rownames(rcb) <- NULL
  colnames(rcb) <- NULL
  
  #generate tax table formatted for biom generation
  tax.tab <- as.data.frame(phyloseq::tax_table(ps))
  tax.tab$taxonomy <- tidyr::unite_(tax.tab, "out", c(colnames(tax.tab)), sep = ";")
  cbt <- as.matrix(cbind(rownames(tax.tab), tax.tab$taxonomy))
  rcbt <- as.matrix(rbind(c("#ASVID", "taxonomy"), cbt))
  rownames(cbt) <- NULL
  colnames(cbt) <- NULL
  
  #generate sampledf table formatted for biom generation
  samdf <- suppressWarnings(as.matrix(phyloseq::sample_data(ps)))
  cbs <- as.matrix(cbind(rownames(samdf), samdf))
  rcbs <- as.matrix(rbind(colnames(cbs), cbs))
  rcbs[1,1] <- "#SampleID"
  rownames(rcbs) <- NULL
  colnames(rcbs) <- NULL
  
  #create output string
  if (isTRUE(writeFASTA)){
    fa <- print(paste0(filePATH, filePREFIX, "_ASVs.fasta"))
  }
  otb <- print(paste0(filePATH, filePREFIX, "_ASV_table.txt"))
  ttb <- print(paste0(filePATH, filePREFIX, "_ASV_taxonomy.txt"))
  stb <- print(paste0(filePATH, filePREFIX, "_sample_data.txt"))
  
  
  #write output
  #ASV fasta 
  if (isTRUE(writeFASTA)){
    write.table(x = f.onames, file = fa, quote = FALSE, sep = "\n", col.names = FALSE)
  }
  #asv.tab
  write.table(x = rcb, file = otb, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  #tax.tab
  write.table(x = rcbt, file = ttb, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  #sampledf
  write.table(x = rcbs, file = stb, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  #return phyloseq object with taxa renamed to ASV1, etc., if desired
  if (isTRUE(rename)){
    return(ps)
  }
}


