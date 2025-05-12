###########################################
# Figure 2 - Ileal Disease - Genera level #
###########################################
devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)
res <- CD_tissue_padj_fig_data

res$Genus <- as.character(res$Genus)
res$GS <- paste(res$Genus,res$Species,sep=" ")
#Microbiota Palette for Reference
microbiota.palette <- c("Acidobacteria"="#be5400","Actinobacteria"="#d92728",
                        "Armatimonadetes"="#d36ff8","Bacteroidetes"="#0955da",
                        "Candidatus_Saccharibacteria"="#88374d","Chloroflexi"="#99c618",
                        "Cyanobacteria/Chloroplast"="#6adab2","Deferribacteres"= "#93aa00",
                        "Deinococcus-Thermus"="#00984c","Firmicutes"="#3b8b00",
                        "Fusobacteria"="#ff5fa3","Gemmatimonadetes"="#63b6ff",
                        "Lentisphaerae"="#ff8f54","Parcubacteria"="#006045",
                        "Proteobacteria"="#72349e","Spirochaetes"="#fbaee3",
                        "SR1"="#50570e","Synergistetes"="#fabb5e",
                        "Tenericutes"="#7b4601","Verrucomicrobia"="#01cacc")
#Setup TaxaVals
TaxaVals <- rep('black',nrow(res))
names(TaxaVals) <- rep('Mid',nrow(res))

TaxaVals[which(res$Phylum == "Actinobacteria")] <- "#d92728"
names(TaxaVals)[which(res$Phylum == "Actinobacteria")] <- "Actinobacteria"

TaxaVals[which(res$Phylum == "Bacteroidetes")] <- "#0955da"
names(TaxaVals)[which(res$Phylum == "Bacteroidetes")] <- "Bacteroidetes"

TaxaVals[which(res$Phylum == "Candidatus_Saccharibacteria")] <- "#88374d"
names(TaxaVals)[which(res$Phylum == "Candidatus_Saccharibacteria")] <- "Candidatus_Saccharibacteria"

TaxaVals[which(res$Phylum == "Firmicutes")] <- "#3b8b00"
names(TaxaVals)[which(res$Phylum == "Firmicutes")] <- "Firmicutes"

TaxaVals[which(res$Phylum == "Fusobacteria")] <- "#ff5fa3"
names(TaxaVals)[which(res$Phylum == "Fusobacteria")] <- "Fusobacteria"

TaxaVals[which(res$Phylum == "Proteobacteria")] <- "#72349e"
names(TaxaVals)[which(res$Phylum == "Proteobacteria")] <- "Proteobacteria"

TaxaVals[which(res$Phylum == "Verrucomicrobia")] <- "#01cacc"
names(TaxaVals)[which(res$Phylum == "Verrucomicrobia")] <- "Verrucomicrobia"

unique(names(TaxaVals))
TaxaVals[1:20]

#Setup ShapeVal 
#Obligate Anaerobes
OAn.CD.taxa <- c("Bifidobacterium",
                 "Alistipes",
                 "Anaerovorax",
                 "Anaerostipes",
                 "Blautia",
                 "Clostridium_XlVa",
                 "Coprococcus",
                 "Fusicatenibacter",
                 "Genus_1595",
                 "Genus_27047",
                 "Genus_27243",
                 "Genus_3926",
                 "Roseburia",
                 "Ruminococcus2",
                 "Intestinibacter",
                 "Anaerotruncus",
                 "Clostridium_IV",
                 "Faecalibacterium",
                 "Flavonifractor",
                 "Intestinimonas",
                 "Oscillibacter",
                 "Pseudoflavonifractor",
                 "Ruminococcus",
                 "Subdoligranulum",
                 "Akkermansia",
                 "Bacteroides",
                 "Flavobacterium")

OA.CD.taxa <- c("Sphingobium",
                "Sphingomonas",
                "Acinetobacter",
                "Stenotrophomonas")

FAn.CD.taxa<- c("Streptococcus",
                "Rothia",
                "Propionibacterium",
                "Cloacibacterium",
                "Abiotrophia",
                "Granulicatella",
                "Enterococcus",
                "Streptococcus",
                "Haemophilus")

Aero.CD.taxa <- c("Bradyrhizobium",
                  "Acidovorax",
                  "Pseudomonas")

Anaero.CD.taxa <- c("Odoribacter",
                    "Veillonella")

TaxaVals.CD.shape <- rep(0,nrow(res))
names(TaxaVals.CD.shape) <- rep('NS',nrow(res))

TaxaVals.CD.shape[which(res$Genus %in% OAn.CD.taxa)] <- 16
names(TaxaVals.CD.shape)[which(res$Genus %in% OAn.CD.taxa)] <- 'Obligate Anaerobe'

TaxaVals.CD.shape[which(res$Genus %in% OA.CD.taxa)] <- 2
names(TaxaVals.CD.shape)[which(res$Genus %in% OA.CD.taxa)] <- 'Obligate Aerobe'

TaxaVals.CD.shape[which(res$Genus %in% FAn.CD.taxa)] <- 18
names(TaxaVals.CD.shape)[which(res$Genus %in% FAn.CD.taxa)] <- 'Facultative Anaerobe'

TaxaVals.CD.shape[which(res$Genus %in% Aero.CD.taxa)] <- 4
names(TaxaVals.CD.shape)[which(res$Genus %in% Aero.CD.taxa)] <- 'Aerobe'

TaxaVals.CD.shape[which(res$Genus %in% Anaero.CD.taxa)] <- 15
names(TaxaVals.CD.shape)[which(res$Genus %in% Anaero.CD.taxa)] <- 'Anaerobe'

unique(names(TaxaVals.CD.shape))
unique(TaxaVals.CD.shape)
TaxaVals.CD.shape[1:20]



#Rename
res$GS[res$GS == "Genus_27243 NA"] <- "Unknown Lachnospiraceae sp,"
res$GS[res$GS == "Genus_27047 NA"] <- "Unknown Lachnospiraceae sp."
res$GS[res$GS == "Genus_3926 NA"] <- "Unknown Lachnospiraceae sp."
res$GS[res$GS == "Genus_10545 NA"] <- "Unknown Erysipelotrichaceae sp."
res$GS[res$GS == "Genus_1595 NA"] <- "Unknown Lachnospiraceae sp."
res$GS[res$GS == "Genus_15493 NA"] <- "Unknown Neisseriaceae sp,"

res$GS[res$GS == "Blautia Species_1161"] <- "Blautia sp."
res$GS[res$GS == "Blautia Species_1174"] <- "Blautia sp."
res$GS[res$GS == "Blautia Species_715"] <- "Blautia sp."
res$GS[res$GS == "Coprococcus Species_3234"] <- "Coprococcus sp."
res$GS[res$GS == "Coprococcus Species_3674"] <- "Coprococcus sp."

res$GS[res$GS == "Clostridium_XlVa Species_1979"] <- "Clostridium_XlVa sp."
res$GS[res$GS == "Clostridium_XlVa Species_2209"] <- "Clostridium_XlVa sp."
res$GS[res$GS == "Clostridium_IV Species_7174"] <- "Clostridium_IV sp."
res$GS[res$GS == "Roseburia Species_26774"] <- "Roseburia sp."
res$GS[res$GS == "Faecalibacterium Species_4902"] <- "Faecalibacterium sp."
res$GS[res$GS == "Faecalibacterium Species_4686"] <- "Faecalibacterium sp."
res$GS[res$GS == "Faecalibacterium Species_4970"] <- "Faecalibacterium sp."
res$GS[res$GS == "Akkermansia NA"] <- "Akkermansia"
res$GS[res$GS == "Ruminococcus NA"] <- "Ruminococcus"
res$GS[res$GS == "Intestinimonas NA"] <- "Intestinimonas"
res$GS[res$GS == "Subdoligranulum NA"] <- "Subdoligranulum"
res$GS[res$GS == "Anaerovorax NA"] <- "Anaerovorax"
res$GS[res$GS == "Oscillibacter NA"] <- "Oscillibacter"
res$GS[res$GS == "Odoribacter NA"] <- "Odoribacter"
res$GS[res$GS == "Anaerotruncus NA"] <- "Anaerotruncus"
res$GS[res$GS == "Anaerostipes NA"] <- "Anaerostipes"
res$GS[res$GS == "Sphingomonas NA"] <- "Sphingomonas"
res$GS[res$GS == "Veillonella NA"] <- "Veillonella"
res$GS[res$GS == "Acidovorax NA"] <- "Acidovorax"
res$GS[res$GS == "Propionibacterium NA"] <- "Propionibacterium"
res$GS[res$GS == "Flavobacterium NA"] <- "Flavobacterium"
res$GS[res$GS == "Acinetobacter NA"] <- "Acinetobacter"
res$GS[res$GS == "Cloacibacterium NA"] <- "Cloacibacterium"
res$GS[res$GS == "Abiotrophia NA"] <- "Abiotrophia"

res$GS[res$GS == "Fusicatenibacter Species_384"] <- "Fusicatenibacter sp."
res$GS[res$GS == "Intestinibacter Species_23797"] <- "Intestinibacter sp."
res$GS[res$GS == "Pseudoflavonifractor Species_6225"] <- "Pseudoflavonifractor sp."
res$GS[res$GS == "Streptococcus Species_9137"] <- "Streptococcus sp."
res$GS[res$GS == "Pseudomonas Species_14795"] <- "Pseudomonas sp."

sigTaxa.CD<- c("Blautia sp.",
               "Coprococcus sp.",
               "Unknown Erysipelotrichaceae species",
               "Clostridium_XlVa sp.",
               "Faecalibacterium sp.",
               "Alistipes finegoldii",
               "Roseburia sp",
               "Bifidobacterium bifidum",
               "Faecalibacterium sp.",
               "Fusicatenibacter sp.",
               "Akkermansia sp.",
               "Unknown Lachnospiraceae species",
               "Unknown Lachnospiraceae species",
               "Ruminococcus sp.",
               "Faecalibacterium sp.",
               "Intestinimonas sp.",
               "Subdoligranulum sp.",
               "Anaerovorax sp.",
               "Oscillibacter sp.",
               "Intestinibacter sp.",
               "Odoribacter sp.",
               "Anaerotruncus sp.",
               "Anaerostipes sp.",
               "Blautia sp.",
               "Pseudoflavonifractor sp.",
               "Clostridium_IV sp.",
               "Blautia hydrogenotrophica",
               "Blautia sp.",
               "Coprococcus sp.",
               "Sphingomonas sp.",
               "Haemophilus parainfluenzae",
               "Rothia dentocariosa",
               "Streptococcus sp.",
               "Bifidobacterium pseudocatenulatum",
               "Coprococcus sp.",
               "Veillonella sp.",
               "Acidovorax sp.",
               "Bradyrhizobium group",
               "Unknown Neisseriaceae species",
               "Propionibacterium sp.",
               "Flavobacterium sp.",
               "Acinetobacter sp.",
               "Bacteroides ovatus",
               "Cloacibacterium sp.",
               "Pseudomonas Species_14795",
               "Granulicatella elegans",
               "Abiotrophia sp.",
               "Granulicatella adiacens",
               "Clostridium_XlVa sp.")

sigTaxa.CD.2 <- c("Blautia sp.",
               "Unknown Lachnospiraceae species",
               "Subdoligranulum",
               "Haemophilus parainfluenzae",
               "Rothia dentocariosa",
               "Streptococcus sp.",
               "Bacteroides ovatus",
               "Cloacibacterium sp.",
               "Pseudomonas sp.",
               "Clostridium_XlVa sp.")

res$log2FoldChange <- (-1)*res$log2FoldChange
CD.da <- EnhancedVolcano(res,
                             lab = res$GS,
                             x = 'log2FoldChange',
                             y = 'global.pvalue',
                             selectLab = sigTaxa.CD.2,
                             xlim = c(-6.5,6.5),
                             title = 'Alterations in the CD Mucosa',
                             ylab = bquote(~-Log[10]~global~italic(P)),
                             pCutoff = 0.045133334,
                             transcriptPointSize = res$log_grp_mra,
                             transcriptLabSize = 3.0,
                             transcriptLabCol = 'black',
                             transcriptLabFace = 'italic',
                             cutoffLineType = "twodash",
                             cutoffLineCol = "red4",
                             cutoffLineWidth = 1.0,
                             colAlpha = 1,
                             colCustom = TaxaVals,
                             shapeCustom = TaxaVals.CD.shape,
                             boxedlabels = TRUE,
                             drawConnectors = TRUE,
                             widthConnectors = 0.2,
                             colConnectors = 'black',
                             legendPosition = 'right',
                             legendLabSize = 16,
                             legendIconSize = 5.0,
                             gridlines.major = TRUE,
                             gridlines.minor = TRUE,
                             border = 'partial',
                             borderWidth = 1.5,
                             borderColour = 'black')

CD.da + coord_flip() 
