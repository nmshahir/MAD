###########################################
# Figure 2 - Ileal Disease - Genera level #
###########################################
devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)
res <- colon_disease_merge_fig_data_sp_gen_aero

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
OAn.colon.taxa <- c("Bacteroides",
                    "Blautia",
                    "Clostridium_XlVa",
                    "Coprococcus",
                    "Dialister",
                    "Genus_27047",
                    "Ruminococcus2",
                    "Alloprevotella",
                    "Barnesiella",
                    "Clostridium_IV",
                    "Clostridium_XlVb",
                    "Clostridium_XVIII",
                    "Dorea",
                    "Faecalibacterium",
                    "Fusicatenibacter",
                    "Genus_23841",
                    "Genus_25334",
                    "Genus_25335",
                    "Genus_25348",
                    "Genus_27241",
                    "Genus_6584",
                    "Mogibacterium",
                    "Oscillibacter",
                    "Romboutsia",
                    "Alistipes",
                    "Flavonifractor",
                    "Genus_26476",
                    "Genus_27153",
                    "Genus_27177",
                    "Hathewaya",
                    "Intestinibacter",
                    "Phascolarctobacterium",
                    "Ruminococcus")

OA.colon.taxa <- c("Janthinobacterium")

FAn.colon.taxa<- c("Exiguobacterium",
                   "Klebsiella",
                   "Actinomyces",
                   "Citrobacter",
                   "Genus_16629")

Aero.colon.taxa <- c("Ralstonia")

Anaero.colon.taxa <- c("Genus_12037",
                       "Senegalimassilia",
                       "Turicibacter")

TaxaVals.colon.shape <- rep(0,nrow(res))
names(TaxaVals.colon.shape) <- rep('NS',nrow(res))

TaxaVals.colon.shape[which(res$Genus %in% OAn.colon.taxa)] <- 16
names(TaxaVals.colon.shape)[which(res$Genus %in% OAn.colon.taxa)] <- 'Obligate Anaerobe'

TaxaVals.colon.shape[which(res$Genus %in% OA.colon.taxa)] <- 2
names(TaxaVals.colon.shape)[which(res$Genus %in% OA.colon.taxa)] <- 'Obligate Aerobe'

TaxaVals.colon.shape[which(res$Genus %in% FAn.colon.taxa)] <- 18
names(TaxaVals.colon.shape)[which(res$Genus %in% FAn.colon.taxa)] <- 'Facultative Anaerobe'

TaxaVals.colon.shape[which(res$Genus %in% Aero.colon.taxa)] <- 4
names(TaxaVals.colon.shape)[which(res$Genus %in% Aero.colon.taxa)] <- 'Aerobe'

TaxaVals.colon.shape[which(res$Genus %in% Anaero.colon.taxa)] <- 15
names(TaxaVals.colon.shape)[which(res$Genus %in% Anaero.colon.taxa)] <- 'Anaerobe'

unique(names(TaxaVals.colon.shape))
unique(TaxaVals.colon.shape)
TaxaVals.colon.shape[1:20]



#Rename
res$GS[res$GS == "Alloprevotella NA"] <- "Alloprevotella"
res$GS[res$GS == "Bacteroides Species_19535"] <- "Bacteroides sp."
res$GS[res$GS == "Barnesiella NA"] <- "Barnesiella"
res$GS[res$GS == "Blautia Species_1005"] <- "Blautia sp."
res$GS[res$GS == "Blautia Species_1038"] <- "Blautia sp."
res$GS[res$GS == "Blautia Species_336"] <- "Blautia sp."
res$GS[res$GS == "Blautia Species_641"] <- "Blautia sp."
res$GS[res$GS == "Blautia Species_721"] <- "Blautia sp."
res$GS[res$GS == "Clostridium_IV NA"] <- "Clostridium_IV"
res$GS[res$GS == "Clostridium_XlVa Species_2410"] <- "Clostridium_XlVa sp."
res$GS[res$GS == "Clostridium_XlVb Species_24049"] <- "Clostridium_XlVb sp."
res$GS[res$GS == "Clostridium_XVIII Species_10550"] <- "Clostridium_XVIII sp."
res$GS[res$GS == "Coprococcus Species_3251"] <- "Coprococcus sp."
res$GS[res$GS == "Coprococcus Species_3657"] <- "Coprococcus sp."
res$GS[res$GS == "Dialister Species_22017"] <- "Dialister sp."
res$GS[res$GS == "Dorea Species_26435"] <- "Dorea sp."
res$GS[res$GS == "Exiguobacterium NA"] <- "Exiguobacterium"
res$GS[res$GS == "Faecalibacterium Species_4580"] <- "Faecalibacterium sp."
res$GS[res$GS == "Faecalibacterium Species_4602"] <- "Faecalibacterium sp."
res$GS[res$GS == "Faecalibacterium Species_4606"] <- "Faecalibacteriums sp."
res$GS[res$GS == "Faecalibacterium Species_4928"] <- "Faecalibacterium sp."
res$GS[res$GS == "Faecalibacterium Species_4991"] <- "Faecalibacterium sp."
res$GS[res$GS == "Fusicatenibacter Species_396"] <- "Fusicatenibacter sp."
res$GS[res$GS == "Genus_12037 NA"] <- "Unknown Ruminococcaceae species"
res$GS[res$GS == "Genus_23841 NA"] <- "Unknown Ruminococcaceae species"
res$GS[res$GS == "Genus_25334 NA"] <- "Unknown Ruminococcaceae species"
res$GS[res$GS == "Genus_25335 NA"] <- "Unknown Lachnospiraceae species"
res$GS[res$GS == "Genus_25348 NA"] <- "Unknown Coriobacteriaceae species"
res$GS[res$GS == "Genus_27047 NA"] <- "Unknown Ruminococcaceae species"
res$GS[res$GS == "Genus_27241 NA"] <- "Unknown Lachnospiraceae species"
res$GS[res$GS == "Genus_6584 NA"] <- "Unknown Ruminococcaceae species"
res$GS[res$GS == "Klebsiella Species_16814"] <- "Klebsiella sp."
res$GS[res$GS == "Mogibacterium NA"] <- "Mogibacterium"
res$GS[res$GS == "Oscillibacter NA"] <- "Oscillibacter"
res$GS[res$GS == "Romboutsia Species_23587"] <- "Romboutsia sp."
res$GS[res$GS == "Romboutsia Species_23649"] <- "Romboutsia sp."
res$GS[res$GS == "Ruminococcus2 Species_25620"] <- "Ruminococcus2 sp."
res$GS[res$GS == "Ruminococcus2 Species_25636"] <- "Ruminococcus2 sp."
res$GS[res$GS == "Ruminococcus2 Species_26546"] <- "Ruminococcus2 sp."
res$GS[res$GS == "Senegalimassilia NA"] <- "Senegalimassilia"




sigTaxa.colon<- c("Alloprevotella",
                  "Bacteroides faecis",
                  "Bacteroides fragilis",
                  "Bacteroides ovatus",
                  "Bacteroides sp.",
                  "Barnesiella",
                  "Blautia sp.",
                  "Clostridium_IV",
                  "Clostridium_XlVa sp.",
                  "Clostridium_XlVb sp.",
                  "Clostridium_XVIII sp.",
                  "Coprococcus sp.",
                  "Dialister sp.",
                  "Dorea sp.",
                  "Exiguobacterium",
                  "Faecalibacterium sp.",
                  "Unknown Lachnospiraceae species",
                  "Unknown Ruminococcaceae species",
                  "Unknown Coriobacteriaceae species",
                  "Klebsiella sp.",
                  "Mogibacterium",
                  "Oscillibacter",
                  "Romboutsia sp.",
                  "Ruminococcus2 sp.",
                  "Senegalimassilia")

sigTaxa.colon.2 <- c("Alloprevotella",
                      "Bacteroides faecis",
                      "Bacteroides fragilis",
                      "Bacteroides ovatus",
                      "Bacteroides sp.",
                      "Barnesiella",
                      "Blautia sp.",
                      "Clostridium_IV",
                      "Clostridium_XlVa sp.",
                      "Clostridium_XlVb sp.",
                      "Clostridium_XVIII sp.",
                      "Faecalibacterium sp.",
                      "Unknown Lachnospiraceae species",
                      "Unknown Ruminococcaceae species",
                      "Unknown Coriobacteriaceae species",
                      "Klebsiella sp.",
                      "Ruminococcus2 sp.",
                      "Senegalimassilia")

sigTaxa.colon.3 <- c("Alloprevotella",
                     "Bacteroides faecis",
                     "Bacteroides fragilis",
                     "Bacteroides ovatus",
                     "Bacteroides sp.",
                     "Barnesiella",
                     "Blautia sp.",
                     "Clostridium_IV",
                     "Clostridium_XlVa sp.",
                     "Clostridium_XlVb sp.",
                     "Clostridium_XVIII sp.",
                     "Faecalibacterium sp.",
                     "Klebsiella sp.",
                     "Ruminococcus2 sp.")

res$log2FoldChange <- (-1)*res$log2FoldChange
colon.da <- EnhancedVolcano(res,
                             lab = res$GS,
                             x = 'log2FoldChange',
                             y = 'global.pvalue',
                             selectLab = sigTaxa.colon.3,
                             xlim = c(-6.5,6.5),
                             title = 'Alterations in the Colonic Mucosa',
                             ylab = bquote(~-Log[10]~global~italic(P)),
                             pCutoff = 0.030809408,
                             transcriptPointSize = res$log_grp_mra,
                             transcriptLabSize = 3.0,
                             transcriptLabCol = 'black',
                             transcriptLabFace = 'italic',
                             cutoffLineType = "twodash",
                             cutoffLineCol = "red4",
                             cutoffLineWidth = 1.0,
                             colAlpha = 1,
                             colCustom = TaxaVals,
                             shapeCustom = TaxaVals.colon.shape,
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

colon.da + coord_flip() 
