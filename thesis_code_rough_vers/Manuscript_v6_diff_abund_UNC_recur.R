###########################################
# Figure 2 - Ileal Disease - Genera level #
###########################################
devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)
res <- UNC_CD_recur_padj_mra_all_aero

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
OAn.UNC.taxa <- c("Bacteroides",
                  "Parabacteroides",
                  "Parasutterella",
                  "Anaerostipes",
                  "Blautia",
                  "Clostridium_IV",
                  "Clostridium_XlVa",
                  "Coprococcus",
                  "Dorea",
                  "Faecalibacterium",
                  "Genus_27050",
                  "Propionispira",
                  "Intestinibacter",
                  "Roseburia",
                  "Genus_55",
                  "Clostridium_sensu_stricto",
                  "Bifidobacterium",
                  "Prevotella",
                  "Anaerococcus",
                  "Butyricicoccus",
                  "Clostridium_XI",
                  "Flavonifractor",
                  "Fusicatenibacter",
                  "Genus_140",
                  "Genus_150",
                  "Genus_2035",
                  "Genus_25466",
                  "Genus_26240",
                  "Genus_27044",
                  "Genus_27125",
                  "Genus_27409",
                  "Genus_3202",
                  "Genus_342",
                  "Hungatella",
                  "Parvimonas",
                  "Peptostreptococcus",
                  "Pseudoflavonifractor",
                  "Romboutsia",
                  "Ruminococcus2",
                  "Subdoligranulum",
                  "Terrisporobacter",
                  "Bilophila",
                  "Solobacterium",
                  "Fusobacterium",
                  "Dialister",
                  "Megasphaera",
                  "Phascolarctobacterium")

OA.UNC.taxa <- c("Novosphingobium",
                 "Sphingobium",
                 "Sphingomonas",
                 "Chryseobacterium",
                 "Flavobacterium",
                 "Acinetobacter",
                 "Pedobacter")

FAn.UNC.taxa<- c("Actinomyces",
                 "Streptococcus",
                 "Citrobacter",
                 "Campylobacter",
                 "Corynebacterium",
                 "Gemella",
                 "Granulicatella",
                 "Genus_10876",
                 "Genus_10941",
                 "Enterobacter",
                 "Genus_16672",
                 "Haemophilus",
                 "Klebsiella")

Aero.UNC.taxa <- c("Rhodoluna",
                   "Rothia",
                   "Rhizobium",
                   "Bacillus",
                   "Acidovorax",
                   "Pseudomonas")

Anaero.UNC.taxa <- c("Sutterella",
                     "Veillonella",
                     "Propionibacterium",
                     "Turicibacter")

TaxaVals.UNC.shape <- rep(0,nrow(res))
names(TaxaVals.UNC.shape) <- rep('NS',nrow(res))

TaxaVals.UNC.shape[which(res$Genus %in% OAn.UNC.taxa)] <- 16
names(TaxaVals.UNC.shape)[which(res$Genus %in% OAn.UNC.taxa)] <- 'Obligate Anaerobe'

TaxaVals.UNC.shape[which(res$Genus %in% OA.UNC.taxa)] <- 2
names(TaxaVals.UNC.shape)[which(res$Genus %in% OA.UNC.taxa)] <- 'Obligate Aerobe'

TaxaVals.UNC.shape[which(res$Genus %in% FAn.UNC.taxa)] <- 18
names(TaxaVals.UNC.shape)[which(res$Genus %in% FAn.UNC.taxa)] <- 'Facultative Anaerobe'

TaxaVals.UNC.shape[which(res$Genus %in% Aero.UNC.taxa)] <- 4
names(TaxaVals.UNC.shape)[which(res$Genus %in% Aero.UNC.taxa)] <- 'Aerobe'

TaxaVals.UNC.shape[which(res$Genus %in% Anaero.UNC.taxa)] <- 15
names(TaxaVals.UNC.shape)[which(res$Genus %in% Anaero.UNC.taxa)] <- 'Anaerobe'

unique(names(TaxaVals.UNC.shape))
unique(TaxaVals.UNC.shape)
TaxaVals.UNC.shape[1:20]



#Rename

res$GS[res$GS == "Streptococcus Species_9162"] <- "Sterptococcus sp."
res$GS[res$GS == "Bacteroides Species_19765"] <- "Bacteroides sp."
res$GS[res$GS == "Bacteroides Species_19458"] <- "Bacteroides sp."
res$GS[res$GS == "Bacteroides Species_19761"] <- "Bacteroides sp."
res$GS[res$GS == "Bacteroides Species_19874"] <- "Bacteroides sp."	
res$GS[res$GS == "Bacteroides Species_20723"] <- "Bacteroides sp."
res$GS[res$GS == "Bacteroides Species_19625"] <- "Bacteroides sp."	
res$GS[res$GS == "Bacteroides Species_19555"] <- "Bacteroides sp."	
res$GS[res$GS == "Bacteroides Species_19435"] <- "Bacteroides sp."

res$GS[res$GS == "Sutterella NA"] <- "Sutterella"
res$GS[res$GS == "Anaerostipes Species_25777"] <- "Anaerostipes sp."
res$GS[res$GS == "Blautia Species_1036"] <- "Blautia sp."
res$GS[res$GS == "Clostridium_IV Species_6885"] <- "Clostridium_IV sp."
res$GS[res$GS == "Clostridium_XlVa Species_2623"] <- "Clostridium_XlVa sp."
res$GS[res$GS == "Clostridium_XlVa Species_2137"] <- "Clostridium_XlVa sp."
res$GS[res$GS == "Clostridium_XlVa Species_3534"] <- "Clostridium_XlVa sp."
res$GS[res$GS == "Clostridium_XlVa Species_2161"] <- "Clostridium_XlVa sp."
res$GS[res$GS == "Clostridium_XlVa Species_2227"] <- "Clostridium_XlVa sp."
res$GS[res$GS == "Coprococcus Species_3763"] <- "Coprococcus sp."
res$GS[res$GS == "Dorea Species_26325"] <- "Dorea sp."

"Genus_12037 NA"
"Propionibacterium NA"
"Pseudomonas NA"
"Methylobacterium NA"
"Delftia NA"
"Acidovorax NA"
"Comamonas NA"
"Acidocella NA"
"Rothia NA"
"Bradyrhizobium NA"
"Senegalimassilia NA"
"Veillonella NA"
"Genus_10447 NA"
"Genus_16621 NA"
"Genus_16796 NA"
"Shewanella NA"
"Enterobacter NA"
"Genus_10449 NA"
"Citrobacter NA"
"Genus_10911 NA"
"Staphylococcus NA"
"Enterococcus NA"
"Lactobacillus NA"
"Brevundimonas NA"
"Klebsiella NA"
"Haemophilus NA"
"Gemella NA"
"Streptococcus NA"
"Corynebacterium NA"
"Leuconostoc NA"
"Weissella NA"
"Granulicatella NA"
"Cellulosimicrobium NA"
"Genus_23869 NA"
"Genus_23849 NA"
"Genus_20845 NA"
"Pelomonas NA"
"Acinetobacter NA"
"Stenotrophomonas NA"
"Sphingobium NA"
"Massilia NA"
"Xanthomonas NA"
"Rhodanobacter NA"
"Genus_3203 NA"
"Peptococcus NA"
"Alloprevotella NA"
"Genus_5419 NA"
"Genus_23841 NA"
"Cellulosilyticum NA"
"Genus_5599 NA"
"Genus_3836 NA"
"Genus_3826 NA"
"Genus_27241 NA"
"Genus_6417 NA"
"Genus_25479 NA"
"Genus_24894 NA"
"Adlercreutzia NA"
"Parabacteroides NA"
"Alistipes NA"
"Clostridium_IV NA"
"Subdoligranulum NA"
"Collinsella NA"
"Genus_27125 NA"
"Ruminococcus NA"
"Genus_27276 NA"
"Fusicatenibacter NA"
"Pseudoflavonifractor NA"
"Genus_27153 NA"
"Akkermansia NA"
"Genus_25348 NA"
"Intestinimonas NA"
"Phascolarctobacterium NA"
"Genus_85 NA"
"Butyricimonas NA"
"Hathewaya NA"
"Genus_3886 NA"
"Genus_25334 NA"
"Genus_25335 NA"
"Anaerotruncus NA"
"Romboutsia NA"
"Paraprevotella NA"
"Genus_3107 NA"
"Genus_27244 NA"
"Genus_3202 NA"
"Clostridium_XlVa NA"
"Butyricicoccus NA"
"Genus_3948 NA"
"Genus_3835 NA"
"Porphyromonas NA"
"Prevotella NA"
"Eisenbergiella NA"
"Genus_3802 NA"
"Clostridium_XVIII NA"
"Anaerostipes NA"
"Gordonibacter NA"
"Megasphaera NA"
"Bifidobacterium NA"
"Genus_1426 NA"
"Dialister NA"
"Genus_1598 NA"
"Genus_1601 NA"
"Anaerococcus NA"
"Genus_27409 NA"
"Genus_118 NA"
"Peptoniphilus NA"
"Coprobacillus NA"
"Solobacterium NA"
"Terrisporobacter NA"
"Eggerthella NA"


sigTaxa.UNC<- c("Sterptococcus sp.",
                "Bacteroides sp.",
                "Bacteroides massiliensis",
                "Bacteroides xylanisolvens",
                "Anaerostipes sp.",
                "Clostridium_IV sp.",
                "Clostridium_XlVa sp.",
                "Coprococcus sp.",
                "Dorea sp.",
                "Blautia sp.",
                "Sutterella")

res$log2FoldChange <- (-1)*res$log2FoldChange
UNC.da <- EnhancedVolcano(res,
                            lab = res$GS,
                            x = 'log2FoldChange',
                            y = 'global.pvalue',
                            xlim = c(-6.5,6.5),
                            title = 'Alterations in the Recurrence - UNC',
                            ylab = bquote(~-Log[10]~global~italic(P)),
                            pCutoff = 0.018482336,
                            transcriptPointSize = res$log_grp_mra,
                            transcriptLabSize = 3.0,
                            transcriptLabCol = 'black',
                            transcriptLabFace = 'italic',
                            cutoffLineType = "twodash",
                            cutoffLineCol = "red4",
                            cutoffLineWidth = 1.0,
                            colAlpha = 1,
                            colCustom = TaxaVals,
                            shapeCustom = TaxaVals.UNC.shape,
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

UNC.da + coord_flip() 
