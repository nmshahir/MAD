###########################################
# Figure 2 - Ileal Disease - Genera level #
###########################################
devtools::install_github('kevinblighe/EnhancedVolcano')
res <- ileal_disease_padj_grp_mean_genera_fig_data

res$Genus <- as.character(res$Genus)
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
OAn.colon.taxa <- c("Akkermansia",
                    "Alistipes",
                    "Alloprevotella",
                    "Anaerococcus",
                    "Anaerosporobacter",
                    "Anaerotruncus",
                    "Anaerovorax",
                    "Bacteroides",
                    "Barnesiella",
                    "Bilophila",
                    "Blautia",
                    "Butyricicoccus",
                    "Butyricimonas",
                    "Butyrivibrio",
                    "Catabacter",
                    "Catonella",
                    "Cellulosilyticum",
                    "Clostridium_IV",
                    "Clostridium_sensu_stricto",
                    "Clostridium_XI",
                    "Clostridium_XlVa",
                    "Clostridium_XlVb",
                    "Clostridium_XVIII",
                    "Collinsella",
                    "Coprococcus",
                    "Dialister",
                    "Eggerthella",
                    "Eisenbergiella",
                    "Eubacterium",
                    "Faecalibacterium",
                    "Flavonifractor",
                    "Fusicatenibacter",
                    "Fusobacterium",
                    "Genus_140",
                    "Genus_150",
                    "Genus_24894",
                    "Genus_25348",
                    "Genus_25506",
                    "Genus_27125",
                    "Genus_27243",
                    "Genus_27309",
                    "Genus_27409",
                    "Genus_3926",
                    "Genus_5500",
                    "Genus_6417",
                    "Genus_7741",
                    "Genus_95",
                    "Gordonibacter",
                    "Hathewaya",
                    "Holdemanella",
                    "Holdemania",
                    "Hungatella",
                    "Intestinibacter",
                    "Intestinimonas",
                    "Lachnoanaerobaculum",
                    "Megasphaera",
                    "Mogibacterium",
                    "Oribacterium",
                    "Oscillibacter",
                    "Parabacteroides",
                    "Paraprevotella",
                    "Parasutterella",
                    "Parvimonas",
                    "Peptostreptococcus",
                    "Phascolarctobacterium",
                    "Porphyromonas",
                    "Prevotella",
                    "Pseudoflavonifractor",
                    "Romboutsia",
                    "Roseburia",
                    "Ruminococcus",
                    "Ruminococcus2",
                    "Slackia",
                    "Solobacterium",
                    "Stomatobaculum",
                    "Subdoligranulum",
                    "Terrisporobacter")
OA.colon.taxa <- c("Acinetobacter",
                   "Burkholderia",
                   "Elizabethkingia",
                   "Flavobacterium",
                   "Janthinobacterium",
                   "Massilia",
                   "Microbacterium",
                   "Stenotrophomonas",
                   "Xanthomonas")
FAn.colon.taxa<- c("Anoxybacillus",
                   "Brevundimonas",
                   "Propionibacterium",
                   "Lactobacillus",
                   "Citrobacter",
                   "Enterobacter",
                   "Enterococcus",
                   "Escherichia/Shigella",
                   "Gemella",
                   "Genus_16796",
                   "Haemophilus",
                   "Klebsiella",
                   "Lactococcus",
                   "Leptotrichia",
                   "Shewanella",
                   "Streptococcus",
                   "Campylobacter",
                   "Genus_10911",
                   "Genus_10933")
Aero.colon.taxa <- c("Acidovorax",
                     "Comamonas",
                     "Delftia",
                     "Legionella",
                     "Methylobacterium",
                     "Methylophilus",
                     "Ochrobactrum",
                     "Pedobacter",
                     "Polaromonas",
                     "Pseudarcicella",
                     "Pseudomonas",
                     "Rhodoluna",
                     "Undibacterium",
                     "Rothia")
Anaero.colon.taxa <- c("Anaeroglobus",
                       "Odoribacter",
                       "Turicibacter",
                       "Genus_12029",
                       "Atopobium")

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
res$Genus[res$Genus == "Genus_3203"] <- "Unknown Lachnospiraceae species"
res$Genus[res$Genus == "Genus_10911"] <- "Unknown Erysipelotrichaceae species"
res$Genus[res$Genus == "Genus_23841"] <- "Unknown Lachnospiraceae species"
res$Genus[res$Genus == "Genus_5613"] <- "Unknown Lachnospiraceae species"
res$Genus[res$Genus == "Genus_3826"] <- "Unknown Lachnospiraceae species"
res$Genus[res$Genus == "Genus_27125"] <-"Unknown Lachnospiraceae species"
res$Genus[res$Genus == "Genus_27153"] <-"Unknown Lachnospiraceae species"

sigTaxa.nonIBD <- c("Unknown Lachnospiraceae species",
                    "Unknown Erysipelotrichaceae species",
                    "Unknown Lachnospiraceae species",
                    "Unknown Lachnospiraceae species",
                    "Unknown Lachnospiraceae species",
                    "Anaerofilum",
                    "Roseburia",
                    "Blautia",
                    "Pseudoflavonifractor",
                    "Anaerostipes",
                    "Bacteroides",
                    "Unknown Lachnospiraceae species",
                    "Clostridium_XlVb",
                    "Flavonifractor",
                    "Unknown Lachnospiraceae species",
                    "Clostridium_XlVa")

res$log2FoldChange <- (-1)*res$log2FoldChange
nonIBD.da <- EnhancedVolcano(res,
                            lab = res$Genus,
                            x = 'log2FoldChange',
                            y = 'global.pvalue',
                            selectLab = sigTaxa.nonIBD,
                            xlim = c(-6.5,6.5),
                            title = 'Alterations in the NonIBD Mucosa - Genera Level',
                            ylab = bquote(~-Log[10]~global~italic(P)),
                            pCutoff = 0.009846661,
                            pointSize = res$log_grp_mean,
                            labSize = 3.0,
                            labCol = 'black',
                            labFace = 'italic',
                            cutoffLineType = "twodash",
                            cutoffLineCol = "red4",
                            cutoffLineWidth = 1.0,
                            colAlpha = 1,
                            colCustom = TaxaVals,
                            shapeCustom = TaxaVals.colon.shape,
                            boxedLabels = TRUE,
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

ileum.da + coord_flip() 
