###########################################
# Figure 2 - Colon Disease - Genera level #
###########################################
devtools::install_github('kevinblighe/EnhancedVolcano')
res <- colon_disease_padj_grp_mean

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

TaxaVals[which(res$Phylum == "Firmicutes")] <- "#3b8b00"
names(TaxaVals)[which(res$Phylum == "Firmicutes")] <- "Firmicutes"

TaxaVals[which(res$Phylum == "Proteobacteria")] <- "#72349e"
names(TaxaVals)[which(res$Phylum == "Proteobacteria")] <- "Proteobacteria"
TaxaVals[which(res$Phylum == "Verrucomicrobia")] <- "#01cacc"
names(TaxaVals)[which(res$Phylum == "Verrucomicrobia")] <- "Verrucomicrobia"

unique(names(TaxaVals))
TaxaVals[1:20]

#Setup ShapeVal 
#Obligate Anaerobes
OAn.colon.taxa <- c("Adlercreutzia",
                     "Collinsella",
                     "Eggerthella",
                     "Gordonibacter",
                     "Olsenella",
                     "Slackia",
                     "Bacteroides",
                     "Barnesiella",
                     "Coprobacter",
                     "Porphyromonas",
                     "Alloprevotella",
                     "Paraprevotella",
                     "Prevotella",
                     "Alistipes",
                     "Parasutterella",
                     "Catabacter",
                     "Clostridium_sensu_stricto",
                     "Hathewaya",
                     "Parvimonas",
                     "Anaerovorax",
                     "Mogibacterium",
                     "Anaerofustis",
                     "Eubacterium",
                     "Anaerostipes",
                     "Blautia",
                     "Cellulosilyticum",
                     "Clostridium_XlVa",
                     "Clostridium_XlVb",
                     "Coprococcus",
                     "Dorea",
                     "Eisenbergiella",
                     "Fusicatenibacter",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Unknown Lachnospiraceae genus",
                     "Hungatella",
                     "Murimonas",
                     "Ruminococcus2",
                     "Unknown Clostridiales genus",
                     "Unknown Clostridiales genus",
                     "Peptococcus",
                     "Peptoniphilus",
                     "Intestinibacter",
                     "Peptostreptococcus",
                     "Romboutsia",
                     "Terrisporobacter",
                     "Anaerofilum",
                     "Anaerotruncus",
                     "Butyricicoccus",
                     "Clostridium_IV",
                     "Ethanoligenens",
                     "Faecalibacterium",
                     "Flavonifractor",
                     "Unknown Ruminococcaceae genus",
                     "Unknown Ruminococcaceae genus",
                     "Unknown Ruminococcaceae genus",
                     "Unknown Ruminococcaceae genus",
                     "Unknown Ruminococcaceae genus",
                     "Unknown Ruminococcaceae genus",
                     "Unknown Ruminococcaceae genus",
                     "Unknown Ruminococcaceae genus",
                     "Unknown Ruminococcaceae genus",
                     "Unknown Ruminococcaceae genus",
                     "Intestinimonas",
                     "Oscillibacter",
                     "Ruminococcus",
                     "Subdoligranulum",
                     "Desulfovibrio",
                     "Holdemania",
                     "Solobacterium",
                     "Phascolarctobacterium",
                     "Dialister",
                     "Propionispira",
                     "Akkermansia")
OA.colon.taxa <- c("Burkholderia",
                    "Janthinobacterium",
                    "Oxalobacter",
                    "Acinetobacter",
                    "Stenotrophomonas",
                    "Xanthomonas")
FAn.colon.taxa<- c("Gemella",
                    "Granulicatella",
                    "Enterococcus",
                    "Leuconostoc",
                    "Lactococcus",
                    "Streptococcus",
                    "Citrobacter",
                    "Escherichia/Shigella",
                    "Unknown Enterobacteriaceae genus",
                    "Klebsiella",
                    "Haemophilus",
                    "Actinomyces",
                    "Unknown Erysipelotrichaceae genus",
                    "Unknown Erysipelotrichaceae genus",
                    "Unknown Erysipelotrichaceae genus",
                    "Unknown Erysipelotrichaceae genus")
Aero.colon.taxa <- c("Rothia",
                     "Rhodococcus",
                     "Pseudomonas")
Anaero.colon.taxa <- c("Bifidobacterium",
                       "Unknown Coriobacteriaceae genus",
                       "Unknown Coriobacteriaceae genus",
                       "Unknown Coriobacteriaceae genus",
                       "Odoribacter",
                       "Clostridium_XVIII",
                       "Coprobacillus",
                       "Turicibacter",
                       "Acidaminococcus",
                       "Veillonella")
TaxaVals.colon.shape <- rep(0,nrow(res))
names(TaxaVals.colon.shape) <- rep('NS',nrow(res))

TaxaVals.colon.shape[which(res$GenusLabel %in% OAn.colon.taxa)] <- 16
names(TaxaVals.colon.shape)[which(res$GenusLabel %in% OAn.colon.taxa)] <- 'Obligate Anaerobe'

TaxaVals.colon.shape[which(res$GenusLabel %in% OA.colon.taxa)] <- 2
names(TaxaVals.colon.shape)[which(res$GenusLabel %in% OA.colon.taxa)] <- 'Obligate Aerobe'

TaxaVals.colon.shape[which(res$GenusLabel %in% FAn.colon.taxa)] <- 18
names(TaxaVals.colon.shape)[which(res$GenusLabel %in% FAn.colon.taxa)] <- 'Facultative Anaerobe'

TaxaVals.colon.shape[which(res$GenusLabel %in% Aero.colon.taxa)] <- 4
names(TaxaVals.colon.shape)[which(res$GenusLabel %in% Aero.colon.taxa)] <- 'Aerobe'

TaxaVals.colon.shape[which(res$GenusLabel %in% Anaero.colon.taxa)] <- 15
names(TaxaVals.colon.shape)[which(res$GenusLabel %in% Anaero.colon.taxa)] <- 'Anaerobe'

unique(names(TaxaVals.colon.shape))
unique(TaxaVals.colon.shape)
TaxaVals.colon.shape[1:20]

res$Genus <- as.character(res$Genus)
res$Genus[res$Genus == "Genus_23841"] <- "Unknown Ruminococcaceae genus"
res$Genus[res$Genus == "Genus_25348"] <- "Unknown Lachnospiraceae genus"
res$Genus[res$Genus == "Genus_6584"] <- "Unknown Ruminococcaceae genus"
res$Genus[res$Genus == "Genus_12037"] <- "Unknown Coriobacteriaceae genus"
#Identify Signifcant Taxa
sigTaxa.colon <- c("Anaerofilum",
                   "Peptococcus",
                   "Unknown Ruminococcaceae genus",
                   "Unknown Lachnospiraceae genus",
                   "Alloprevotella",
                   "Unknown Ruminococcaceae genus",
                   "Unknown Coriobacteriaceae genus",
                   "Mogibacterium")
res$log2FoldChange <- (-1)*res$log2FoldChange
colon.da <- EnhancedVolcano(res,
                            lab = res$Genus,
                            x = 'log2FoldChange',
                            y = 'global.pvalue',
                            selectLab = sigTaxa.colon,
                            xlim = c(-6.5,6.5),
                            title = 'Alterations in the Colonic Mucosa - Genera Level',
                            ylab = bquote(~-Log[10]~global~italic(P)),
                            pCutoff = 0.005984966,
                            transcriptPointSize = res$log_grp_mean,
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
