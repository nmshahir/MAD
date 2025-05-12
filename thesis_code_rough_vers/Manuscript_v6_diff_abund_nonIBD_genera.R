###########################################
# Figure 2 - Ileal Disease - Genera level #
###########################################
devtools::install_github('kevinblighe/EnhancedVolcano')
res <- nonIBD_tissue_padj_grp_genera_mean_fig_data

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
OAn.colon.taxa <- c("Oribacterium",
                    "Genus_3203",
                    "Genus_23841",
                    "Genus_5613",
                    "Genus_3826",
                    "Anaerofilum",
                    "Roseburia",
                    "Blautia",
                    "Pseudoflavonifractor",
                    "Anaerostipes",
                    "Bacteroides",
                    "Genus_27125",
                    "Clostridium_XlVb",
                    "Flavonifractor",
                    "Genus_27153",
                    "Fusobacterium",
                    "Clostridium_XlVa",
                    "Anaerosporobacter",
                    "Anaerotruncus",
                    "Clostridium_IV",
                    "Genus_150",
                    "Genus_7287",
                    "Genus_6584",
                    "Eggerthella",
                    "Genus_25458",
                    "Alloprevotella",
                    "Adlercreutzia",
                    "Genus_95",
                    "Ruminococcus",
                    "Anaerofustis",
                    "Butyricicoccus",
                    "Parabacteroides",
                    "Genus_5599",
                    "Clostridium_XVIII",
                    "Holdemania",
                    "Genus_6661",
                    "Genus_25506",
                    "Lachnoanaerobaculum",
                    "Stomatobaculum",
                    "Parvimonas",
                    "Genus_1595",
                    "Ruminococcus2",
                    "Eisenbergiella",
                    "Genus_3344",
                    "Genus_5576",
                    "Genus_5662",
                    "Genus_26491",
                    "Intestinimonas",
                    "Solobacterium",
                    "Intestinibacter",
                    "Genus_3886",
                    "Parasutterella",
                    "Genus_3202",
                    "Phascolarctobacterium",
                    "Alistipes",
                    "Coprococcus",
                    "Bifidobacterium",
                    "Genus_25334",
                    "Murimonas",
                    "Collinsella",
                    "Genus_3387",
                    "Genus_3845",
                    "Prevotella",
                    "Fusicatenibacter",
                    "Genus_6159",
                    "Faecalibacterium",
                    "Romboutsia",
                    "Genus_27309",
                    "Oscillibacter",
                    "Genus_26218",
                    "Genus_7432",
                    "Genus_24894",
                    "Hungatella",
                    "Genus_3871",
                    "Holdemanella",
                    "Dorea",
                    "Genus_25466",
                    "Genus_3116",
                    "Peptococcus",
                    "Subdoligranulum",
                    "Genus_27276",
                    "Cellulosilyticum",
                    "Genus_27248",
                    "Butyrivibrio",
                    "Genus_6417",
                    "Barnesiella",
                    "Anaerovorax",
                    "Coprobacter",
                    "Butyricimonas",
                    "Hathewaya",
                    "Mogibacterium",
                    "Genus_5511",
                    "Gordonibacter",
                    "Genus_3802",
                    "Genus_27243",
                    "Genus_25348",
                    "Clostridium_sensu_stricto",
                    "Ethanoligenens",
                    "Genus_1601",
                    "Paraprevotella",
                    "Olsenella",
                    "Eubacterium",
                    "Dialister",
                    "Desulfovibrio",
                    "Catabacter",
                    "Slackia")
OA.colon.taxa <- c("Flavobacterium",
                   "Elizabethkingia",
                   "Chryseobacterium",
                   "Acinetobacter",
                   "Burkholderia",
                   "Janthinobacterium",
                   "Sphingomonas",
                   "Massilia",
                   "Stenotrophomonas",
                   "Rhodanobacter")

FAn.colon.taxa<- c("Genus_10911",
                   "Granulicatella",
                   "Staphylococcus",
                   "Streptococcus",
                   "Propionibacterium",
                   "Gemella",
                   "Genus_10447",
                   "Genus_10933",
                   "Genus_16796",
                   "Genus_10545",
                   "Actinomyces",
                   "Enterobacter",
                   "Genus_10449",
                   "Leptotrichia",
                   "Brevundimonas",
                   "Escherichia/Shigella",
                   "Lactobacillus",
                   "Corynebacterium",
                   "Citrobacter",
                   "Haemophilus")

Aero.colon.taxa <- c("Pseudomonas",
                     "Legionella",
                     "Pedobacter",
                     "Methylobacterium",
                     "Pseudarcicella",
                     "Acidovorax",
                     "Rhodoluna",
                     "Ochrobactrum",
                     "Methylophilus",
                     "Undibacterium",
                     "Rothia")

Anaero.colon.taxa <- c("Veillonella",
                       "Atopobium",
                       "Genus_12029",
                       "Sutterella",
                       "Odoribacter")

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

nonIBD.da + coord_flip() 
