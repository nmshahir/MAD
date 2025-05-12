###########################################
# Figure 2 - Ileal Disease - Genera level #
###########################################
devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)
res <- ileum_disease_merge_fig_data_gen_sp_aero

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
OAn.ileum.taxa <- c("Bacteroides",
                    "Bifidobacterium",
                    "Blautia",
                    "Clostridium_IV",
                    "Clostridium_XlVa",
                    "Clostridium_XVIII",
                    "Coprococcus",
                    "Flavonifractor",
                    "Genus_150",
                    "Genus_27409",
                    "Genus_95",
                    "Intestinibacter",
                    "Parabacteroides",
                    "Roseburia",
                    "Ruminococcus2",
                    "Alistipes",
                    "Alloprevotella",
                    "Anaerovorax",
                    "Butyricimonas",
                    "Catabacter",
                    "Clostridium_sensu_stricto",
                    "Clostridium_XlVb",
                    "Faecalibacterium",
                    "Genus_24885",
                    "Genus_24894",
                    "Genus_25131",
                    "Genus_25506",
                    "Genus_27276",
                    "Genus_27309",
                    "Genus_3107",
                    "Genus_6417",
                    "Megasphaera",
                    "Oscillibacter",
                    "Prevotella",
                    "Romboutsia",
                    "Ruminococcus",
                    "Stomatobaculum",
                    "Subdoligranulum",
                    "Barnesiella",
                    "Paraprevotella",
                    "Porphyromonas",
                    "Parasutterella",
                    "Anaerostipes",
                    "Anaerotruncus",
                    "Butyrivibrio",
                    "Clostridium_XI",
                    "Dorea",
                    "Eisenbergiella",
                    "Fusicatenibacter",
                    "Genus_140",
                    "Genus_2035",
                    "Genus_25335",
                    "Genus_25348",
                    "Genus_25466",
                    "Genus_25471",
                    "Genus_25598",
                    "Genus_27044",
                    "Genus_27177",
                    "Genus_27243",
                    "Genus_3203",
                    "Genus_3867",
                    "Genus_3926",
                    "Genus_49",
                    "Genus_55",
                    "Hungatella",
                    "Oribacterium",
                    "Parvimonas",
                    "Peptostreptococcus",
                    "Terrisporobacter",
                    "Atopobium",
                    "Collinsella",
                    "Eggerthella",
                    "Bilophila",
                    "Desulfovibrio",
                    "Holdemanella",
                    "Fusobacterium",
                    "Dialister",
                    "Phascolarctobacterium",
                    "Propionispira",
                    "Selenomonas",
                    "Akkermansia")

OA.ileum.taxa <- c("Stenotrophomonas",
                   "Burkholderia",
                   "Chryseobacterium",
                   "Flavobacterium",
                   "Janthinobacterium",
                   "Massilia",
                   "Elizabethkingia",
                   "Acinetobacter",
                   "Novosphingobium",
                   "Sphingobium",
                   "Stenotrophomonas")

FAn.ileum.taxa<- c("Enterococcus",
                   "Genus_10911",
                   "Genus_16672",
                   "Haemophilus",
                   "Klebsiella",
                   "Raoultella",
                   "Actinomyces",
                   "Anoxybacillus",
                   "Genus_16796",
                   "Lactobacillus",
                   "Staphylococcus",
                   "Streptococcus",
                   "Gemella",
                   "Citrobacter",
                   "Enterobacter",
                   "Genus_16629",
                   "Genus_16873",
                   "Genus_10941",
                   "Granulicatella",
                   "Allisonella")

Aero.ileum.taxa <- c("Delftia",
                     "Duganella",
                     "Legionella",
                     "Methylobacterium",
                     "Methylotenera",
                     "Ochrobactrum",
                     "Pseudarcicella",
                     "Rhodoluna",
                     "Mycobacterium",
                     "Rothia",
                     "Paenibacillus",
                     "Acidovorax",
                     "Pseudomonas")

Anaero.ileum.taxa <- c("Odoribacter",
                       "Sutterella",
                       "Veillonella",
                       "Turicibacter")

TaxaVals.ileum.shape <- rep(0,nrow(res))
names(TaxaVals.ileum.shape) <- rep('NS',nrow(res))

TaxaVals.ileum.shape[which(res$Genus %in% OAn.ileum.taxa)] <- 16
names(TaxaVals.ileum.shape)[which(res$Genus %in% OAn.ileum.taxa)] <- 'Obligate Anaerobe'

TaxaVals.ileum.shape[which(res$Genus %in% OA.ileum.taxa)] <- 2
names(TaxaVals.ileum.shape)[which(res$Genus %in% OA.ileum.taxa)] <- 'Obligate Aerobe'

TaxaVals.ileum.shape[which(res$Genus %in% FAn.ileum.taxa)] <- 18
names(TaxaVals.ileum.shape)[which(res$Genus %in% FAn.ileum.taxa)] <- 'Facultative Anaerobe'

TaxaVals.ileum.shape[which(res$Genus %in% Aero.ileum.taxa)] <- 4
names(TaxaVals.ileum.shape)[which(res$Genus %in% Aero.ileum.taxa)] <- 'Aerobe'

TaxaVals.ileum.shape[which(res$Genus %in% Anaero.ileum.taxa)] <- 15
names(TaxaVals.ileum.shape)[which(res$Genus %in% Anaero.ileum.taxa)] <- 'Anaerobe'

unique(names(TaxaVals.ileum.shape))
unique(TaxaVals.ileum.shape)
TaxaVals.ileum.shape[1:20]



#Rename

res$GS[res$GS == "Bacteroides NA"] <- "Bacteroides"
res$GS[res$GS == "Blautia Species_542"] <- "Blautia sp."
res$GS[res$GS == "Blautia Species_502"] <- "Blautia sp."
res$GS[res$GS == "Blautia Species_3455"] <- "Blautia sp."
res$GS[res$GS == "Blautia Species_3459"] <- "Blautia sp."
res$GS[res$GS == "Clostridium_IV Species_6888"] <- "Clostridium_IV sp."
res$GS[res$GS == "Clostridium_XlVa NA"] <- "Clostridium_XlVa"
res$GS[res$GS == "Clostridium_XVIII NA"] <- "Clostridium_XVIII"
res$GS[res$GS == "Coprococcus Species_3674"] <- "Coprococcus sp."
res$GS[res$GS == "Coprococcus Species_3634"] <- "Coprococcus sp."
res$GS[res$GS == "Coprococcus Species_3710"] <- "Coprococcus sp."
res$GS[res$GS == "Coprococcus Species_3716"] <- "Coprococcus sp."
res$GS[res$GS == "Coprococcus Species_3726"] <- "Coprococcus sp."
res$GS[res$GS == "Delftia NA"] <- "Delftia"
res$GS[res$GS == "Enterococcus Species_9849"] <- "Enterococcus sp."
res$GS[res$GS == "Genus_10911 NA"] <- "Unknown Erysipelotrichaceae species"
res$GS[res$GS == "Genus_150 NA"] <- "Unknown Lachnospiraceae species"
res$GS[res$GS == "Genus_15493 NA"] <- "Unknown Neisseriaceae species"
res$GS[res$GS == "Genus_16672 NA"] <- "Unknown Enterobacteriaceae species"
res$GS[res$GS == "Genus_27409 NA"] <- "Unknown Lachnospiraceae species"
res$GS[res$GS == "Genus_95 NA"] <- "Unknown Lachnospiraceae species"
res$GS[res$GS == "Haemophilus NA"] <- "Haemophilus"
res$GS[res$GS == "Intestinibacter Species_23381"] <- "Intestinibacter sp."
res$GS[res$GS == "Klebsiella NA"] <- "Klebsiella"
res$GS[res$GS == "Parabacteroides Species_19022"] <- "Parabacteroides sp."
res$GS[res$GS == "Parabacteroides Species_19004"] <- "Parabacteroides sp."
res$GS[res$GS == "Raoultella NA"] <- "Raoultella"
res$GS[res$GS == "Roseburia NA"] <- "Roseburia"
res$GS[res$GS == "Ruminococcus2 Species_26538"] <- "Ruminococcus2 sp."
res$GS[res$GS == "Actinomyces Species_13071"] <- "Actinomyces sp."
res$GS[res$GS == "Alistipes Species_18607"] <- "Alistipes sp."
res$GS[res$GS == "Alloprevotella NA"] <- "Alloprevotella"
res$GS[res$GS == "Anaerovorax NA"] <- "Anaerovorax"
res$GS[res$GS == "Anoxybacillus NA"] <- "Anoxybacillus"
res$GS[res$GS == "Bacteroides Species_20706"] <- "Bacteroides sp."
res$GS[res$GS == "Bacteroides Species_19164"] <- "Bacteroides sp."
res$GS[res$GS == "Blautia Species_1161"] <- "Blautia sp."
res$GS[res$GS == "Blautia Species_890"] <- "Blautia sp."
res$GS[res$GS == "Burkholderia NA"] <- "Burkholderia"
res$GS[res$GS == "Butyricimonas NA"] <- "Butyricimonas"
res$GS[res$GS == "Catabacter NA"] <- "Catabacter"
res$GS[res$GS == "Chryseobacterium Species_21434"] <- "Chryseobacterium sp."
res$GS[res$GS == "Clostridium_sensu_stricto Species_22884"] <- "Clostridium_sensu_stricto sp."
res$GS[res$GS == "Clostridium_sensu_stricto Species_23287"] <- "Clostridium_sensu_stricto sp."
res$GS[res$GS == "Clostridium_XlVb Species_23989"] <- "Clostridium_XlVb sp."
res$GS[res$GS == "Duganella NA"] <- "Duganella"
res$GS[res$GS == "Faecalibacterium Species_4610"] <- "Faecalibacterium sp."
res$GS[res$GS == "Faecalibacterium Species_4947"] <- "Faecalibacterium sp."
res$GS[res$GS == "Flavobacterium NA"] <- "Flavobacterium"
res$GS[res$GS == "Genus_16796 NA"] <- "Unknown Enterobacteriaceae species"
res$GS[res$GS == "Genus_20845 NA"] <- "Unknown Rikenellaceae species"
res$GS[res$GS == "Genus_24885 NA"] <- "Unknown Clostridiales species"
res$GS[res$GS == "Genus_24894 NA"] <- "Unknown Clostridiales species"
res$GS[res$GS == "Genus_25131 NA"] <- "Unknown Clostridiales species"
res$GS[res$GS == "Genus_25506 NA"] <- "Unknown Lachnospiraceae species"
res$GS[res$GS == "Genus_27276 NA"] <- "Unknown Lachnospiraceae species"
res$GS[res$GS == "Genus_27309 NA"] <- "Unknown Lachnospiraceae species"
res$GS[res$GS == "Genus_3107 NA"] <- "Unknown Lachnospiraceae species"
res$GS[res$GS == "Genus_6417 NA"] <- "Unknown Ruminococcaceae species"
res$GS[res$GS == "Janthinobacterium Species_15402"] <- "Janthinobacterium sp."
res$GS[res$GS == "Legionella NA"] <- "Legionella"
res$GS[res$GS == "Massilia Species_15388"] <- "Massilia sp."
res$GS[res$GS == "Megamonas NA"] <- "Megamonas"
res$GS[res$GS == "Megasphaera NA"] <- "Megasphaera"
res$GS[res$GS == "Methylobacterium Species_13962"] <- "Methylobacterium sp."
res$GS[res$GS == "Methylobacterium Species_13960"] <- "Methylobacterium sp."
res$GS[res$GS == "Odoribacter NA"] <- "Odoribacter"
res$GS[res$GS == "Oscillibacter NA"] <- "Oscillibacter"
res$GS[res$GS == "Prevotella Species_20272"] <- "Prevotella sp."
res$GS[res$GS == "Prevotella Species_20317"] <- "Prevotella sp."
res$GS[res$GS == "Prevotella Species_20381"] <- "Prevotella sp."
res$GS[res$GS == "Pseudarcicella NA"] <- "Pseudarcicella"
res$GS[res$GS == "Rhodoluna NA"] <- "Rhodoluna"
res$GS[res$GS == "Romboutsia Species_23651"] <- "Romboutsia sp."
res$GS[res$GS == "Romboutsia Species_23564"] <- "Romboutsia sp."
res$GS[res$GS == "Ruminococcus Species_6614"] <- "Ruminococcus sp."
res$GS[res$GS == "Staphylococcus Species_11527"] <- "Staphylococcus sp."
res$GS[res$GS == "Staphylococcus Species_11519"] <- "Staphylococcus sp."
res$GS[res$GS == "Stomatobaculum NA"] <- "Stomatobaculum"
res$GS[res$GS == "Streptococcus Species_8594"] <- "Streptococcus sp."
res$GS[res$GS == "Streptococcus Species_9170"] <- "Streptococcus sp."
res$GS[res$GS == "Subdoligranulum NA"] <- "Subdoligranulum"
res$GS[res$GS == "Sutterella Species_14902"] <- "Sutterella sp."
res$GS[res$GS == "Veillonella Species_21697"] <- "Veillonella sp."




sigTaxa.ileum<- c("Actinomyces sp.",
                  "Alistipes finegoldii",
                  "Alistipes putredinis",
                  "Alistipes sp.",
                  "Alloprevotella",
                  "Anaerovorax",
                  "Anoxybacillus",
                  "Bacteroides fragilis",
                  "Bacteroides",
                  "Bacteroides sp.",
                  "Bifidobacterium pseudocatenulatum",
                  "Blautia sp.",
                  "Burkholderia",
                  "Butyricimonas",
                  "Catabacter",
                  "Chryseobacterium sp.",
                  "Clostridium_IV sp.",
                  "Clostridium_sensu_stricto sp.",
                  "Clostridium_XlVa",
                  "Clostridium_XlVb sp.",
                  "Clostridium_XVIII",
                  "Coprococcus sp.",
                  "Delftia",
                  "Duganella",
                  "Enterococcus faecium",
                  "Enterococcus sp.",
                  "Faecalibacterium sp.",
                  "Flavobacterium",
                  "Flavonifractor plautii",
                  "Unknown Lachnospiraceae species",
                  "Unknown Ruminococcaceae species",
                  "Unknown Clostridiales species",
                  "Unknown Enterobacteriaceae species",
                  "Unknown Rikenellaceae species",
                  "Haemophilus",
                  "Intestinibacter sp.",
                  "Janthinobacterium sp.",
                  "Klebsiella",
                  "Lactobacillus fermentum",
                  "Legionella",
                  "Massilia sp.",
                  "Megamonas",
                  "Megasphaera",
                  "Methylobacterium sp.",
                  "Methylotenera versatilis",
                  "Ochrobactrum tritici",
                  "Odoribacter",
                  "Oscillibacter",
                  "Parabacteroides sp.",
                  "Prevotella sp.",
                  "Pseudarcicella",
                  "Raoultella",
                  "Rhodoluna",
                  "Romboutsia sp.",
                  "Roseburia",
                  "Ruminococcus sp.",
                  "Ruminococcus2 sp.",
                  "Staphylococcus sp.",
                  "Stenotrophomonas maltophilia",
                  "Stomatobaculum",
                  "Streptococcus anginosus",
                  "Streptococcus sp.",
                  "Subdoligranulum",
                  "Sutterella sp.",
                  "Veillonella sp.")

sigTaxa.ileum.2 <- c("Bacteroides fragilis",
                  "Bacteroides",
                  "Clostridium_XlVa",
                  "Clostridium_XVIII",
                  "Coprococcus sp.",
                  "Enterococcus faecium",
                  "Enterococcus sp.",
                  "Faecalibacterium sp.",
                  "Unknown Lachnospiraceae species",
                  "Unknown Ruminococcaceae species",
                  "Unknown Clostridiales species",
                  "Unknown Enterobacteriaceae species",
                  "Haemophilus",
                  "Klebsiella")

res$log2FoldChange <- (-1)*res$log2FoldChange
ileum.da <- EnhancedVolcano(res,
                            lab = res$GS,
                            x = 'log2FoldChange',
                            y = 'global.pvalue',
                            selectLab = sigTaxa.ileum.2,
                            xlim = c(-6.5,6.5),
                            title = 'Alterations in the Ileal Mucosa',
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
                            shapeCustom = TaxaVals.ileum.shape,
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

ileum.da + coord_flip() 
