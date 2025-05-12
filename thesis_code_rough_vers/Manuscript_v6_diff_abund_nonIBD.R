###########################################
# Figure 2 - Ileal Disease - Genera level #
###########################################
devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)
res <- nonIBD_tissue_fig_data_gen_sp_fin

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
OAn.nonIBD.taxa <- c("Collinsella",
                 "Bacteroides",
                 "Barnesiella",
                 "Alistipes",
                 "Eubacterium",
                 "Anaerostipes",
                 "Blautia",
                 "Clostridium_XlVa",
                 "Clostridium_XlVb",
                 "Coprococcus",
                 "Dorea",
                 "Eisenbergiella",
                 "Fusicatenibacter",
                 "Genus_150",
                 "Genus_25479",
                 "Genus_27125",
                 "Genus_27153",
                 "Genus_27241",
                 "Genus_3203",
                 "Genus_3826",
                 "Genus_3836",
                 "Roseburia",
                 "Ruminococcus2",
                 "Anaerotruncus",
                 "Clostridium_IV",
                 "Faecalibacterium",
                 "Flavonifractor",
                 "Genus_23841",
                 "Oscillibacter",
                 "Pseudoflavonifractor",
                 "Ruminococcus",
                 "Parasutterella",
                 "Prevotella",
                 "Clostridium_sensu_stricto",
                 "Parvimonas",
                 "Oribacterium",
                 "Romboutsia",
                 "Megasphaera",
                 "Fusobacterium",
                 "Bifidobacterium",
                 "Adlercreutzia",
                 "Eggerthella",
                 "Parabacteroides",
                 "Mogibacterium",
                 "Cellulosilyticum",
                 "Genus_118",
                 "Genus_1595",
                 "Genus_1598",
                 "Genus_1601",
                 "Genus_27223",
                 "Genus_27248",
                 "Genus_27318",
                 "Genus_3202",
                 "Genus_3886",
                 "Hungatella",
                 "Peptococcus",
                 "Subdoligranulum",
                 "Clostridium_XVIII",
                 "Holdemania",
                 "Phascolarctobacterium",
                 "Bilophila",
                 "Butyricimonas",
                 "Alloprevotella",
                 "Hathewaya",
                 "Butyrivibrio",
                 "Lachnoanaerobaculum",
                 "Stomatobaculum",
                 "Genus_24885",
                 "Solobacterium",
                 "Dialister")

OA.nonIBD.taxa <- c("Chryseobacterium",
                "Elizabethkingia",
                "Flavobacterium",
                "Burkholderia",
                "Janthinobacterium",
                "Sphingomonas",
                "Massilia",
                "Rhodanobacter",
                "Stenotrophomonas")

FAn.nonIBD.taxa<- c("Genus_10911",
                "Actinomyces",
                "Propionibacterium",
                "Anoxybacillus",
                "Gemella",
                "Staphylococcus",
                "Granulicatella",
                "Streptococcus",
                "Enterobacter",
                "Genus_10447",
                "Genus_10933",
                "Genus_10941",
                "Citrobacter",
                "Rothia",
                "Lactobacillus",
                "Genus_16796",
                "Haemophilus")

Aero.nonIBD.taxa <- c("Mycobacterium",
                  "Legionella",
                  "Pseudomonas",
                  "Rhodoluna",
                  "Ochrobactrum",
                  "Methylobacterium")

Anaero.nonIBD.taxa <- c("Atopobium",
                    "Odoribacter",
                    "Veillonella",
                    "Senegalimassilia",
                    "Sutterella")

TaxaVals.nonIBD.shape <- rep(0,nrow(res))
names(TaxaVals.nonIBD.shape) <- rep('NS',nrow(res))

TaxaVals.nonIBD.shape[which(res$Genus %in% OAn.nonIBD.taxa)] <- 16
names(TaxaVals.nonIBD.shape)[which(res$Genus %in% OAn.nonIBD.taxa)] <- 'Obligate Anaerobe'

TaxaVals.nonIBD.shape[which(res$Genus %in% OA.nonIBD.taxa)] <- 2
names(TaxaVals.nonIBD.shape)[which(res$Genus %in% OA.nonIBD.taxa)] <- 'Obligate Aerobe'

TaxaVals.nonIBD.shape[which(res$Genus %in% FAn.nonIBD.taxa)] <- 18
names(TaxaVals.nonIBD.shape)[which(res$Genus %in% FAn.nonIBD.taxa)] <- 'Facultative Anaerobe'

TaxaVals.nonIBD.shape[which(res$Genus %in% Aero.nonIBD.taxa)] <- 4
names(TaxaVals.nonIBD.shape)[which(res$Genus %in% Aero.nonIBD.taxa)] <- 'Aerobe'

TaxaVals.nonIBD.shape[which(res$Genus %in% Anaero.nonIBD.taxa)] <- 15
names(TaxaVals.nonIBD.shape)[which(res$Genus %in% Anaero.nonIBD.taxa)] <- 'Anaerobe'

unique(names(TaxaVals.nonIBD.shape))
unique(TaxaVals.nonIBD.shape)
TaxaVals.nonIBD.shape[1:20]



#Rename
res$GS[res$GS == "Genus_10911 NA"] <- "Unknown Erysipelotrichaceae sp."
res$GS[res$GS == "Genus_150 NA"] <- "Unknown Lachnospiraceae sp."
res$GS[res$GS == "Genus_23841 NA"] <- "Unknown Ruminococcaceae sp."
res$GS[res$GS == "Genus_25479 NA"] <- "Unknown Lachnospiraceae sp."
res$GS[res$GS == "Genus_27125 NA"] <- "Unknown Lachnospiraceae sp."
res$GS[res$GS == "Genus_27153 NA"] <- "Unknown Lachnospiraceae sp."
res$GS[res$GS == "Genus_27241 NA"] <- "Unknown Lachnospiraceae sp."
res$GS[res$GS == "Genus_3203 NA"] <- "Unknown Lachnospiraceae sp."
res$GS[res$GS == "Genus_3826 NA"] <- "Unknown Lachnospiraceae sp."
res$GS[res$GS == "Genus_3836 NA"] <- "Unknown Lachnospiraceae sp."

res$GS[res$GS =="Bacteroides NA"] <- "Bacteroides"
res$GS[res$GS =="Anaerostipes NA"] <- "Anaerostipes"
res$GS[res$GS =="Blautia NA"] <- "Blautia"
res$GS[res$GS =="Clostridium_XlVa NA"] <- "Clostridium_XlVa"
res$GS[res$GS =="Clostridium_XlVb NA"] <- "Clostridium_XlVb"
res$GS[res$GS =="Anaerotruncus NA"] <- "Anaerotruncus"
res$GS[res$GS =="Clostridium_IV NA"] <- "Clostridium_IV"
res$GS[res$GS =="Flavonifractor NA"] <- "Flavonifractor"
res$GS[res$GS =="Pseudoflavonifractor NA"] <- "Pseudoflavonifractor"
res$GS[res$GS =="Mycobacterium NA"] <- "Mycobacterium"
res$GS[res$GS =="Propionibacterium NA"] <- "Propionibacterium"
res$GS[res$GS =="Chryseobacterium NA"] <- "Chryseobacterium"
res$GS[res$GS =="Elizabethkingia NA"] <- "Elizabethkingia"
res$GS[res$GS =="Flavobacterium NA"] <- "Flavobacterium"
res$GS[res$GS =="Anoxybacillus NA"] <- "Anoxybacillus"
res$GS[res$GS =="Gemella NA"] <- "Gemella"
res$GS[res$GS =="Staphylococcus NA"] <- "Staphylococcus"
res$GS[res$GS =="Granulicatella NA"] <- "Granulicatella"
res$GS[res$GS =="Streptococcus NA"] <- "Streptococcus"
res$GS[res$GS =="Parvimonas NA"] <- "Parvimonas"
res$GS[res$GS =="Oribacterium NA"] <- "Oribacterium"
res$GS[res$GS =="Megamonas NA"] <- "Megamonas"
res$GS[res$GS =="Veillonella NA"] <- "Veillonella"
res$GS[res$GS =="Fusobacterium NA"] <- "Fusobacterium"
res$GS[res$GS =="Burkholderia NA"] <- "Burkholderia"
res$GS[res$GS =="Janthinobacterium NA"] <- "Janthinobacterium"
res$GS[res$GS =="Legionella NA"] <- "Legionella"
res$GS[res$GS =="Pseudomonas NA"] <- "Pseudomonas"

res$GS[res$GS =="Actinomyces Species_13071"] <- "Actinomyces sp."
res$GS[res$GS =="Actinomyces Species_13171"] <- "Actinomyces sp."
res$GS[res$GS =="Actinomyces Species_13303"] <- "Actinomyces sp."
res$GS[res$GS =="Actinomyces Species_13328"] <- "Actinomyces sp."
res$GS[res$GS =="Enterobacter Species_16784"] <- "Enterobacter sp."
res$GS[res$GS =="Bacteroides Species_19535"] <- "Bacteroides sp."
res$GS[res$GS =="Prevotella Species_20272"] <- "Prevotella sp."
res$GS[res$GS =="Prevotella Species_20317"] <- "Prevotella sp."
res$GS[res$GS =="Clostridium_sensu_stricto Species_23287"] <- "Clostridium_sensu_stricto sp."
res$GS[res$GS =="Romboutsia Species_23651"] <- "Romboutsia sp."
res$GS[res$GS =="Collinsella Species_11855"] <- "Collinsella sp."
res$GS[res$GS =="Parasutterella Species_15000"] <- "Parasutterella sp."
res$GS[res$GS =="Alistipes Species_18579"] <- "Alistipes sp."
res$GS[res$GS =="Barnesiella Species_20818"] <- "Barnesiella sp."
res$GS[res$GS =="Odoribacter Species_21279"] <- "Odoribacter sp."
res$GS[res$GS =="Dorea Species_26410"] <- "Dorea sp."
res$GS[res$GS =="Dorea Species_26435"] <- "Dorea sp."
res$GS[res$GS =="Ruminococcus2 Species_26537"] <- "Ruminococcus2 sp."
res$GS[res$GS =="Ruminococcus2 Species_26538"] <- "Ruminococcus2 sp."
res$GS[res$GS =="Ruminococcus2 Species_26546"] <- "Ruminococcus2 sp."
res$GS[res$GS =="Roseburia Species_27431"] <- "Roseburia sp."
res$GS[res$GS =="Coprococcus Species_3251"] <- "Coprococcus sp."
res$GS[res$GS =="Fusicatenibacter Species_389"] <- "Fusicatenibacter sp."
res$GS[res$GS =="Eisenbergiella Species_4174"] <- "Eisenbergiella sp."
res$GS[res$GS =="Faecalibacterium Species_4664"] <- "Faecalibacterium sp."
res$GS[res$GS =="Faecalibacterium Species_4902"] <- "Faecalibacterium sp."
res$GS[res$GS =="Faecalibacterium Species_4916" ] <- "Faecalibacterium sp."
res$GS[res$GS =="Faecalibacterium Species_4932"] <- "Faecalibacterium sp."
res$GS[res$GS =="Faecalibacterium Species_4991"] <- "Faecalibacterium sp."
res$GS[res$GS =="Oscillibacter Species_5838"] <- "Oscillibacter sp."
res$GS[res$GS =="Roseburia Species_59"] <- "Unknown Roseburia sp."
res$GS[res$GS =="Oscillibacter Species_6003"] <- "Oscillibactere sp."
res$GS[res$GS =="Oscillibacter Species_6024"] <- "Oscillibacter sp."
res$GS[res$GS =="Ruminococcus Species_7793"] <- "Ruminococcus sp."
res$GS[res$GS =="Ruminococcus Species_7801"] <- "Ruminococcus sp."
res$GS[res$GS =="Ruminococcus Species_7809"] <- "Ruminococcus sp."



sigTaxa.nonIBD<- c("Eubacterium limosum",
                   "Bacteroides",
                   "Anaerostipes",
                   "Blautia",
                   "Clostridium_XlVa",
                   "Clostridium_XlVb",
                   "Unknown Lachnospiraceae sp.",
                   "Unknown Erysipelotrichaceae sp.",
                   "Unknown Ruminococcaceae sp.",
                   "Anaerotruncus",
                   "Clostridium_IV",
                   "Flavonifractor",
                   "Pseudoflavonifractor",
                   "Collinsella sp.",
                   "Parasutterella sp.",
                   "Alistipes sp.",
                   "Barnesiella sp.",
                   "Odoribacter sp.",
                   "Dorea sp.",
                   "Ruminococcus2 sp.",
                   "Coprococcus sp.",
                   "Fusicatenibacter sp.",
                   "Eisenbergiella sp.",
                   "Faecalibacterium sp.",
                   "Roseburia sp.",
                   "Oscillibacter sp.",
                   "Ruminococcus sp.",
                   "Megasphaera micronuciformis",
                   "Mycobacterium",
                   "Propionibacterium",
                   "Chryseobacterium",
                   "Elizabethkingia",
                   "Flavobacterium",
                   "Anoxybacillus",
                   "Gemella",
                   "Staphylococcus",
                   "Granulicatella",
                   "Streptococcus",
                   "Parvimonas",
                   "Oribacterium",
                   "Megamonas",
                   "Veillonella",
                   "Fusobacterium",
                   "Burkholderia",
                   "Janthinobacterium",
                   "Legionella",
                   "Pseudomonas",
                   "Actinomyces sp.",
                   "Enterobacter sp.",
                   "Bacteroides sp.",
                   "Prevotella sp.",
                   "Clostridium_sensu_stricto sp.",
                   "Romboutsia sp.")

sigTaxa.nonIBD.2<- c("Bacteroides",
                     "Blautia",
                     "Clostridium_XlVa",
                     "Unknown Lachnospiraceae sp.",
                     "Faecalibacterium sp.",
                     "Roseburia sp.",
                     "Burkholderia")

res$log2FoldChange <- (-1)*res$log2FoldChange
nonIBD.da <- EnhancedVolcano(res,
                         lab = res$GS,
                         x = 'log2FoldChange',
                         y = 'global.pvalue',
                         selectLab = sigTaxa.nonIBD.2,
                         xlim = c(-6.5,6.5),
                         title = 'Alterations in the nonIBD Mucosa',
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
                         shapeCustom = TaxaVals.nonIBD.shape,
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
