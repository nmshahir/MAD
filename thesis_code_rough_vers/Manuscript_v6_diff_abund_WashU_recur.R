###########################################
# Figure 2 - Ileal Disease - Genera level #
###########################################
devtools::install_github('kevinblighe/EnhancedVolcano')
library(EnhancedVolcano)
res <- WashU_CD_recur_padj_mra_aero

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
OAn.WashU.taxa <- c("Alistipes",
                    "Bacteroides",
                    "Prevotella",
                    "Anaerovorax",
                    "Catabacter",
                    "Clostridium_XlVa",
                    "Genus_27125",
                    "Romboutsia",
                    "Roseburia",
                    "Dialister",
                    "Pyramidobacter",
                    "Barnesiella",
                    "Clostridium_XI",
                    "Coprobacter",
                    "Parabacteroides",
                    "Clostridium_sensu_stricto")

OA.WashU.taxa <- c("Sphingobium",
                   "Vulcaniibacterium")

FAn.WashU.taxa<- c("Propionibacterium",
                   "Brevundimonas",
                   "Gemella",
                   "Streptococcus",
                   "Citrobacter",
                   "Caulobacter",
                   "Enterococcus",
                   "Genus_16622",
                   "Haemophilus",
                   "Corynebacterium")

Aero.WashU.taxa <- c("Mycobacterium",
                     "Neisseria")

Anaero.WashU.taxa <- c("Sutterella",
                     "Odoribacter")

TaxaVals.WashU.shape <- rep(0,nrow(res))
names(TaxaVals.WashU.shape) <- rep('NS',nrow(res))

TaxaVals.WashU.shape[which(res$Genus %in% OAn.WashU.taxa)] <- 16
names(TaxaVals.WashU.shape)[which(res$Genus %in% OAn.WashU.taxa)] <- 'Obligate Anaerobe'

TaxaVals.WashU.shape[which(res$Genus %in% OA.WashU.taxa)] <- 2
names(TaxaVals.WashU.shape)[which(res$Genus %in% OA.WashU.taxa)] <- 'Obligate Aerobe'

TaxaVals.WashU.shape[which(res$Genus %in% FAn.WashU.taxa)] <- 18
names(TaxaVals.WashU.shape)[which(res$Genus %in% FAn.WashU.taxa)] <- 'Facultative Anaerobe'

TaxaVals.WashU.shape[which(res$Genus %in% Aero.WashU.taxa)] <- 4
names(TaxaVals.WashU.shape)[which(res$Genus %in% Aero.WashU.taxa)] <- 'Aerobe'

TaxaVals.WashU.shape[which(res$Genus %in% Anaero.WashU.taxa)] <- 15
names(TaxaVals.WashU.shape)[which(res$Genus %in% Anaero.WashU.taxa)] <- 'Anaerobe'

unique(names(TaxaVals.WashU.shape))
unique(TaxaVals.WashU.shape)
TaxaVals.WashU.shape[1:20]



#Rename

res$GS[res$GS == "Streptococcus Species_8647"] <- "Streptococcus sp."
res$GS[res$GS == "Bacteroides Species_19917"] <- "Bacteroides sp."
res$GS[res$GS == "Bacteroides Species_19531"] <- "Bacteroides sp."
res$GS[res$GS =="Bacteroides Species_19874"] <- "Bacteroides sp."
res$GS[res$GS =="Bacteroides Species_20723"	] <- "Bacteroides sp."
res$GS[res$GS =="Bacteroides Species_19625"	] <- "Bacteroides sp."
res$GS[res$GS =="Bacteroides Species_19555"	] <- "Bacteroides sp."
res$GS[res$GS =="Bacteroides Species_19435"	] <- "Bacteroides sp."

res$GS[res$GS == "Sutterella Species_14902"] <- "Sutterella sp."
res$GS[res$GS =="Clostridium_XlVa Species_27419"] <- "Clostridium_XlVa sp."
res$GS[res$GS =="Genus_27125 NA"] <- "Unknown Lachnospiraceae species"
res$GS[res$GS =="Roseburia Species_26687"] <- "Roseburia sp."



sigTaxa.WashU<- c("Streptococcus sp.",
                "Bacteroides sp.",
                "Roseburia intestinalis",
                "Roseburia sp.",
                "Clostridium_XlVa sp.",
                "Blautia sp.",
                "Sutterella sp.")

res$log2FoldChange <- (-1)*res$log2FoldChange
WashU.da <- EnhancedVolcano(res,
                          lab = res$GS,
                          x = 'log2FoldChange',
                          y = 'global.pvalue',
                          #selectLab = sigTaxa.WashU,
                          xlim = c(-6.5,6.5),
                          title = 'Alterations in the Recurrence - WashU',
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
                          shapeCustom = TaxaVals.WashU.shape,
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

WashU.da + coord_flip() 
