######################################################################
# abundance_comparison.R
# Goal: Implementation of Terry's 
# Bacterial High Abundance Tracker
# Input: Phyloseq object, string of subgroup A, string of subgroup B
# Step(s):
# 1 - Get all taxa that are > 5% relative abundance in group A
# 2 - For each taxa found in step 1, check if it is 
#    a) < 5% relative abundance in group B AND 
#    b) 3-fold less than the greatest relative abundance in group B
# 3 - Return list of taxa that fulfill criterion in step 2
#######################################################################

test <- readRDS("nonIBD_only_UNC.rds")
test <- readRDS("nonIBD_species.rds")
abundTrack <- function(phyla_object, grpA,grpB){
  
}

nPercent <- function(phyla_object, n){
  
}

find.top.taxa <- function(physeq,taxa){
  require(phyloseq)
  top.taxa <- transform_sample_counts(physeq,function(x) x/sum(x)) #tax_glom(x, taxa)
  otu <- otu_table(top.taxa) # gets the otu count table
  tax <- tax_table(top.taxa) # gets taxa table
  j<-apply(otu,1,which.max) # goes through each row and finds the OTU that is the max for that sample
  k <- j[!duplicated(j)] # removes duplicate OTUs
  l <- data.frame(tax@.Data[k,]) # gets data frame with taxa from
  m <- data.frame(otu@.Data[,k]) # gets otu table with only most abundant taxa
  #s <- as.name(taxa) #???
  colnames(m) = l[,taxa] #makes taxa level picked (i.e. phylum, genus, etc) the column name for m
  n <- colnames(m)[apply(m,1,which.max)] #makes new data frame with 
  m[,taxa] <- n
  return(m)
}
find.top.taxa2 <- function(physeq){
  require(phyloseq)
  top.taxa <- transform_sample_counts(physeq,function(x) x/sum(x)) #tax_glom(x, taxa)
  otu <- otu_table(top.taxa) # gets the otu count table
  tax <- tax_table(top.taxa) # gets taxa table
  j<-apply(otu,1,which.max) # goes through each row and finds the OTU that is the max for that sample
  k <- j[!duplicated(j)] # removes duplicate OTUs
  l <- data.frame(tax@.Data[k,]) # gets data frame with taxa from
  m <- data.frame(otu@.Data[,k]) # gets otu table with only most abundant taxa
  l$GS <- paste(l$Genus,l$Species)
  colnames(m) = l[,"GS"] #makes taxa level picked (i.e. phylum, genus, etc) the column name for m
  n <- colnames(m)[apply(m,1,which.max)] #makes new data frame with 
  m[,"GS"] <- n
  return(m)
}

find.top.taxa(test,"Genus")
find.top.taxa(test,"Species")

top.tax <- transform_sample_counts(test,function(x) x/sum(x))
top.tax.c <- subset_samples(top.tax,sample_data(top.tax)$AnatomSite == "colon")
otu <- otu_table(top.tax.c)
tax <- tax_table(top.tax.c)
j<-apply(otu@.Data,1,which.max)
k <- j[!duplicated(j)]
l <- data.frame(tax@.Data[k,])
m <- data.frame(otu@.Data[,k])
l$GS <- paste(l$Genus,l$Species)
colnames(m) = l[,"GS"]
n <- colnames(m)[apply(m,1,function(x) x > 0.1)]
m[,"GS"] <- n

##############################################################################
top.tax <- transform_sample_counts(test,function(x) x/sum(x))
top.tax.c <- subset_samples(top.tax,sample_data(top.tax)$AnatomSite == "colon")
otu <- otu_table(top.tax.c)
tax <- tax_table(top.tax.c)
j1<-apply(otu@.Data,1,function(x) which(x > 0.05))
j2 <-unlist(j1) # select for index
k1 <- j2[!duplicated(j2)]
l1 <- data.frame(tax@.Data[k1,])
m1 <- data.frame(otu@.Data[,k1])
l1$GS <- paste(l1$Genus,l1$Species)
colnames(m1) = l1[,"GS"]
###############################################################################

find.top.taxa2 <- function(x,taxa,num){
  require(phyloseq)
  require(magrittr)
  
  top.taxa <- tax_glom(x, taxa)
  otu <- otu_table(top.taxa) # remove the transformation if using a merge_sample object
  tax <- tax_table(top.taxa)
  j1 <- apply(otu,1,sort,index.return=T, decreasing=T) # modifying which.max to return a list of sorted index
  j2 <- lapply(j1,'[[',"ix") # select for index
  
  #k <- j[!duplicated(j)] # Replaced with unique() below
  l <- data.frame(unique(tax@.Data[unlist(j2),]))
  m <- data.frame(otu@.Data[,unique(unlist(j2))])
  #s <- as.name(taxa) # This isn't necessary!
  colnames(m) = l[,taxa]
  n <- apply(m,1,sort,index.return=T, decreasing=T) %>%
    lapply('[[',"ix") %>%  # Extract index
    lapply(head,n=num) # This to returns the top x tax
  # I want to apply a list of list of index to obtains a list of list of translated names
  
  
  # https://stackoverflow.com/questions/31561238/lapply-function-loops-on-list-of-lists-r
  p <- list()
  for(i in 1:length(n)){
    p[[i]]<- colnames(m)[n[[i]]]
  }
  m$taxa <- p # replacing [,taxa], but now the new column is called "taxa" instead of the inputted taxonomic rank
  return(m)
}

###############################################################################
UC.match <- subset_samples(UC.patients, sample_data(UC.patients)$PatientNo %in% c("52","46","43800","43500","43000","42700","42300","41200","41000","40600","40400","40"))
UC.match <- prune_taxa(taxa_sums(UC.match)>0,UC.match)
