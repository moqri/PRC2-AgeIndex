
########## Using AnnotatR to annotate the top1000 LMRs 01/13/23 ##########

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("annotatr")

BiocManager::install("GenomicFeatures")

BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

BiocManager::install("org.Hs.eg.db") 

BiocManager::install("rtracklayer")


library(annotatr)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)
library(GenomicRanges)

txdb.hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene


######### Loading Blood top 1000 LMRs #########

BloodLMRs1000 = read.table('Blood_top1000_LMRs_Redit.bed')



#turning the bed file into GR ranges object
BloodLMRs1000_gr <- GRanges(seqnames = BloodLMRs1000$V1,
                            ranges = IRanges(start = BloodLMRs1000$V2, end = BloodLMRs1000$V4),
                            strand = BloodLMRs1000$V3,
                            pr = BloodLMRs1000$V6,
                            delta = BloodLMRs1000$V7,
                            number = seq(1:length(BloodLMRs1000$V1))) #corresponds to row number 


#how many hypermethylate with age

length(which(BloodLMRs1000_gr$delta > 0))/1000
#0.945 - 95%

############# Histogram of LMR sizes, blood #############

pdf(file = "BloodLMRs_histo.pdf",width = 8,height = 6, useDingbats=FALSE)

hist(as.data.frame(BloodLMRs1000_gr)$width, breaks = 25, col = "#C03F3F", border= "white",
     xlab = "LMR Size (bp)", ylab = "No. of LMRs",main = "Blood LMR Size Distribution")

hist(log10(as.data.frame(BloodLMRs1000_gr)$width), breaks = 25, col = "#C03F3F", border= "white",
     xlab = "Log10 LMR Size (bp)", ylab = "No. of LMRs", main = "Blood LMR Size Distribution")

dev.off()


hg38_annos <- builtin_annotations()[grep("hg38", builtin_annotations())]

# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg38', annotations = hg38_annos)

#now to annotate our LMRs.
BloodLMRs_annotated = annotate_regions(
  regions = BloodLMRs1000_gr,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)


#Turning to a dataframe to subset and use more easily
BloodLMRs_annotated.df <- data.frame(BloodLMRs_annotated)



#Filtering out some of info 

table(BloodLMRs_annotated.df$annot.type)

delsies <- c("hg38_genes_1to5kb","hg38_genes_exonintronboundaries","hg38_genes_intronexonboundaries", "hg38_genes_cds")

BloodLMRs_annotated.df_filt <- BloodLMRs_annotated.df[!BloodLMRs_annotated.df$annot.type %in% delsies,]

#cleaning up the terms
BloodLMRs_annotated.df_filt$annot.type <- gsub("hg38_", "", BloodLMRs_annotated.df_filt$annot.type)

BloodLMRs_annotated.df_filt$annot.type <- gsub("genes_", "", BloodLMRs_annotated.df_filt$annot.type)


###### Making blood annot table final ###### 


#Summarising by category to then sum into a table
BldLMR_Catsum_filt = summarize_categorical(
  annotated_regions = BloodLMRs_annotated.df_filt,
  by =  c("number", "annot.type", "annot.symbol"),
  quiet = TRUE)
print(BldLMR_Catsum_filt)


BldLMR_Catsum_filt <- as.data.frame(BldLMR_Catsum_filt)

#clean annotation table, easier to use gr formatted version
BloodLMRs1000_cleanAnno <- as.data.frame(BloodLMRs1000_gr)

#Summarizing genes and features into final table

BloodLMR_Summed_annos <- list()
BloodLMR_Summed_genes <- list()
for(r in 1:length(rownames(BloodLMRs1000_cleanAnno))){
  print(r)
  i <- rownames(BloodLMRs1000_cleanAnno)[r]
  print(i)
  temp <- BldLMR_Catsum_filt[BldLMR_Catsum_filt$number == i,]
  
  stuff <- toString(unique(na.omit(temp$annot.type))) #merges all into one entry, useful
  
  BloodLMR_Summed_annos[[r]] <- c(i,toString(stuff))
  
  stuff2 <- toString(unique(na.omit(temp$annot.symbol))) #merges all into one entry, useful
  
  BloodLMR_Summed_genes[[r]] <- c(i,toString(stuff2))
}


BloodLMR_Summed_genes <-as.data.frame(t(as.data.frame(BloodLMR_Summed_genes)))



BloodLMR_Summed_annos <-as.data.frame(t(as.data.frame(BloodLMR_Summed_annos)))


#since its the same order and number, should just be fine to append on
BloodLMRs1000_cleanAnno$Annotated.gene.symbols <-  BloodLMR_Summed_genes$V2



#since its the same order and number, should just be fine to append on
BloodLMRs1000_cleanAnno$Detected.Features <-  BloodLMR_Summed_annos$V2


############## Writing out clean Bld anno table ##############

write.table(BloodLMRs1000_cleanAnno,file = "Blood_Top1000LRMs_finalTable.txt",sep = "\t", quote = F, row.names = F)



############## Removing absolute duplicates within each category  ##############

# i.e. removing annotations with exactly the same coordiantes within each category, such as multiple repeat exons

#extracting just the annotations, so removing the LMR annotations for now
BloodLMRs_annotonly_Df <- BloodLMRs_annotated.df_filt[,-(1:8)]


#making unified names to make the matching easier
BloodLMRs_annotonly_Df$coordname <- paste(BloodLMRs_annotonly_Df$annot.seqnames,(paste(BloodLMRs_annotonly_Df$annot.start, BloodLMRs_annotonly_Df$annot.end)))

BldLMR_Unique_annos <- list()
for(i in levels(as.factor(BloodLMRs_annotonly_Df$annot.type))){
  anno <- BloodLMRs_annotonly_Df[BloodLMRs_annotonly_Df$annot.type == i,]
  anno <- anno[!duplicated(anno$coordname),]
  BldLMR_Unique_annos[[i]] <- anno
  
}


BldLMR_Unique_annos <- as.data.frame(do.call(rbind, BldLMR_Unique_annos))


#removing first exons present in the "exon" category

#extracting from the exon coordinate names the ones that are not present in the first introns
uniqueExons <- BldLMR_Unique_annos$coordname[BldLMR_Unique_annos$annot.type == "exons"][which(!BldLMR_Unique_annos$coordname[BldLMR_Unique_annos$annot.type == "exons"] %in% BldLMR_Unique_annos$coordname[BldLMR_Unique_annos$annot.type == "firstexons"])]

#creating a separate exons table to edit and get the other exons
BldLMR_Unique_AllExons <- BldLMR_Unique_annos[BldLMR_Unique_annos$annot.type == "exons",]

BldLMR_Unique_AllExons <- BldLMR_Unique_AllExons[BldLMR_Unique_AllExons$coordname %in% uniqueExons,]

BldLMR_Unique_AllExons$annot.type <- rep("other.exons", length(BldLMR_Unique_AllExons$annot.type))


#removing exons from main table and adding in new one

BldLMR_Unique_Annos_ex <- BldLMR_Unique_annos[!BldLMR_Unique_annos$annot.type == "exons",]

BldLMR_Unique_Annos_ex <- rbind(BldLMR_Unique_Annos_ex, BldLMR_Unique_AllExons)

table(BldLMR_Unique_Annos_ex$annot.type)


#Okay looks good, now moving onto the graph part, need to do proportions first
library(reshape2)
library(dplyr)
library(tidyr)

#filtering just for bp and annotation type
BldLMR_Unique_Annos_ex_fiddle <- BldLMR_Unique_Annos_ex[,c("annot.width","annot.type")]

BldLMR_Unique_Annos_ex_fiddle <- BldLMR_Unique_Annos_ex_fiddle %>%
  pivot_wider(names_from = annot.type, values_from = annot.width, values_fn = sum)


BldLMR_Unique_Annos_ex_fiddle <- melt(BldLMR_Unique_Annos_ex_fiddle)


####### Sorting Blood terms for plotting later ####### 

#Creating plots with proportion relative to itself, Ie divide total bp of each feature by total bp of all LMR annotation features 
#(rather than total bp of all LMRs) that should be relatively fair, especially if I'm comparing the same thing with the genome.

#sort into categories

#cpgs only
Bld_LMR_CpGanno <- BldLMR_Unique_Annos_ex_fiddle[3:6,]
#genomic only
Bld_LMR_Genomic.anno <- BldLMR_Unique_Annos_ex_fiddle[-(3:6),]


#dividing 
Bld_LMR_CpGanno$Proportion <- (Bld_LMR_CpGanno$value/sum(Bld_LMR_CpGanno$value))*100
Bld_LMR_CpGanno$Group <- rep("Blood LMR CpGs", length(Bld_LMR_CpGanno$variable))

Bld_LMR_Genomic.anno$Proportion <- (Bld_LMR_Genomic.anno$value/sum(Bld_LMR_Genomic.anno$value))*100
Bld_LMR_Genomic.anno$Group <- rep("Blood LMR Genomic", length(Bld_LMR_Genomic.anno$variable))





###### Extracting blood gene features for GO analysis (Prob GREAT...maybe webgestalt but i doubt it) #############


#GREAT is used for annotating cis-acting regulatory regions (cis-acting non-coding DNA regions that regulate the transcription of genes eg. Promoters, enhancers, silencers, and insulators),usually found by chip-seq experiments. 

#Extracting a list of first exons, promoters and enhancers.

Blood_CREs <- BldLMR_Unique_Annos_ex

#need specific columns only for GREAT bed; chr, start, end, name

Blood_CREs <- Blood_CREs[,1:3]

Blood_CREs$names <- BldLMR_Unique_Annos_ex$annot.tx_id

Blood_CREs <- Blood_CREs[BldLMR_Unique_Annos_ex$annot.type == "promoters" | BldLMR_Unique_Annos_ex$annot.type == "enhancers_fantom" | BldLMR_Unique_Annos_ex$annot.type == "firstexons",]

write.table(Blood_CREs, "Blood_LMR_CisRegs_Fexons_filt.bed",quote = F,row.names = F)



############################################################################################


######### Loading skin top 1000 LMRs #########

SkinLMRs1000 = read.table('Skin_top1000_LMRs_Rver.bed')

#turning the bed file into GR ranges object
SkinLMRs1000_gr <- GRanges(seqnames = SkinLMRs1000$V1,
                           ranges = IRanges(start = SkinLMRs1000$V2, end = SkinLMRs1000$V4),
                           strand = SkinLMRs1000$V3,
                           pr = SkinLMRs1000$V6,
                           delta = SkinLMRs1000$V7,
                           number = seq(1:length(SkinLMRs1000$V1))) #corresponds to row number 


#How many top 1000 LMRs hypermethylate with age
length(which(SkinLMRs1000_gr$delta > 0))/1000
#0.83 - 83%


############# Histogram of LMR sizes, skin #############

pdf(file = "SkinLMRs_histo.pdf",width = 8,height = 6, useDingbats=FALSE)

hist(as.data.frame(SkinLMRs1000_gr)$width, breaks = 25, col = "#3671A8",border= "white", 
     xlab = "LMR Size (bp)", ylab = "No. of LMRs",main = "Epidermis LMR Size Distribution")

hist(log10(as.data.frame(SkinLMRs1000_gr)$width), breaks = 25, col = "#3671A8",border= "white", 
     xlab = "Log10 LMR Size (bp)", ylab = "No. of LMRs", main = "Epidermis LMR Size Distribution")

dev.off()



############# Annotating Skin LMRS #############

#now to annotate our LMRs.
SkinLMRs_annotated = annotate_regions(
  regions = SkinLMRs1000_gr,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

#Turning to a dataframe to subset and use more easily
SkinLMRs_annotated.df <- data.frame(SkinLMRs_annotated)


#Filtering out some of info that doesn't seem as useful

table(SkinLMRs_annotated.df$annot.type)

delsies <- c("hg38_genes_1to5kb","hg38_genes_exonintronboundaries","hg38_genes_intronexonboundaries", "hg38_genes_cds")

SkinLMRs_annotated.df_filt <- SkinLMRs_annotated.df[!SkinLMRs_annotated.df$annot.type %in% delsies,]

#cleaning up the terms
SkinLMRs_annotated.df_filt$annot.type <- gsub("hg38_", "", SkinLMRs_annotated.df_filt$annot.type)

SkinLMRs_annotated.df_filt$annot.type <- gsub("genes_", "", SkinLMRs_annotated.df_filt$annot.type)


############ Making Skin Annot table final  ############

#Summarising by category to then sum into a table
SkinLMR_Catsum_filt = summarize_categorical(
  annotated_regions = SkinLMRs_annotated.df_filt,
  by =  c("number", "annot.type", "annot.symbol"),
  quiet = TRUE)
print(SkinLMR_Catsum_filt)


SkinLMR_Catsum_filt <- as.data.frame(SkinLMR_Catsum_filt)



##### clean Skin annotation table ###### 
#easier to use gr formatted version
SkinLMRs1000_cleanAnno <- as.data.frame(SkinLMRs1000_gr)



SkinLMR_Summed_annos <- list()
SkinLMR_Summed_genes <- list()
for(r in 1:length(rownames(SkinLMRs1000_cleanAnno))){
  print(r)
  i <- rownames(SkinLMRs1000_cleanAnno)[r]
  print(i)
  temp <- SkinLMR_Catsum_filt[SkinLMR_Catsum_filt$number == i,]
  
  stuff <- toString(unique(na.omit(temp$annot.type))) #merges all into one entry, useful
  
  SkinLMR_Summed_annos[[r]] <- c(i,toString(stuff))
  
  stuff2 <- toString(unique(na.omit(temp$annot.symbol))) #merges all into one entry, useful
  
  SkinLMR_Summed_genes[[r]] <- c(i,toString(stuff2))
}

SkinLMR_Summed_genes <-as.data.frame(t(as.data.frame(SkinLMR_Summed_genes)))

SkinLMR_Summed_annos <-as.data.frame(t(as.data.frame(SkinLMR_Summed_annos)))


#since its the same order and number, should just be fine to append on
SkinLMRs1000_cleanAnno$Annotated.gene.symbols <-  SkinLMR_Summed_genes$V2



#since its the same order and number, should just be fine to append on
SkinLMRs1000_cleanAnno$Detected.Features <-  SkinLMR_Summed_annos$V2


#### Writing out clean skin anno table ##############

write.table(SkinLMRs1000_cleanAnno,file = "Skin_Top1000LRMs_finalTable.txt",sep = "\t", quote = F, row.names = F)

############## Removing absolute duplicates within each category  ##############

# i.e. removing annotations with exactly the same coordiantes within each category, such as multiple repeat exons

#extracting just the annotations, so removing the LMR annotations for now
SkinLMRs_annotonly_Df <- SkinLMRs_annotated.df_filt[,-(1:8)]

#making unified names to make the matching easier
SkinLMRs_annotonly_Df$coordname <- paste(SkinLMRs_annotonly_Df$annot.seqnames,(paste(SkinLMRs_annotonly_Df$annot.start, SkinLMRs_annotonly_Df$annot.end)))

SkinLMR_Unique_annos <- list()
for(i in levels(as.factor(SkinLMRs_annotonly_Df$annot.type))){
  anno <- SkinLMRs_annotonly_Df[SkinLMRs_annotonly_Df$annot.type == i,]
  anno <- anno[!duplicated(anno$coordname),]
  SkinLMR_Unique_annos[[i]] <- anno
  
}


SkinLMR_Unique_annos <- as.data.frame(do.call(rbind, SkinLMR_Unique_annos))


#removing first exons from "exons" and renaming to "other.exons"

#extracting from the exon coordinate names the ones that are not present in the first introns
uniqueExons <- SkinLMR_Unique_annos$coordname[SkinLMR_Unique_annos$annot.type == "exons"][which(!SkinLMR_Unique_annos$coordname[SkinLMR_Unique_annos$annot.type == "exons"] %in% SkinLMR_Unique_annos$coordname[SkinLMR_Unique_annos$annot.type == "firstexons"])]

#creating a separate exons table to edit and get the other exons
SkinLMR_Unique_AllExons <- SkinLMR_Unique_annos[SkinLMR_Unique_annos$annot.type == "exons",]

SkinLMR_Unique_AllExons <- SkinLMR_Unique_AllExons[SkinLMR_Unique_AllExons$coordname %in% uniqueExons,]

SkinLMR_Unique_AllExons$annot.type <- rep("other.exons", length(SkinLMR_Unique_AllExons$annot.type))


#removing exons from main table and adding in new one

SkinLMR_Unique_Annos_ex <- SkinLMR_Unique_annos[!SkinLMR_Unique_annos$annot.type == "exons",]

SkinLMR_Unique_Annos_ex <- rbind(SkinLMR_Unique_Annos_ex, SkinLMR_Unique_AllExons)

table(SkinLMR_Unique_Annos_ex$annot.type)


####### Sorting Blood terms for plotting later ####### 

library(reshape2)
library(dplyr)
library(tidyr)

#filtering just for bp and annotation type
SkinLMR_Unique_Annos_ex_fiddle <- SkinLMR_Unique_Annos_ex[,c("annot.width","annot.type")]

SkinLMR_Unique_Annos_ex_fiddle <- SkinLMR_Unique_Annos_ex_fiddle %>%
  pivot_wider(names_from = annot.type, values_from = annot.width, values_fn = sum)


SkinLMR_Unique_Annos_ex_fiddle <- melt(SkinLMR_Unique_Annos_ex_fiddle)

#sort into categories

#CpG only
Skin_LMR_CpGanno <- SkinLMR_Unique_Annos_ex_fiddle[3:6,]
#genomic only
Skin_LMR_Genomic.anno <- SkinLMR_Unique_Annos_ex_fiddle[-(3:6),]


#dividing 
Skin_LMR_CpGanno$Proportion <- (Skin_LMR_CpGanno$value/sum(Skin_LMR_CpGanno$value))*100
Skin_LMR_CpGanno$Group <- rep("Skin LMR CpGs", length(Skin_LMR_CpGanno$variable))

Skin_LMR_Genomic.anno$Proportion <- (Skin_LMR_Genomic.anno$value/sum(Skin_LMR_Genomic.anno$value))*100
Skin_LMR_Genomic.anno$Group <- rep("Skin LMR Genomic", length(Skin_LMR_Genomic.anno$variable))




#### writing out skin features for GREAT ####
Skin_CREs <- SkinLMR_Unique_Annos_ex


Skin_CREs <- Skin_CREs[,1:3]

Skin_CREs$names <- SkinLMR_Unique_Annos_ex$annot.tx_id

Skin_CREs <- Skin_CREs[SkinLMR_Unique_Annos_ex$annot.type == "promoters" | SkinLMR_Unique_Annos_ex$annot.type == "enhancers_fantom" | SkinLMR_Unique_Annos_ex$annot.type == "firstexons",]



write.table(Skin_CREs, "Skin_LMR_CisRegs_Fexons_filt.bed",quote = F,row.names = F)




############################################################################################


######### Loading Fibro top 1000 LMRs #########

#reading in Fibro (fv) LMRs
FibroLMRs1000 = read.table('Fibro_fv_Top1000.bed')

#turning the bed file into GR ranges object
FibroLMRs1000_gr <- GRanges(seqnames = FibroLMRs1000$V1,
                            ranges = IRanges(start = FibroLMRs1000$V2, end = FibroLMRs1000$V3),
                            strand = rep("*", length(FibroLMRs1000$V1)),
                            pr = FibroLMRs1000$V5,
                            delta = FibroLMRs1000$V6,
                            number = seq(1:length(FibroLMRs1000$V1))) #corresponds to row number 


###Histogram of LMR sizes, Fibro #############

pdf(file = "FibroLMRs_histo.pdf",width = 8,height = 6, useDingbats=FALSE)

hist(as.data.frame(FibroLMRs1000_gr)$width, breaks = 25, col = "#E2812C", border= "white",
     xlab = "LMR Size (bp)", ylab = "No. of LMRs",main = "Fibroblast LMR Size Distribution")

hist(log10(as.data.frame(FibroLMRs1000_gr)$width), breaks = 25, col = "#E2812C",border= "white",  xlim = c(2.5,5),
     xlab = "Log10 LMR Size (bp)", ylab = "No. of LMRs", main = "Fibroblast LMR Size Distribution")

dev.off()



################### Annotating Fibro LMRS ###################

#Note - Fibro LMRs used here are using the PRC2 LMR regions id'd in hESCs, same as Blood and Epidermis/Fibro in the rest of the paper

FibroLMRs_annotated = annotate_regions(
  regions = FibroLMRs1000_gr,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)


FibroLMRs_annotated.df <- data.frame(FibroLMRs_annotated)

#Filtering out some of info that doesn't seem as useful

table(FibroLMRs_annotated.df$annot.type)

delsies <- c("hg38_genes_1to5kb","hg38_genes_exonintronboundaries","hg38_genes_intronexonboundaries", "hg38_genes_cds")

FibroLMRs_annotated.df_filt <- FibroLMRs_annotated.df[!FibroLMRs_annotated.df$annot.type %in% delsies,]

#cleaning up the terms
FibroLMRs_annotated.df_filt$annot.type <- gsub("hg38_", "", FibroLMRs_annotated.df_filt$annot.type)

FibroLMRs_annotated.df_filt$annot.type <- gsub("genes_", "", FibroLMRs_annotated.df_filt$annot.type)



############ Making Fibro Annot table final  ############

#Summarising by category to then sum into a table
FibroLMR_Catsum_filt = summarize_categorical(
  annotated_regions = FibroLMRs_annotated.df_filt,
  by =  c("number", "annot.type", "annot.symbol"),
  quiet = TRUE)
print(FibroLMR_Catsum_filt)


FibroLMR_Catsum_filt <- as.data.frame(FibroLMR_Catsum_filt)



######  clean Fibro annotation table ###### 

FibroLMRs1000_cleanAnno <- as.data.frame(FibroLMRs1000_gr)

FibroLMR_Summed_annos <- list()
FibroLMR_Summed_genes <- list()
for(r in 1:length(rownames(FibroLMRs1000_cleanAnno))){
  print(r)
  i <- rownames(FibroLMRs1000_cleanAnno)[r]
  print(i)
  temp <- FibroLMR_Catsum_filt[FibroLMR_Catsum_filt$number == i,]
  
  stuff <- toString(unique(na.omit(temp$annot.type))) #merges all into one entry, useful
  
  FibroLMR_Summed_annos[[r]] <- c(i,toString(stuff))
  
  stuff2 <- toString(unique(na.omit(temp$annot.symbol))) #merges all into one entry, useful
  
  FibroLMR_Summed_genes[[r]] <- c(i,toString(stuff2))
}


FibroLMR_Summed_genes <-as.data.frame(t(as.data.frame(FibroLMR_Summed_genes)))

FibroLMR_Summed_annos <-as.data.frame(t(as.data.frame(FibroLMR_Summed_annos)))


#since its the same order and number, should just be fine to append on
FibroLMRs1000_cleanAnno$Annotated.gene.symbols <-  FibroLMR_Summed_genes$V2



#since its the same order and number, should just be fine to append on
FibroLMRs1000_cleanAnno$Detected.Features <-  FibroLMR_Summed_annos$V2


############## Writing out clean Fibro anno table ##############

write.table(FibroLMRs1000_cleanAnno,file = "Fibro_Top1000LRMs_finalTable.txt",sep = "\t", quote = F, row.names = F)



############## Uniquing regions for Plotting later Fibro ##############

#Isolating annotations only of LMRs
FibroLMRs_annotonly_Df <- FibroLMRs_annotated.df_filt[,-(1:8)]

#making unified names to make the matching easier
FibroLMRs_annotonly_Df$coordname <- paste(FibroLMRs_annotonly_Df$annot.seqnames,(paste(FibroLMRs_annotonly_Df$annot.start, FibroLMRs_annotonly_Df$annot.end)))


FibroLMR_Unique_annos <- list()
for(i in levels(as.factor(FibroLMRs_annotonly_Df$annot.type))){
  anno <- FibroLMRs_annotonly_Df[FibroLMRs_annotonly_Df$annot.type == i,]
  anno <- anno[!duplicated(anno$coordname),]
  FibroLMR_Unique_annos[[i]] <- anno
  
}


FibroLMR_Unique_annos <- as.data.frame(do.call(rbind, FibroLMR_Unique_annos))


table(FibroLMR_Unique_annos$annot.type)


#removing duplicate exons 

#extracting from the exon coordinate names the ones that are not present in the first introns
uniqueExons <- FibroLMR_Unique_annos$coordname[FibroLMR_Unique_annos$annot.type == "exons"][which(!FibroLMR_Unique_annos$coordname[FibroLMR_Unique_annos$annot.type == "exons"] %in% FibroLMR_Unique_annos$coordname[FibroLMR_Unique_annos$annot.type == "firstexons"])]

#creating a separate exons table to edit and get the other exons
FibroLMR_Unique_AllExons <- FibroLMR_Unique_annos[FibroLMR_Unique_annos$annot.type == "exons",]

FibroLMR_Unique_AllExons <- FibroLMR_Unique_AllExons[FibroLMR_Unique_AllExons$coordname %in% uniqueExons,]

FibroLMR_Unique_AllExons$annot.type <- rep("other.exons", length(FibroLMR_Unique_AllExons$annot.type))


#removing exons from main table and adding in new one

FibroLMR_Unique_Annos_ex <- FibroLMR_Unique_annos[!FibroLMR_Unique_annos$annot.type == "exons",]

FibroLMR_Unique_Annos_ex <- rbind(FibroLMR_Unique_Annos_ex, FibroLMR_Unique_AllExons)

table(FibroLMR_Unique_Annos_ex$annot.type)





####writing out Fibro features for GREAT ####
Fibro_CREs <- FibroLMR_Unique_Annos_ex


Fibro_CREs <- Fibro_CREs[,1:3]

Fibro_CREs$names <- FibroLMR_Unique_Annos_ex$annot.tx_id


Fibro_CREs <- Fibro_CREs[FibroLMR_Unique_Annos_ex$annot.type == "promoters" | FibroLMR_Unique_Annos_ex$annot.type == "enhancers_fantom" | FibroLMR_Unique_Annos_ex$annot.type == "firstexons" ,]



write.table(Fibro_CREs, "Fibro_LMR_CisRegs_Fexon_filt.bed",quote = F,row.names = F)


##### Preparing Fibro data to plot #####

#filtering just for bp and annotation type
FibroLMR_Unique_Annos_ex_fiddle <- FibroLMR_Unique_Annos_ex[,c("annot.width","annot.type")]

FibroLMR_Unique_Annos_ex_fiddle <- FibroLMR_Unique_Annos_ex_fiddle %>%
  pivot_wider(names_from = annot.type, values_from = annot.width, values_fn = sum)


FibroLMR_Unique_Annos_ex_fiddle <- melt(FibroLMR_Unique_Annos_ex_fiddle)

#sort into categories

#cpg only
Fibro_LMR_CpGanno <- FibroLMR_Unique_Annos_ex_fiddle[3:6,]
#genomic only
Fibro_LMR_Genomic.anno <- FibroLMR_Unique_Annos_ex_fiddle[-(3:6),]


#dividing 
Fibro_LMR_CpGanno$Proportion <- (Fibro_LMR_CpGanno$value/sum(Fibro_LMR_CpGanno$value))*100
Fibro_LMR_CpGanno$Group <- rep("Fibro LMR CpGs", length(Fibro_LMR_CpGanno$variable))

Fibro_LMR_Genomic.anno$Proportion <- (Fibro_LMR_Genomic.anno$value/sum(Fibro_LMR_Genomic.anno$value))*100
Fibro_LMR_Genomic.anno$Group <- rep("Fibro LMR Genomic", length(Fibro_LMR_Genomic.anno$variable))




################ using all Blood LMRs as background for GREAT and plotting ################ 


#reading in blood LMRs
AllBldLMRs = read.table('Blood_LRM_CoordinatesOnly.bed')

#turning the bed file into GR ranges object
AllBldLMRs_gr <- GRanges(seqnames = AllBldLMRs$V1,
                         ranges = IRanges(start = AllBldLMRs$V2, end = AllBldLMRs$V3),
                         strand = rep("*", length(AllBldLMRs$V1)),
                         number = seq(1:length(AllBldLMRs$V1))) #corresponds to row number 





#now to annotate our LMRs.
AllBldLMRs_annotated = annotate_regions(
  regions = AllBldLMRs_gr,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)


AllBldLMRs_annotated.df <- data.frame(AllBldLMRs_annotated)

#removing unwanted terms as before
AllBldLMRs_annotated.df_filt <- AllBldLMRs_annotated.df[!AllBldLMRs_annotated.df$annot.type %in% delsies,]

#cleaning up the terms
AllBldLMRs_annotated.df_filt$annot.type <- gsub("hg38_", "", AllBldLMRs_annotated.df_filt$annot.type)

AllBldLMRs_annotated.df_filt$annot.type <- gsub("genes_", "", AllBldLMRs_annotated.df_filt$annot.type)


#isolating annotations only
AllBldLMRs_annotonly_Df <- AllBldLMRs_annotated.df_filt[,-(1:6)]

#making unified names to make the matching easier
AllBldLMRs_annotonly_Df$coordname <- paste(AllBldLMRs_annotonly_Df$annot.seqnames,(paste(AllBldLMRs_annotonly_Df$annot.start, AllBldLMRs_annotonly_Df$annot.end)))


AllBldLMR_Unique_annos <- list()
for(i in levels(as.factor(AllBldLMRs_annotonly_Df$annot.type))){
  anno <- AllBldLMRs_annotonly_Df[AllBldLMRs_annotonly_Df$annot.type == i,]
  anno <- anno[!duplicated(anno$coordname),]
  AllBldLMR_Unique_annos[[i]] <- anno
  
}


AllBldLMR_Unique_annos <- as.data.frame(do.call(rbind, AllBldLMR_Unique_annos))

#removing duplicate exons 

#extracting from the exon coordinate names the ones that are not present in the first introns
uniqueExons <- AllBldLMR_Unique_annos$coordname[AllBldLMR_Unique_annos$annot.type == "exons"][which(!AllBldLMR_Unique_annos$coordname[AllBldLMR_Unique_annos$annot.type == "exons"] %in% AllBldLMR_Unique_annos$coordname[AllBldLMR_Unique_annos$annot.type == "firstexons"])]

#creating a separate exons table to edit and get the other exons
BldLMR_Unique_AllExons <- AllBldLMR_Unique_annos[AllBldLMR_Unique_annos$annot.type == "exons",]

BldLMR_Unique_AllExons <- BldLMR_Unique_AllExons[BldLMR_Unique_AllExons$coordname %in% uniqueExons,]

BldLMR_Unique_AllExons$annot.type <- rep("other.exons", length(BldLMR_Unique_AllExons$annot.type))


#removing exons from main table and adding in new one

AllBldLMR_Unique_annos_ex <- AllBldLMR_Unique_annos[!AllBldLMR_Unique_annos$annot.type == "exons",]

AllBldLMR_Unique_annos_ex <- rbind(AllBldLMR_Unique_annos_ex, BldLMR_Unique_AllExons)

table(AllBldLMR_Unique_annos_ex$annot.type)


####writing out AllBld features for GREAT background use ####

#(GREAT requires the input samples to also be in the background, so the background needs to be all the annos of the LMRs, not just all LMR coords)
AllBld_CREs <- AllBldLMR_Unique_annos_ex

AllBld_CREs <- AllBld_CREs[,1:3]

AllBld_CREs$names <- AllBldLMR_Unique_annos_ex$annot.tx_id

AllBld_CREs <- AllBld_CREs[AllBldLMR_Unique_annos_ex$annot.type == "promoters" | AllBldLMR_Unique_annos_ex$annot.type == "enhancers_fantom" | AllBldLMR_Unique_annos_ex$annot.type == "firstexons",]




write.table(AllBld_CREs, "AllBLood_LMR_CisRegs_Fexons_filt.bed",quote = F,row.names = F)



#finding proportions of genomic annotations in all Blood LMRs

#filtering just for bp and annotation type
AllBldLMR_Unique_Annos_ex_fiddle <- AllBldLMR_Unique_annos_ex[,c("annot.width","annot.type")]

AllBldLMR_Unique_Annos_ex_fiddle <- AllBldLMR_Unique_Annos_ex_fiddle %>%
  pivot_wider(names_from = annot.type, values_from = annot.width, values_fn = sum)


AllBldLMR_Unique_Annos_ex_fiddle <- melt(AllBldLMR_Unique_Annos_ex_fiddle)

#Separating 

#Cpg only
AllBld_LMR_CpGanno <- AllBldLMR_Unique_Annos_ex_fiddle[3:6,]

#Genomic only
AllBld_LMR_Genomic.anno <- AllBldLMR_Unique_Annos_ex_fiddle[-(3:6),]


#dividing 
AllBld_LMR_CpGanno$Proportion <- (AllBld_LMR_CpGanno$value/sum(AllBld_LMR_CpGanno$value))*100
AllBld_LMR_CpGanno$Group <- rep("Blood LMR CpGs", length(AllBld_LMR_CpGanno$variable))

AllBld_LMR_Genomic.anno$Proportion <- (AllBld_LMR_Genomic.anno$value/sum(AllBld_LMR_Genomic.anno$value))*100
AllBld_LMR_Genomic.anno$Group <- rep("Blood LMR Genomic", length(AllBld_LMR_Genomic.anno$variable))


################ using all Skin LMRs as background for GREAT and plotting ################ 



#reading in Skin LMRs
AllSkinLMRs = read.table('skin_epid_LMRcoord_only.bed')

#turning the bed file into GR ranges object
AllSkinLMRs_gr <- GRanges(seqnames = AllSkinLMRs$V1,
                          ranges = IRanges(start = AllSkinLMRs$V2, end = AllSkinLMRs$V3),
                          strand = rep("*", length(AllSkinLMRs$V1)),
                          number = seq(1:length(AllSkinLMRs$V1))) #corresponds to row number 





#now to annotate our LMRs.
AllSkinLMRs_annotated = annotate_regions(
  regions = AllSkinLMRs_gr,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)


AllSkinLMRs_annotated.df <- data.frame(AllSkinLMRs_annotated)


AllSkinLMRs_annotated.df_filt <- AllSkinLMRs_annotated.df[!AllSkinLMRs_annotated.df$annot.type %in% delsies,]

#cleaning up the terms
AllSkinLMRs_annotated.df_filt$annot.type <- gsub("hg38_", "", AllSkinLMRs_annotated.df_filt$annot.type)

AllSkinLMRs_annotated.df_filt$annot.type <- gsub("genes_", "", AllSkinLMRs_annotated.df_filt$annot.type)



AllSkinLMRs_annotonly_Df <- AllSkinLMRs_annotated.df_filt[,-(1:6)]

#making unified names to make the matching easier
AllSkinLMRs_annotonly_Df$coordname <- paste(AllSkinLMRs_annotonly_Df$annot.seqnames,(paste(AllSkinLMRs_annotonly_Df$annot.start, AllSkinLMRs_annotonly_Df$annot.end)))


AllSkinLMR_Unique_annos <- list()
for(i in levels(as.factor(AllSkinLMRs_annotonly_Df$annot.type))){
  anno <- AllSkinLMRs_annotonly_Df[AllSkinLMRs_annotonly_Df$annot.type == i,]
  anno <- anno[!duplicated(anno$coordname),]
  AllSkinLMR_Unique_annos[[i]] <- anno
  
}


AllSkinLMR_Unique_annos <- as.data.frame(do.call(rbind, AllSkinLMR_Unique_annos))

table(AllSkinLMR_Unique_annos$annot.type)



#removing duplicate exons 

#extracting from the exon coordinate names the ones that are not present in the first introns
uniqueExons <- AllSkinLMR_Unique_annos$coordname[AllSkinLMR_Unique_annos$annot.type == "exons"][which(!AllSkinLMR_Unique_annos$coordname[AllSkinLMR_Unique_annos$annot.type == "exons"] %in% AllSkinLMR_Unique_annos$coordname[AllSkinLMR_Unique_annos$annot.type == "firstexons"])]

#creating a separate exons table to edit and get the other exons
SkinLMR_Unique_AllExons <- AllSkinLMR_Unique_annos[AllSkinLMR_Unique_annos$annot.type == "exons",]

SkinLMR_Unique_AllExons <- SkinLMR_Unique_AllExons[SkinLMR_Unique_AllExons$coordname %in% uniqueExons,]

SkinLMR_Unique_AllExons$annot.type <- rep("other.exons", length(SkinLMR_Unique_AllExons$annot.type))


#removing exons from main table and adding in new one

AllSkinLMR_Unique_annos_ex <- AllSkinLMR_Unique_annos[!AllSkinLMR_Unique_annos$annot.type == "exons",]

AllSkinLMR_Unique_annos_ex <- rbind(AllSkinLMR_Unique_annos_ex, SkinLMR_Unique_AllExons)

table(AllSkinLMR_Unique_annos_ex$annot.type)



####writing out AllSkin features for GREAT background use ####

#(GREAT requires the input samples to also be in the background, so the background needs to be all the annos of the LMRs, not just all LMR coords)
AllSkin_CREs <- AllSkinLMR_Unique_annos_ex

AllSkin_CREs <- AllSkin_CREs[,1:3]

AllSkin_CREs$names <- AllSkinLMR_Unique_annos_ex$annot.tx_id

AllSkin_CREs <- AllSkin_CREs[AllSkinLMR_Unique_annos_ex$annot.type == "promoters" | AllSkinLMR_Unique_annos_ex$annot.type == "enhancers_fantom" | AllSkinLMR_Unique_annos_ex$annot.type == "firstexons",]




write.table(AllSkin_CREs, "AllSKin_LMR_CisRegs_Fexons_filt.bed",quote = F,row.names = F)






#### Finding proportions for plot later ####

#filtering just for bp and annotation type
AllSkinLMR_Unique_Annos_ex_fiddle <- AllSkinLMR_Unique_annos_ex[,c("annot.width","annot.type")]

AllSkinLMR_Unique_Annos_ex_fiddle <- AllSkinLMR_Unique_Annos_ex_fiddle %>%
  pivot_wider(names_from = annot.type, values_from = annot.width, values_fn = sum)


AllSkinLMR_Unique_Annos_ex_fiddle <- melt(AllSkinLMR_Unique_Annos_ex_fiddle)

#sort into categories

#CpGs only
AllSkin_LMR_CpGanno <- AllSkinLMR_Unique_Annos_ex_fiddle[3:6,]
#Genomic only
AllSkin_LMR_Genomic.anno <- AllSkinLMR_Unique_Annos_ex_fiddle[-(3:6),]


#dividing 
AllSkin_LMR_CpGanno$Proportion <- (AllSkin_LMR_CpGanno$value/sum(AllSkin_LMR_CpGanno$value))*100
AllSkin_LMR_CpGanno$Group <- rep("Skin LMR CpGs", length(AllSkin_LMR_CpGanno$variable))

AllSkin_LMR_Genomic.anno$Proportion <- (AllSkin_LMR_Genomic.anno$value/sum(AllSkin_LMR_Genomic.anno$value))*100
AllSkin_LMR_Genomic.anno$Group <- rep("Skin LMR Genomic", length(AllSkin_LMR_Genomic.anno$variable))





####Now all fibro LMRs for background ###############

#reading in all Fibro LMRs
AllFibroLMRs = read.csv('Fibros_fv_summary.csv')

#turning the bed file into GR ranges object
AllFibroLMRs_gr <- GRanges(seqnames = AllFibroLMRs$ch,
                           ranges = IRanges(start = AllFibroLMRs$b, end = AllFibroLMRs$e),
                           strand = rep("*", length(AllFibroLMRs$ch)),
                           number = seq(1:length(AllFibroLMRs$ch))) #corresponds to row number 





#now to annotate our LMRs.
AllFibroLMRs_annotated = annotate_regions(
  regions = AllFibroLMRs_gr,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)


AllFibroLMRs_annotated.df <- data.frame(AllFibroLMRs_annotated)


AllFibroLMRs_annotated.df_filt <- AllFibroLMRs_annotated.df[!AllFibroLMRs_annotated.df$annot.type %in% delsies,]

#cleaning up the terms
AllFibroLMRs_annotated.df_filt$annot.type <- gsub("hg38_", "", AllFibroLMRs_annotated.df_filt$annot.type)

AllFibroLMRs_annotated.df_filt$annot.type <- gsub("genes_", "", AllFibroLMRs_annotated.df_filt$annot.type)



AllFibroLMRs_annotonly_Df <- AllFibroLMRs_annotated.df_filt[,-(1:6)]

#making unified names (nothing fancy) to make the matching easier
AllFibroLMRs_annotonly_Df$coordname <- paste(AllFibroLMRs_annotonly_Df$annot.seqnames,(paste(AllFibroLMRs_annotonly_Df$annot.start, AllFibroLMRs_annotonly_Df$annot.end)))


AllFibroLMR_Unique_annos <- list()
for(i in levels(as.factor(AllFibroLMRs_annotonly_Df$annot.type))){
  anno <- AllFibroLMRs_annotonly_Df[AllFibroLMRs_annotonly_Df$annot.type == i,]
  anno <- anno[!duplicated(anno$coordname),]
  AllFibroLMR_Unique_annos[[i]] <- anno
  
}


AllFibroLMR_Unique_annos <- as.data.frame(do.call(rbind, AllFibroLMR_Unique_annos))

table(AllFibroLMR_Unique_annos$annot.type)



#removing duplicate exons 

#extracting from the exon coordinate names the ones that are not present in the first introns
uniqueExons <- AllFibroLMR_Unique_annos$coordname[AllFibroLMR_Unique_annos$annot.type == "exons"][which(!AllFibroLMR_Unique_annos$coordname[AllFibroLMR_Unique_annos$annot.type == "exons"] %in% AllFibroLMR_Unique_annos$coordname[AllFibroLMR_Unique_annos$annot.type == "firstexons"])]

#creating a separate exons table to edit and get the other exons
FibroLMR_Unique_AllExons <- AllFibroLMR_Unique_annos[AllFibroLMR_Unique_annos$annot.type == "exons",]

FibroLMR_Unique_AllExons <- FibroLMR_Unique_AllExons[FibroLMR_Unique_AllExons$coordname %in% uniqueExons,]

FibroLMR_Unique_AllExons$annot.type <- rep("other.exons", length(FibroLMR_Unique_AllExons$annot.type))


#removing exons from main table and adding in new one

AllFibroLMR_Unique_annos_ex <- AllFibroLMR_Unique_annos[!AllFibroLMR_Unique_annos$annot.type == "exons",]

AllFibroLMR_Unique_annos_ex <- rbind(AllFibroLMR_Unique_annos_ex, FibroLMR_Unique_AllExons)

table(AllFibroLMR_Unique_annos_ex$annot.type)



####writing out AllFibro features for GREAT background use ####

#(GREAT requires the input samples to also be in the background, so the background needs to be all the annos of the LMRs, not just all LMR coords)
AllFibro_CREs <- AllFibroLMR_Unique_annos_ex


AllFibro_CREs <- AllFibro_CREs[,1:3]

AllFibro_CREs$names <- AllFibroLMR_Unique_annos_ex$annot.tx_id

AllFibro_CREs <- AllFibro_CREs[AllFibroLMR_Unique_annos_ex$annot.type == "promoters" | AllFibroLMR_Unique_annos_ex$annot.type == "enhancers_fantom" | AllFibroLMR_Unique_annos_ex$annot.type == "firstexons",]




write.table(AllFibro_CREs, "AllFibro_LMR_CisRegs_Fexons_filt.bed",quote = F,row.names = F)




###### Finding proportions for plotting later ######

#filtering just for bp and annotation type
AllFibroLMR_Unique_Annos_ex_fiddle <- AllFibroLMR_Unique_annos_ex[,c("annot.width","annot.type")]

AllFibroLMR_Unique_Annos_ex_fiddle <- AllFibroLMR_Unique_Annos_ex_fiddle %>%
  pivot_wider(names_from = annot.type, values_from = annot.width, values_fn = sum)


AllFibroLMR_Unique_Annos_ex_fiddle <- melt(AllFibroLMR_Unique_Annos_ex_fiddle)

#sort into categories

#CpG only
AllFibro_LMR_CpGanno <- AllFibroLMR_Unique_Annos_ex_fiddle[3:6,]
#Genomic only
AllFibro_LMR_Genomic.anno <- AllFibroLMR_Unique_Annos_ex_fiddle[-(3:6),]


#dividing 
AllFibro_LMR_CpGanno$Proportion <- (AllFibro_LMR_CpGanno$value/sum(AllFibro_LMR_CpGanno$value))*100
AllFibro_LMR_CpGanno$Group <- rep("Fibro LMR CpGs", length(AllFibro_LMR_CpGanno$variable))

AllFibro_LMR_Genomic.anno$Proportion <- (AllFibro_LMR_Genomic.anno$value/sum(AllFibro_LMR_Genomic.anno$value))*100
AllFibro_LMR_Genomic.anno$Group <- rep("Fibro LMR Genomic", length(AllFibro_LMR_Genomic.anno$variable))





############ Assembling all background LMRs for plotting ############


#doing CpGs first

#bld
AllLMR_Proport_CpG <- AllBld_LMR_CpGanno[,-4]


colnames(AllLMR_Proport_CpG)[2:3] <- c("Bld_AllLMR_bp", "Bld_AllLMR_Proportion")

#skin

AllLMR_Proport_CpG <- cbind(AllLMR_Proport_CpG,AllSkin_LMR_CpGanno[,2:3])

colnames(AllLMR_Proport_CpG)[4:5] <- c("Skin_AllLMR_bp", "Skin_AllLMR_Proportion")


#Fibro

AllLMR_Proport_CpG <- cbind(AllLMR_Proport_CpG,AllFibro_LMR_CpGanno[,2:3])

colnames(AllLMR_Proport_CpG)[6:7] <- c("Fibro_AllLMR_bp", "Fibro_AllLMR_Proportion")




#Adding in corresponding top1000 LMRs proportions
AllLMR_Proport_CpG$BloodLMR_prop <- Bld_LMR_CpGanno$Proportion[match(AllLMR_Proport_CpG$variable,Bld_LMR_CpGanno$variable)]

AllLMR_Proport_CpG$EpidLMR_prop <- Skin_LMR_CpGanno$Proportion[match(AllLMR_Proport_CpG$variable,Skin_LMR_CpGanno$variable)]

AllLMR_Proport_CpG$FibroLMR_prop <- Fibro_LMR_CpGanno$Proportion[match(AllLMR_Proport_CpG$variable,Fibro_LMR_CpGanno$variable)]


##Geno now


AllLMR_Proport_Geno<- AllBld_LMR_Genomic.anno[,-4]


colnames(AllLMR_Proport_Geno)[2:3] <- c("Bld_AllLMR_bp", "Bld_AllLMR_Proportion")

#skin

AllLMR_Proport_Geno <- cbind(AllLMR_Proport_Geno,AllSkin_LMR_Genomic.anno[,2:3])

colnames(AllLMR_Proport_Geno)[4:5] <- c("Skin_AllLMR_bp", "Skin_AllLMR_Proportion")


#Fibro

AllLMR_Proport_Geno <- cbind(AllLMR_Proport_Geno,AllFibro_LMR_Genomic.anno[,2:3])

colnames(AllLMR_Proport_Geno)[6:7] <- c("Fibro_AllLMR_bp", "Fibro_AllLMR_Proportion")


#Adding in corresponding top1000 LMRs
AllLMR_Proport_Geno$BloodLMR_prop <- Bld_LMR_Genomic.anno$Proportion[match(AllLMR_Proport_Geno$variable,Bld_LMR_Genomic.anno$variable)]

AllLMR_Proport_Geno$EpidLMR_prop <- Skin_LMR_Genomic.anno$Proportion[match(AllLMR_Proport_Geno$variable,Skin_LMR_Genomic.anno$variable)]

AllLMR_Proport_Geno$FibroLMR_prop <- Fibro_LMR_Genomic.anno$Proportion[match(AllLMR_Proport_Geno$variable,Fibro_LMR_Genomic.anno$variable)]






#Merging both here
AllLMR_Proport_All <- rbind(AllLMR_Proport_Geno,AllLMR_Proport_CpG)



#Calculating fold change as eitehr hyper or hypo divided by the corresponding chromatin value from the whole genome.
#if hyper/hypo is bigger than genome, then fold change is greater than 1
#doing log2 of this results in a clear positive vs negative, so when hyper/hypo is less than 1 ie less than the genome, it's negative


#Blood
AllLMR_Proport_All$Blood_Log_Fold_Change <- log2(AllLMR_Proport_All$BloodLMR_prop/AllLMR_Proport_All$Bld_AllLMR_Proportion)

#For the purposes of plotting, all infinity values are changed to 0
AllLMR_Proport_All$Blood_Log_Fold_Change[AllLMR_Proport_All$Blood_Log_Fold_Change == "-Inf"] <- 0


#Epiderm
AllLMR_Proport_All$Epid_Log_Fold_Change <- log2(AllLMR_Proport_All$EpidLMR_prop/AllLMR_Proport_All$Skin_AllLMR_Proportion)


#Fibro
AllLMR_Proport_All$Fibro_Log_Fold_Change <- log2(AllLMR_Proport_All$FibroLMR_prop/AllLMR_Proport_All$Fibro_AllLMR_Proportion)




#getting just the values we want to plot

AllLMR_Proport_All_plt <- AllLMR_Proport_All[,c("variable","Blood_Log_Fold_Change","Epid_Log_Fold_Change" ,"Fibro_Log_Fold_Change" )]

colnames(AllLMR_Proport_All_plt) <- c("Feature","Blood LMRs","Epidermis LMRs", "Fibroblasts LMRs" )

library(reshape2)
AllLMR_Proport_All_plt <- melt(AllLMR_Proport_All_plt)

#redoing the levels cos CpG stuff still there

#reordering genomic features by genic, regulatory/non coding

ordered <- c("firstexons","other.exons","introns", "promoters","enhancers_fantom",
             "lncrna_gencode","intergenic","5UTRs","3UTRs","cpg_inter", "cpg_islands","cpg_shelves","cpg_shores")


#Reordering after talking to Andrea


ordered <- c("cpg_islands","cpg_shelves","cpg_shores", "cpg_inter", #CpG
             "5UTRs","3UTRs","firstexons","other.exons","introns","lncrna_gencode", #Genic
             "promoters","enhancers_fantom", "intergenic") #regulatory/non-coding




AllLMR_Proport_All_plt$Feature <- factor(AllLMR_Proport_All_plt$Feature, levels = ordered)


############ Plotting Annotations all ############

library(ggplot2)

#Plotting just Blood and Epiderm
pdf(file = "Blood_Epid_LMR_logGeom_CpG_LMR-BACKG.pdf",width = 10,height = 6, useDingbats=FALSE)
ggplot(AllLMR_Proport_All_plt[1:26,], aes(fill=variable, y=value, x=Feature)) + 
  geom_bar(position="dodge", stat="identity") + theme_classic() +theme(axis.text = element_text(size = 14),axis.title = element_text(size = 14), legend.text = element_text(size = 10),axis.text.x=element_text(colour="black")) +
  scale_fill_manual(values = c("#C03F3F","#3671A8")) +theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ylab("Log2 Percentage Fold Change")

dev.off()



#Plotting blood, epidermis and fibro
pdf(file = "Blood_Epid_Fibro_LMR_logGeom_CpG_LMR-BACKGver3.pdf",width = 10,height = 6, useDingbats=FALSE)
ggplot(AllLMR_Proport_All_plt, aes(fill=variable, y=value, x=Feature)) + 
  geom_bar(position="dodge", stat="identity") + theme_classic() +theme(axis.text = element_text(size = 14),axis.title = element_text(size = 14), legend.text = element_text(size = 10),axis.text.x=element_text(colour="black")) +
  scale_fill_manual(values = c("#C03F3F","#3671A8","#E2812C")) +theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ylab("Log2 Percentage Fold Change")

dev.off()







################## Looking at overlaps between LMR datasets ##################



############## Annotating skin and blood with individual CpGs  ##############

#https://support.bioconductor.org/p/95239/

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

library(BSgenome.Hsapiens.UCSC.hg38)  
chrs <- names(Hsapiens)[1:24]
cgs <- lapply(chrs, function(x) start(matchPattern("CG", Hsapiens[[x]])))

cpgr <- do.call(c, lapply(1:24, function(x) GRanges(names(Hsapiens)[x], IRanges(cgs[[x]], width = 2))))

#now to annotate our LMRs. _ Blood
BloodLMRs_CpGanno = annotate_regions(
  regions = BloodLMRs1000_gr,
  annotations = cpgr,
  ignore.strand = TRUE,
  quiet = FALSE)


#Skin CpG anno
SkinLMRs_CpGanno = annotate_regions(
  regions = SkinLMRs1000_gr,
  annotations = cpgr,
  ignore.strand = TRUE,
  quiet = FALSE)



#Now to get the overlaps between CpGs, I think the best way is to extract the annotations only, then calculate the overlaps


BloodLMRs_CpGanno_df <- as.data.frame(BloodLMRs_CpGanno)

BloodLMRs_CpGanno_df <- BloodLMRs_CpGanno_df[,-(1:8)]

BloodLMRs_CpGanno_df <- GRanges(seqnames = BloodLMRs_CpGanno_df$annot.seqnames,
                                ranges = IRanges(start = BloodLMRs_CpGanno_df$annot.start, 
                                                 end = BloodLMRs_CpGanno_df$annot.end),
                                strand = BloodLMRs_CpGanno_df$annot.strand) #corresponds to row number 





SkinLMRs_CpGanno_df <- as.data.frame(SkinLMRs_CpGanno)

SkinLMRs_CpGanno_df <- SkinLMRs_CpGanno_df[,-(1:8)]

SkinLMRs_CpGanno_df <- GRanges(seqnames = SkinLMRs_CpGanno_df$annot.seqnames,
                               ranges = IRanges(start = SkinLMRs_CpGanno_df$annot.start, end = SkinLMRs_CpGanno_df$annot.end),
                               strand = SkinLMRs_CpGanno_df$annot.strand) #corresponds to row number 



# Creating venn overlap
BiocManager::install("ChIPpeakAnno")
install.packages("devtools"); library(devtools)

#YOu also have to download Rtools separately and install for your version of R (4.2 in this case)

install_github("js229/Vennerable"); library(Vennerable)

library(ChIPpeakAnno)

BvS_venn <- makeVennDiagram(Peaks=list(BloodLMRs_CpGanno_df, SkinLMRs_CpGanno_df),
                            NameOfPeaks=c("T cells LMRs", "Epidermis LMRs"))


venn_cnt2venn <- function(venn_cnt){
  n <- which(colnames(venn_cnt)=="Counts") - 1
  SetNames=colnames(venn_cnt)[1:n]
  Weight=venn_cnt[,"Counts"]
  names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
  compute.Venn(Venn(SetNames=SetNames, Weight=Weight))
}

v <- venn_cnt2venn(BvS_venn$vennCounts)


#setting up colours manually https://stackoverflow.com/questions/45341383/manually-set-color-for-venn-diagrams-in-vennerable-package
gp <- VennThemes(v,"sequential")

#old colours
#gp[["Face"]][["11"]]$fill <-  "#B5D66C" #green overlap
#gp[["Face"]][["01"]]$fill <-  "#3095E5" #blue (epid)
#gp[["Face"]][["10"]]$fill <-  "#E0EB43" #yellow tcells


# "#C03F3F","#3671A8","#E2812C" all the colours we're using so far

gp[["Face"]][["11"]]$fill <-  "#8C88A3" #purple overlap
gp[["Face"]][["01"]]$fill <-  "#6AA4D3" #blue (epid)
gp[["Face"]][["10"]]$fill <-  "#C03F3F" #red tcells



#setting text colours
gp$SetText$Set1$col <- "black"
gp$SetText$Set2$col <- "black"

#outlines removed
gp$Set$Set1$lty <- 0
gp$Set$Set2$lty <- 0


pdf(file = "Bld_SkinLMRs_Venn2.pdf",width = 6,height = 6, useDingbats=FALSE)

plot(v, gp= gp)

dev.off()









##### annotating Fibro with CpGs to check for overlaps ##############


#Fibro CpG anno
FibroLMRs_CpGanno = annotate_regions(
  regions = FibroLMRs1000_gr,
  annotations = cpgr,
  ignore.strand = TRUE,
  quiet = FALSE)




FibroLMRs_CpGanno_df <- as.data.frame(FibroLMRs_CpGanno)

FibroLMRs_CpGanno_df <- FibroLMRs_CpGanno_df[,-(1:8)]

FibroLMRs_CpGanno_df <- GRanges(seqnames = FibroLMRs_CpGanno_df$annot.seqnames,
                                ranges = IRanges(start = FibroLMRs_CpGanno_df$annot.start, 
                                                 end = FibroLMRs_CpGanno_df$annot.end),
                                strand = FibroLMRs_CpGanno_df$annot.strand) #corresponds to row number 





###Proportional plot then manually colour later



BvSvF_venn <- makeVennDiagram(Peaks=list(BloodLMRs_CpGanno_df, SkinLMRs_CpGanno_df, FibroLMRs_CpGanno_df),
                              NameOfPeaks=c("T cells LMRs", "Epidermis LMRs", "Passaged Fibroblast LMRs"),
                              scaled = T,
                              fill=c("#E0EB43", "#3095E5","#fa7c1b"), # circle fill color
                              #col=c("#000000FF", "#000000FF","#000000FF")#circle border color
                              lty = 'blank' #no outline
                              # sub.fontfamily = "ariel", main.fontfamily = "ariel" # this didnt work
)

#BvS_venn <- makeVennDiagram(Peaks=list(BloodLMRs_CpGanno_df, SkinLMRs_CpGanno_df),
#                          NameOfPeaks=c("T cells LMRs", "Epidermis LMRs"))




v2 <- venn_cnt2venn(BvSvF_venn$vennCounts)


#setting up colours manually https://stackoverflow.com/questions/45341383/manually-set-color-for-venn-diagrams-in-vennerable-package
gp2 <- VennThemes(v2,"sequential")

#old colours
#gp[["Face"]][["11"]]$fill <-  "#B5D66C" #green overlap
#gp[["Face"]][["01"]]$fill <-  "#3095E5" #blue (epid)
#gp[["Face"]][["10"]]$fill <-  "#E0EB43" #yellow tcells


# "#C03F3F","#3671A8","#E2812C" all the colours we're using so far

gp2[["Face"]][["11"]]$fill <-  "#993399" #green overlap
gp2[["Face"]][["01"]]$fill <-  "#3671A8" #blue (epid)
gp2[["Face"]][["10"]]$fill <-  "#C03F3F" #red tcells



#setting text colours
gp2$SetText$Set1$col <- "black"
gp2$SetText$Set2$col <- "black"

#outlines removed
gp2$Set$Set1$lty <- 0
gp2$Set$Set2$lty <- 0
gp2$Set$Set3$lty <- 0


pdf(file = "Bld_Skin_fibroLMRs_Venn.pdf",width = 6,height = 6, useDingbats=FALSE)

plot(v2, gp= gp2)

dev.off()

#manually recoloured in illustrator





############## Annotating Fibro_chip with CpGs to check for overlaps ##############

#reading in Fibro chip (fv) LMRs
Fibro_chip_all = read.csv('fv_FibChIP_summary.csv',header = T)




Fibro_chip_all <- GRanges(seqnames = Fibro_chip_all$ch,
                          ranges = IRanges(start = Fibro_chip_all$b, end = Fibro_chip_all$e),
                          strand = rep("*", length(Fibro_chip_all$ch)),
                          pr = Fibro_chip_all$pr,
                          delta = Fibro_chip_all$d,
                          number = seq(1:length(Fibro_chip_all$ch))) #corresponds to row number 


#taking top 1000, already sorted for PRC2 binding
Fibro_chip_1000 <- Fibro_chip_all[1:1000,]


#Fibro_chip CpG anno
Fibro_chipLMRs_CpGanno = annotate_regions(
  regions = Fibro_chip_1000,
  annotations = cpgr,
  ignore.strand = TRUE,
  quiet = FALSE)




Fibro_chipLMRs_CpGanno_df <- as.data.frame(Fibro_chipLMRs_CpGanno)

Fibro_chipLMRs_CpGanno_df <- Fibro_chipLMRs_CpGanno_df[,-(1:8)]

Fibro_chipLMRs_CpGanno_df <- GRanges(seqnames = Fibro_chipLMRs_CpGanno_df$annot.seqnames,
                                     ranges = IRanges(start = Fibro_chipLMRs_CpGanno_df$annot.start, 
                                                      end = Fibro_chipLMRs_CpGanno_df$annot.end),
                                     strand = Fibro_chipLMRs_CpGanno_df$annot.strand) #corresponds to row number 




##Lets take on this VENNture....ehehehehe


FvFch_venn <- makeVennDiagram(Peaks=list(FibroLMRs_CpGanno_df, Fibro_chipLMRs_CpGanno_df),
                              NameOfPeaks=c("Fibroblast Emb LMRs", "Fibroblast ChIP LMRs"))


v3 <- venn_cnt2venn(FvFch_venn$vennCounts)


#setting up colours manually https://stackoverflow.com/questions/45341383/manually-set-color-for-venn-diagrams-in-vennerable-package
gp3 <- VennThemes(v3,"sequential")

#old colours
#gp[["Face"]][["11"]]$fill <-  "#B5D66C" #green overlap
#gp[["Face"]][["01"]]$fill <-  "#3095E5" #blue (epid)
#gp[["Face"]][["10"]]$fill <-  "#E0EB43" #yellow tcells


# "#C03F3F","#3671A8","#E2812C" all the colours we're using so far

gp3[["Face"]][["11"]]$fill <-  "#BB5225" 
gp3[["Face"]][["01"]]$fill <-  "#BD6042" 
gp3[["Face"]][["10"]]$fill <-  "#F9BD8D" 



#setting text colours
gp3$SetText$Set1$col <- "black"
gp3$SetText$Set2$col <- "black"

#outlines removed
gp3$Set$Set1$lty <- 0
gp3$Set$Set2$lty <- 0


pdf(file = "Fibro_vsFchip_LMRs_Venn.pdf",width = 6,height = 6, useDingbats=FALSE)

plot(v3, gp= gp3)

dev.off()

#calc how much is shared

#fibro emb
114684/(145972+114684)
#44%

#fibrochip
114684/(120374+114684)
#49%











