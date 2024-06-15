
###unofficial recreate ver #####

#### Diffbind neo vs old ######

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DiffBind")
BiocManager::install("profileplyr")


library(DiffBind)
library(profileplyr)

NEOvsOLD_P2_test <- read.csv("../All_Fibro_Pass.csv")


#mapped, deduplicated ChIP files must be saved in this working directory for All_Fibro_Pass.csv to have the correct file paths for the bamReads and bamControl columns

NEOvsOLD_P2_test <- dba(sampleSheet = NEOvsOLD_P2_test)

#applying black and grey lists
NEOvsOLD_P2_test <- dba.blacklist(NEOvsOLD_P2_test,blacklist = T,greylist = T)

NEOvsOLD_P2_testCnt<-dba.count(NEOvsOLD_P2_test)

NEOvsOLD_P2_testCnt <- dba.normalize(NEOvsOLD_P2_testCnt)

NEOvsOLD_P2_testCnt <- dba.contrast(NEOvsOLD_P2_testCnt, minMembers = 2, categories = DBA_TISSUE)

NEOvsOLD_P2_testCnt <- dba.analyze(NEOvsOLD_P2_testCnt, method = DBA_ALL_METHODS)


dba.show(NEOvsOLD_P2_testCnt, bContrasts=T)



###### NEO & OLD LMRs shared - hESC ordered, Plotting top 1000 LMRs (High-PRC2 LMRs) for both NEO and OLD #######


NEO_OLD_EmbOrdLMRs <- read.csv("../ALL_fibropass_methBWs_ENC_hESC_EZH2suz12_FE.csv")

NEO_OLD_EmbOrdLMRs <- NEO_OLD_EmbOrdLMRs[order(NEO_OLD_EmbOrdLMRs$pr, decreasing = T),]

#remove  empty rows added at end of dataframe
NEO_OLD_EmbOrdLMRs <- na.omit(NEO_OLD_EmbOrdLMRs)

NEO_OLD_EmbOrdLMRs <- GRanges(seqnames = NEO_OLD_EmbOrdLMRs$ch,
                              ranges = IRanges(start = NEO_OLD_EmbOrdLMRs$b, end = NEO_OLD_EmbOrdLMRs$e),
                              strand = rep("*", length(NEO_OLD_EmbOrdLMRs$ch)),
                              pr = NEO_OLD_EmbOrdLMRs$pr,
                              delta = NEO_OLD_EmbOrdLMRs$d,
                              number = seq(1:length(NEO_OLD_EmbOrdLMRs$ch)),
                              index = NEO_OLD_EmbOrdLMRs$index,
                              NEO2P2=NEO_OLD_EmbOrdLMRs$NEO2_P2_meth,
                              NEO2P5=NEO_OLD_EmbOrdLMRs$NEO2_P5_meth,
                              NEO2P8=NEO_OLD_EmbOrdLMRs$NEO2_P8_meth,
                              OLD3P2=NEO_OLD_EmbOrdLMRs$OLD3_P2_meth,
                              OLD3P5=NEO_OLD_EmbOrdLMRs$OLD3_P5_meth,
                              OLD3P8=NEO_OLD_EmbOrdLMRs$OLD3_P8_meth) 



#plotting Emb sorted top 1000 

Neo_oldEMB_Top1000_plot <- dba.plotProfile(NEOvsOLD_P2_testCnt,sites = NEO_OLD_EmbOrdLMRs[1:1000,]) 

png("NEO_OLDLMRs_Embrank_NvOChIPProfile.png", width = 4, height =5, units = "in", res = 400)
dba.plotProfile(Neo_oldEMB_Top1000_plot)
dev.off()

pdf("NEO_OLDLMRs_Embrank_NvOChIPProfile.pdf", width = 4, height =5)
dba.plotProfile(Neo_oldEMB_Top1000_plot)
dev.off()


#### Figure 6 heatmaps #####


#### NEO overlap heatmap 

#### NEO P2 emb ordered LMRs

NEO_EmbOrdLMRs <- read.csv("../results/NEOP2s_methBWs_hESCord_FE.csv")

NEO_EmbOrdLMRs <- NEO_EmbOrdLMRs[order(NEO_EmbOrdLMRs$pr, decreasing = T),]


NEO_EmbOrdLMRs <- GRanges(seqnames = NEO_EmbOrdLMRs$ch,
                          ranges = IRanges(start = NEO_EmbOrdLMRs$b, end = NEO_EmbOrdLMRs$e),
                          strand = rep("*", length(NEO_EmbOrdLMRs$ch)),
                          pr = NEO_EmbOrdLMRs$pr,
                          delta = NEO_EmbOrdLMRs$d,
                          number = seq(1:length(NEO_EmbOrdLMRs$ch)),
                          index = NEO_EmbOrdLMRs$index,
                          NEO2P2=NEO_EmbOrdLMRs$NEO2_P2_meth,
                          NEO2P5=NEO_EmbOrdLMRs$NEO2_P5_meth,
                          NEO2P8=NEO_EmbOrdLMRs$NEO2_P8_meth) 



#### NEOP2 ordered NEO LMRs #############

NEO_NEOP2OrdLMRs <- read.csv("../results/NEOP2s_methBWs_NEOMergeP2_FE.csv")

NEO_NEOP2OrdLMRs <- NEO_NEOP2OrdLMRs[order(NEO_NEOP2OrdLMRs$pr, decreasing = T),]

NEO_NEOP2OrdLMRs <- GRanges(seqnames = NEO_NEOP2OrdLMRs$ch,
                            ranges = IRanges(start = NEO_NEOP2OrdLMRs$b, end = NEO_NEOP2OrdLMRs$e),
                            strand = rep("*", length(NEO_NEOP2OrdLMRs$ch)),
                            pr = NEO_NEOP2OrdLMRs$pr,
                            delta = NEO_NEOP2OrdLMRs$d,
                            number = seq(1:length(NEO_NEOP2OrdLMRs$ch)),
                            index = NEO_NEOP2OrdLMRs$index,
                            NEO2P2=NEO_NEOP2OrdLMRs$NEO2_P2_meth,
                            NEO2P5=NEO_NEOP2OrdLMRs$NEO2_P5_meth,
                            NEO2P8=NEO_NEOP2OrdLMRs$NEO2_P8_meth)




NEOEMB_Top1000_ove <- as.data.frame(findOverlaps(NEO_NEOP2OrdLMRs,NEO_EmbOrdLMRs))

colnames(NEOEMB_Top1000_ove) <- c("NEO_Rank", "NEO_Rank_Emb_Order")


library(pheatmap)
library(RColorBrewer)
#setting breaks for easier colour stuff
breaksList = seq(0, length(rownames(NEOEMB_Top1000_ove)), by = 1)

#plotting full LMR list of LMRs ordered (and coloured by) Neonatal passage 2 EZH2 binding (left column), then plotting same color but ordered by hESC binding (right column)

#png("NEOrank_Embproject.png", width = 4, height =8, units = "in", res = 400)
pdf("NEOrank_Embproject.pdf", width = 4, height =8)
pheatmap(NEOEMB_Top1000_ove,cluster_cols = F,cluster_rows = F,labels_col = "",#c("T-cell LMRs ranked by CD4 PRC2 Binding", "T-Cell CD4 LMRs in hESC Ranking"), 
         labels_row = seq(0,length(rownames(NEOEMB_Top1000_ove)), by = 1000), #labels_row = "",
         color = colorRampPalette(brewer.pal(n = 7, name = "PuOr"))(length(breaksList)),
         border_color = NA)
dev.off()

#plotting top 1000 LMRs only as "zoomed in" version

#png("NEOrank_Embproject_top1000.png", width = 4, height =1.2, units = "in", res = 400)
pdf("NEOrank_Embproject_top1000.pdf", width = 4, height =1.2)
pheatmap(NEOEMB_Top1000_ove[1:1000,],cluster_cols = F,cluster_rows = F,labels_col = "", labels_row = seq(0,1000, by = 500),
         color = colorRampPalette(brewer.pal(n = 7, name = "PuOr"))(length(breaksList)),
         border_color = NA, breaks = breaksList)
dev.off()




####shared heatmap OLD ############

#old p1 LMRs, hESC EZH2 binding ordered
OLD_EmbOrdLMRs <- read.csv("../results/OLDP2s_methBWs_hESCord_FE.csv")

OLD_EmbOrdLMRs <- OLD_EmbOrdLMRs[order(OLD_EmbOrdLMRs$pr, decreasing = T),]


OLD_EmbOrdLMRs <- GRanges(seqnames = OLD_EmbOrdLMRs$ch,
                          ranges = IRanges(start = OLD_EmbOrdLMRs$b, end = OLD_EmbOrdLMRs$e),
                          strand = rep("*", length(OLD_EmbOrdLMRs$ch)),
                          pr = OLD_EmbOrdLMRs$pr,
                          delta = OLD_EmbOrdLMRs$d,
                          number = seq(1:length(OLD_EmbOrdLMRs$ch)),
                          index = OLD_EmbOrdLMRs$index,
                          OLD3P2=OLD_EmbOrdLMRs$OLD3_P2_meth,
                          OLD3P5=OLD_EmbOrdLMRs$OLD3_P5_meth,
                          OLD3P8=OLD_EmbOrdLMRs$OLD3_P8_meth) 

#Old P2 lmrs ordered by old p1 EZH2 binding

OLD_OLDP2OrdLMRs <- read.csv("../results/OLDP2s_methBWs_OLDMergeP2_FE.csv")

OLD_OLDP2OrdLMRs <- OLD_OLDP2OrdLMRs[order(OLD_OLDP2OrdLMRs$pr, decreasing = T),]

OLD_OLDP2OrdLMRs <- GRanges(seqnames = OLD_OLDP2OrdLMRs$ch,
                            ranges = IRanges(start = OLD_OLDP2OrdLMRs$b, end = OLD_OLDP2OrdLMRs$e),
                            strand = rep("*", length(OLD_OLDP2OrdLMRs$ch)),
                            pr = OLD_OLDP2OrdLMRs$pr,
                            delta = OLD_OLDP2OrdLMRs$d,
                            number = seq(1:length(OLD_OLDP2OrdLMRs$ch)),
                            index = OLD_OLDP2OrdLMRs$index,
                            OLD3P2=OLD_OLDP2OrdLMRs$OLD3_P2_meth,
                            OLD3P5=OLD_OLDP2OrdLMRs$OLD3_P5_meth,
                            OLD3P8=OLD_OLDP2OrdLMRs$OLD3_P8_meth)



OLDEMB_Top1000_ove <- as.data.frame(findOverlaps(OLD_OLDP2OrdLMRs,OLD_EmbOrdLMRs))

colnames(OLDEMB_Top1000_ove) <- c("OLD_Rank", "OLD_Rank_Emb_Order")


library(pheatmap)
library(RColorBrewer)
#setting breaks for easier colour stuff
breaksList = seq(0, length(rownames(OLDEMB_Top1000_ove)), by = 1)

#plotting full LMR list of LMRs ordered (and coloured by) Old passage 2 EZH2 binding (left column), then plotting same color but ordered by hESC binding (right column)

#png("OLDrank_Embproject.png", width = 4, height =8, units = "in", res = 400)
pdf("OLDrank_Embproject.png", width = 4, height =8)
pheatmap(OLDEMB_Top1000_ove,cluster_cols = F,cluster_rows = F,labels_col = "",#c("T-cell LMRs ranked by CD4 PRC2 Binding", "T-Cell CD4 LMRs in hESC Ranking"), 
         labels_row = seq(0,length(rownames(OLDEMB_Top1000_ove)), by = 1000), #labels_row = "",
         color = colorRampPalette(brewer.pal(n = 7, name = "PuOr"))(length(breaksList)),
         border_color = NA)
dev.off()

#plotting top 1000 LMRs only as "zoomed in" version

#png("OLDrank_Embproject_top1000.png", width = 4, height =1.2, units = "in", res = 400)
pdf("OLDrank_Embproject_top1000.pdf", width = 4, height =1.2)
pheatmap(OLDEMB_Top1000_ove[1:1000,],cluster_cols = F,cluster_rows = F,labels_col = "", labels_row = seq(0,1000, by = 500),
         color = colorRampPalette(brewer.pal(n = 7, name = "PuOr"))(length(breaksList)),
         border_color = NA, breaks = breaksList)
dev.off()



####shared heatmap CD4 ############

#CD4 p1 LMRs, hESC EZH2 binding ordered
CD4_EmbOrdLMRs <- read.csv("../results/t.csv")

CD4_EmbOrdLMRs <- CD4_EmbOrdLMRs[order(CD4_EmbOrdLMRs$pr, decreasing = T),]


CD4_EmbOrdLMRs <- GRanges(seqnames = CD4_EmbOrdLMRs$ch,
                          ranges = IRanges(start = CD4_EmbOrdLMRs$b, end = CD4_EmbOrdLMRs$e),
                          strand = rep("*", length(CD4_EmbOrdLMRs$ch)),
                          pr = CD4_EmbOrdLMRs$pr,
                          delta = CD4_EmbOrdLMRs$d,
                          number = seq(1:length(CD4_EmbOrdLMRs$ch)),
                          index = CD4_EmbOrdLMRs$index,
                          X0=CD4_EmbOrdLMRs$X0,
                          X18=CD4_EmbOrdLMRs$X18,
                          X25=CD4_EmbOrdLMRs$X25,
                          X82=CD4_EmbOrdLMRs$X82,
                          X86=CD4_EmbOrdLMRs$X86,
                          X103=CD4_EmbOrdLMRs$X103) 

#CD4 P2 lmrs ordered by CD4 p1 EZH2 binding

CD4_CD4P2OrdLMRs <- read.csv("../results/t_cell_lmrs_CD4pooled_prc2_ENCinp_FE_meth.csv")

CD4_CD4P2OrdLMRs <- CD4_CD4P2OrdLMRs[order(CD4_CD4P2OrdLMRs$pr, decreasing = T),]

CD4_CD4P2OrdLMRs <- GRanges(seqnames = CD4_CD4P2OrdLMRs$ch,
                            ranges = IRanges(start = CD4_CD4P2OrdLMRs$b, end = CD4_CD4P2OrdLMRs$e),
                            strand = rep("*", length(CD4_CD4P2OrdLMRs$ch)),
                            pr = CD4_CD4P2OrdLMRs$pr,
                            delta = CD4_CD4P2OrdLMRs$d,
                            number = seq(1:length(CD4_CD4P2OrdLMRs$ch)),
                            index = CD4_CD4P2OrdLMRs$index,
                            X0=CD4_CD4P2OrdLMRs$X0,
                            X18=CD4_CD4P2OrdLMRs$X18,
                            X25=CD4_CD4P2OrdLMRs$X25,
                            X82=CD4_CD4P2OrdLMRs$X82,
                            X86=CD4_CD4P2OrdLMRs$X86,
                            X103=CD4_CD4P2OrdLMRs$X103)



CD4EMB_Top1000_ove <- as.data.frame(findOverlaps(CD4_CD4P2OrdLMRs,CD4_EmbOrdLMRs))

colnames(CD4EMB_Top1000_ove) <- c("CD4_Rank", "CD4_Rank_Emb_Order")


library(pheatmap)
library(RColorBrewer)
#setting breaks for easier colour stuff
breaksList = seq(0, length(rownames(CD4EMB_Top1000_ove)), by = 1)

#plotting full LMR list of LMRs ordered (and coloured by) CD4 passage 2 EZH2 binding (left column), then plotting same color but ordered by hESC binding (right column)

#png("CD4rank_Embproject.png", width = 4, height =8, units = "in", res = 400)
pdf("CD4rank_Embproject.pdf", width = 4, height =8)
pheatmap(CD4EMB_Top1000_ove,cluster_cols = F,cluster_rows = F,labels_col = "",#c("T-cell LMRs ranked by CD4 PRC2 Binding", "T-Cell CD4 LMRs in hESC Ranking"), 
         labels_row = seq(0,length(rownames(CD4EMB_Top1000_ove)), by = 1000), #labels_row = "",
         color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(length(breaksList)),
         border_color = NA)
dev.off()

#plotting top 1000 LMRs only as "zoomed in" version

#png("CD4rank_Embproject_top1000.png", width = 4, height =1.2, units = "in", res = 400)
pdf("CD4rank_Embproject_top1000.pdf", width = 4, height =1.2)
pheatmap(CD4EMB_Top1000_ove[1:1000,],cluster_cols = F,cluster_rows = F,labels_col = "", labels_row = seq(0,1000, by = 500),
         color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(length(breaksList)),
         border_color = NA, breaks = breaksList)
dev.off()

