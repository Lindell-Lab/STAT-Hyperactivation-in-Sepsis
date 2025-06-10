library(Seurat)
library(CellChat)
library(patchwork)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(clusterProfiler)
library(ReactomePA)
library(scplotter)

setwd("/lindell_lab/Data/Sepsis_scRNAseq_TCR/")

###load data
sub=readRDS("Seurat/Tcellsubset_soupx_annotated.RDS")

Idents(sub)="final_mods"

#Subset out HC
healthy=subset(sub,idents="Healthy")
healthy@meta.data$samples=healthy@meta.data$final_mods

#Subset out Sepsis samples
sepsis = subset(sub,idents="Sepsis")
sepsis@meta.data$samples= sepsis@meta.data$final_mods

CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB)


#Run cellchat
runcellchat=function(sub){
  cc <- createCellChat(object = sub, group.by = "sctype_classification", assay = "RNA")
  cc <- setIdent(cc, ident.use = "sctype_classification")
  cc@DB<-CellChatDB.use
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
}

cellchat.hc=runcellchat(healthy)
cellchat.sepsis=runcellchat(sepsis)

#merge cellchat object
object.list <- list(Healthy = cellchat.hc, Sepsis= cellchat.sepsis)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#save cellchat object
save(cellchat.hc, file = "cellchat_HC.RData")
save(cellchat.sepsis, file = "cellchat_Sepsis.RData")
save(cellchat, file = "cellchat_merged_HC_Sepsis.RData")

png("cellchat/interaction.png", width = 7,height=5,units='in',res = 300)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()

#red colored edges represent increased signaling in the second dataset compared to the first one. blue shows lower
#1- healthy #2 sepsis
png("cellchat/net_interaction.png", width = 10,height=10,units='in',res = 300)
p1=netVisual_diffInteraction(cellchat, weight.scale = T)
dev.off()

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 15)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",height=15)
png("plots/signalling.png", width = 10,height=10,units='in',res = 300)
ht1 + ht2
dev.off()


#### LR pairs ####
gg1 <- netVisual_bubble(cellchat, targets.use = c(1:3,5,9,10), sources.use = c(4,12),  comparison = c(1,2), max.dataset = 2, title.name = "Increased signaling in Sepsis", angle.x = 90, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, targets.use = c(1:3,5,9,10), sources.use = c(4,12),  comparison = c(1,2), max.dataset = 1, title.name = "Decreased signaling in Sepsis", angle.x = 90, remove.isolate = T)
png("cellchat_full/commprob_dotplot_allT_target.png", width = 10,height=8,units = "in",res=300)
gg1 + gg2
dev.off()

