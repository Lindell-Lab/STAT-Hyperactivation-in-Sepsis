library(Seurat)
library(scExtras)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(SCP)
library(ComplexHeatmap)
library(purrr)
library(future)
library(SoupX)
setwd("/lindell_lab/Data/Sepsis_scRNAseq_TCR/")

# Create output directories
outdir="processed_files"
dir.create(outdir,recursive = T,showWarnings = F)
plotdir <- 'plots'
dir.create(plotdir,showWarnings = F)
qcdir <-'qc'
dir.create(qcdir,showWarnings = F)

#Define sample names
samples=c("3005","4241","4901","5401","5791","6126","6499","8877","9480","HC-304","HC-581","HC-746")

files=list.files("rawdata/")
files=files[c(-1,-14)]

#Define input files

input10x= paste0("cellranger/",files,"/outs",sep="")
ndim <- 30

RunSoupX = function(input10x,files){
  filt <- Read10X(paste0(input10x,'/per_sample_outs/',files,"/count/sample_filtered_feature_bc_matrix/"))
  raw<- Read10X(paste0(input10x,'/multi/count/raw_feature_bc_matrix/'))
  sc = SoupChannel(raw,filt)
  s <- CreateSeuratObject(counts = filt) %>% 
    SCTransform(verbose = F) %>%
    RunPCA(verbose = F) %>%
    RunUMAP(dims = 1:30, verbose = F) %>% 
    FindNeighbors(dims = 1:30, verbose = F) %>%
    FindClusters(verbose = T)
  meta    <- s@meta.data
  umap    <- s@reductions$umap@cell.embeddings
  sc  <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
  sc  <- setDR(sc, umap)
  sc  <- autoEstCont(sc,forceAccept = F)
  adj.matrix  <- adjustCounts(sc, roundToInt = T)
  return(adj.matrix)
}

# Write function to process seurat obj

Processseurat=function(input10x,samples,files){
  ndim=30
  print(samples)
  print("Running SoupX")
  cts <-  RunSoupX(input10x,files)
  scrna <- CreateSeuratObject(cts,min.cells = 3, min.features = 200)
  #scrna=JoinLayers(scrna)
  print("QC")
  scrna= RunQC(scrna,org="human",dir=qcdir,filter=T,doubletdetection=F,UpperMitoCutoff=10,sample=samples)
  print("SCT")
  scrna <- SCTransform(scrna,vars.to.regress = c('percent.mito','nCount_RNA','nFeature_RNA','percent.ribo'),method = 'glmGamPoi',return.only.var.genes=F)
  print("Clustering")
  scrna <- RunPCA(scrna, features = VariableFeatures(object = scrna))
  scrna <- FindNeighbors(scrna, dims = 1:ndim)
  scrna <- FindClusters(scrna, resolution = seq(0.1,0.8,0.1))
  scrna <- RunUMAP(scrna, dims = 1:ndim)
  print("Computing doublet")
  scrna <- ComputeDoublets(scrna, ndim = 1:ndim, rate = 0.075, sct = TRUE)
  #Set res to 0.4 
  Idents(scrna) <- 'SCT_snn_res.0.4'
  scrna$orig_cluster <- Idents(scrna)
  scrna$orig.ident <-samples
  print("Saving data")
  #scrna@misc[["findallmarkers"]] <- FindAllMarkers(scrna, only.pos = TRUE)
  saveRDS(scrna,paste0(outdir,"/",samples,"_soupx_Seurat.RDS",sep=""))
}

processdata=function(samples){
  print("reading data")
  scrna=readRDS(paste0("processed_files/",samples,"_soupx_Seurat.RDS",sep=""))
  scrna[["percent.ribo"]] <- PercentageFeatureSet(scrna,pattern="^RP[LS]") 
  print("Computing doublet")
  scrna <- ComputeDoublets(scrna, ndim = 1:10, rate = 0.075, sct = TRUE)
  print("Saving data")
  saveRDS(scrna,paste0(outdir,"/",samples,"_soupx_Seurat_doublet.RDS",sep=""))
}
#Run all samps through the function
for(i in 1:length(samples)){
  Processseurat(input10x[i],samples[i],files[i])
}

#Run all samps through the function
for(i in 1:length(samples)){
  processdata(samples[i])
}

#### Integrate data ####
obj.list=list()
obj=list.files("processed_files/",pattern="_soupx_Seurat_doublet.RDS")
obj.list <- lapply(paste0("processed_files/",obj,sep=""), readRDS)
names(obj.list) <- gsub("_soupx_Seurat_doublet.RDS","",obj)

#Add metadata
scrna <- merge(obj.list[[1]],obj.list[2:12])
meta=read.csv("data/metadata.csv")
meta= meta %>% select(sample_id,timepoint,condition,final_mods) 


fullmeta <- scrna@meta.data %>% rownames_to_column("ID") %>% inner_join(., meta,by=c('orig.ident'='sample_id')) %>% column_to_rownames("ID")
scrna <- AddMetaData(scrna,fullmeta)

scrna[["RNA"]] <- split(scrna[["RNA"]], f = scrna$orig.ident)

#Normalize, scale and pca data
DefaultAssay(scrna) <- 'RNA'
scrna <- NormalizeData(scrna)
scrna <- FindVariableFeatures(scrna)
scrna <- ScaleData(scrna, vars.to.regress = c('percent.mito','nCount_RNA','nFeature_RNA','percent.ribo'))
scrna <- RunPCA(scrna)

#Integrate data by RPCA

scrna <- IntegrateLayers(
  object = scrna, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = FALSE
)

scrna <- FindNeighbors(scrna, reduction = "integrated.rpca", dims = 1:30)
scrna <- FindClusters(scrna, resolution = 0.5, cluster.name = "rpca_clusters")
scrna <- RunUMAP(scrna, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
scrna@meta.data$rpca_clusters= as.numeric(as.character(scrna@meta.data$rpca_clusters))

saveRDS(scrna,"Seurat/Final_integration_soupx.RDS")


#### Annotate data ####
##Run Azimuth
DefaultAssay(scrna)="RNA"
#Read Reference
reference= readRDS("pbmc_multimodal_2023.rds")

anchors <- FindTransferAnchors(
  reference = reference,
  query = scrna,
  normalization.method = "SCT",
  reference.reduction = "pca",
  dims = 1:50,
  recompute.residuals=FALSE
)

scrna <- MapQuery(
  anchorset = anchors,
  query = scrna,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    celltype.l3 = "celltype.l3"
  ),
  reference.reduction = "pca"
)

p1=DimPlot(scrna,group.by = "rpca_clusters",reduction="umap.rpca",label=T) 
p2=DimPlot(scrna,group.by = "predicted.celltype.l1",reduction="umap.rpca",label=T,label.size = 4) +NoLegend()
p3=DimPlot(scrna,group.by = "predicted.celltype.l2",reduction="umap.rpca",label=T,label.size = 4) +NoLegend()
p4=DimPlot(scrna,group.by = "predicted.celltype.l3",reduction="umap.rpca",label=T,label.size = 4) +NoLegend()
p=(p1+p2)/(p3+p4)
ggsave(p,filename="plots/RPCA_int_annotated.png",width=15,height=10,units="in",dpi=300)

#Use azimuth annotations as first pass to identify the largere cell types (Bcells, T cells, Monocytes)
### subset other celltypes to annotate ####
#### BCELLS ####
bcell_clusts=c(0,5,13,14,15,16)
sub.b=subset(scrna,idents = bcell_clusts)

sub.b <- RunPCA(sub.b, npcs = 50)
ElbowPlot(scrna, ndims = 50)
dims <- 1:10
sub.b <- FindNeighbors(object = sub.b,dims=dims,k.param = 20)
sub.b <- FindClusters(object = sub.b,resolution=0.3,cluster.name = "reclus_clusters")
sub.b <- RunUMAP(object = sub.b, reduction = "pca", n.neighbors = 25,n.components = 2,dims = dims,min.dist=0.3,reduction.name="reclus_umap")
DimPlot(sub.b,group.by = "reclus_clusters",label=T,reduction = "reclus_umap")
sub.b$var_cluster <- sub.b$reclus_clusters

library(HGNChelper)
#Annotate using scType
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# get cell-type-specific gene sets from our in-built database (DB)
gs_list <- gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R"); 
Idents(sub.b) = "var_cluster"
sub.b <- run_sctype(sub.b, assay = "RNA", scaled = TRUE, known_tissue_type="Immune system",custom_marker_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx",name="sctype_classification")
DimPlot(sub.b,group.by = "sctype_classification",label=T,reduction = "reclus_umap") +DimPlot(sub.b,group.by = "reclus_clusters",label=T,reduction = "reclus_umap")

FeaturePlot(sub.b,features=c("CD19","CD22","MS4A1","IGHG1","IGHA1","IGHE","IGHM","IL4R","CD38"),reduction = "reclus_umap") 
DotPlot(sub.b,features=c("CD19","CD22","MS4A1","IGHG1","IGHA1","IGHE","IGHM","IL4R","CD38"),group.by = 'reclus_clusters')
FeaturePlot(sub.b,features=c("ITGAX","CD83","CD34","HLA-DRA","CD1A"),reduction = "reclus_umap") 
DotPlot(sub.b,features=c("ITGAX","CD83","CD34","HLA-DRA","CD1A"),group.by = 'reclus_clusters')
FeaturePlot(sub.b,features=c("TNFRSF17","PRDM1","CD27","CD38","SDC1","CXCR4"),reduction = "reclus_umap") 
DotPlot(sub.b,features=c("TNFRSF17","PRDM1","CD27","CD38","SDC1","CXCR4"),group.by = 'reclus_clusters') #plasma
FeaturePlot(sub.b,features=c("CD19","MS4A1","CD27","CR2","TNFRSF13C"),reduction = "reclus_umap") 
DotPlot(sub.b,features=c("CD19","MS4A1","CD27","CR2","TNFRSF13C"),group.by = 'reclus_clusters') #MEMORY
FeaturePlot(sub.b,features=c("CD19","MS4A1","TNF","IFNG","CSF2","IL2","IL4","IL17A","IL6"),reduction = "reclus_umap") 
DotPlot(sub.b,features=c("CD19","MS4A1","TNF","IFNG","CSF2","IL2","IL4","IL17A","IL6"),group.by = 'reclus_clusters') #EFFECTOR B
FeaturePlot(sub.b,features=c("CD19","MS4A1","CR2","FCER2","CD24","CD38","CD5"),reduction = "reclus_umap") 
DotPlot(sub.b,features=c("CD19","MS4A1","CR2","FCER2","CD24","CD38","CD5"),group.by = 'reclus_clusters') #MATURE B
FeaturePlot(sub.b,features=c("MME","CD19","MS4A1","CR2","CD24","CD38","CD93","IL4R"),reduction = "reclus_umap") 
DotPlot(sub.b,features=c("MME","CD19","MS4A1","CR2","CD24","CD38","CD93","IL4R"),group.by = 'reclus_clusters') #IMMATURE B

DimPlot(sub.b,group.by = "sctype_classification",label=T,reduction = "reclus_umap") +DimPlot(sub.b,group.by = "reclus_clusters",label=T,reduction = "reclus_umap")
deg_bcells<- FindAllMarkers(object = sub.b, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>% rownames_to_column("Gene")

sub.b$possible_anno=ifelse(sub.b$reclus_clusters %in% c("0","1","2","3"),"Naive B cells",ifelse(sub.b$reclus_clusters =="4","Immature B cell",ifelse(sub.b$reclus_clusters %in% c("5","6"),"Memory B cell","Plasma B cell")))
DimPlot(sub.b,reduction = "reclus_umap",group.by = "possible_anno")

scrna=subset(scrna,percent.mito<5)
saveRDS(scrna,"Seurat/Final_integration_mitofilter.RDS")

####MONOCYTES####
mono_clusts=c(7,11,17)
sub.mono=subset(scrna,idents = mono_clusts)
#Remove t-cell cluster and recluster
sub.mono=subset(sub.mono,CD3E>0,invert=T)
sub.mono2=subset(sub.mono,CD3E>0)
sub.mono <- RunPCA(sub.mono, npcs = 50)
ElbowPlot(scrna, ndims = 50)
dims <- 1:10
sub.mono <- FindNeighbors(object = sub.mono,dims=dims,k.param = 20)
sub.mono <- FindClusters(object = sub.mono,resolution=0.3,cluster.name = "reclus_clusters")
sub.mono <- RunUMAP(object = sub.mono, reduction = "pca", n.neighbors = 25,n.components = 2,dims = dims,min.dist=0.3,reduction.name="reclus_umap")
DimPlot(sub.mono,group.by = "reclus_clusters",label=T,reduction = "reclus_umap")
sub.mono$var_cluster <- sub.mono$reclus_clusters
Idents(sub.mono) = "var_cluster"
sub.mono <- run_sctype(sub.mono, assay = "RNA", scaled = TRUE, known_tissue_type="Immune system",custom_marker_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx",name="sctype_classification")
DimPlot(sub.mono,group.by = "sctype_classification",label=T,reduction = "reclus_umap") +DimPlot(sub.mono,group.by = "reclus_clusters",label=T,reduction = "reclus_umap")

FeaturePlot(sub.mono,features=c("CD14","FCGR3A","NKG7","NCAM1","SELL","CCR2","CCR5","ITGAX","HLA-DRA"),reduction = "reclus_umap") 
DotPlot(sub.mono,features=c("CD14","FCGR3A","NKG7","NCAM1","SELL","CCR2","CCR5","ITGAX","HLA-DRA"),group.by = 'reclus_clusters')
FeaturePlot(sub.mono,features=c("IFNG","CD83","CD34","CD38","CD27","SDC1","TNFRSF17","PRDM1","CDCR4"),reduction = "reclus_umap") 
DotPlot(sub.mono,features=c("ITGAX","CD83","CD34","HLA-DRA","CD1A"),group.by = 'reclus_clusters')

deg_mono<- FindAllMarkers(object = sub.mono, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>% rownames_to_column("Gene")

sub.mono$possible_anno=ifelse(sub.mono$reclus_clusters %in% c("0","1","4","5","6"),"Classical monocytes",ifelse(sub.mono$reclus_clusters =="3","Non Classical monocytes",ifelse(sub.mono$reclus_clusters =="7","Unknown monocyte 6","Naive CD8+ T Cells")))
DimPlot(sub.mono,reduction = "reclus_umap",group.by = "possible_anno",label=T)

sub.mono2=subset(sub.mono, CD3E>0)
sub.mono2$possible_anno="Naive CD8+ T Cells"

#### Subset tcells ####
#Specify T-cell clusters to subset
tcell_clusts=c(1,2,3,4,6,8,9,10,12,18,19)

#subset and recluster
sub=subset(scrna,idents = tcell_clusts)

#Recluster
sub <- RunPCA(sub, npcs = 50)
ElbowPlot(scrna, ndims = 50)
dims <- 1:10
sub <- FindNeighbors(object = sub,dims=dims,k.param = 20)
sub <- FindClusters(object = sub,resolution=0.4,cluster.name = "reclus_clusters")
sub <- RunUMAP(object = sub, reduction = "pca", n.neighbors = 25,n.components = 2,dims = dims,min.dist=0.3,reduction.name="reclus_umap")
DimPlot(sub,group.by = "reclus_clusters",label=T,reduction = "reclus_umap") %>% ggsave(filename="plots/Tcellsubset_umap.png", width = 8, height = 7, units = "in")
sub$var_cluster <- sub$reclus_clusters
Idents(sub) = "var_cluster"

#Find top markers
sub=JoinLayers(sub)
sub@misc[["findallmarkers"]] <- FindAllMarkers(object = sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
deg=sub@misc[["findallmarkers"]]
#get top50 markers
d=deg %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50) %>% select(gene,p_val:cluster)

#Annotate using scType
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# get cell-type-specific gene sets from our in-built database (DB)
gs_list <- gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R"); 
sub <- run_sctype(sub, assay = "RNA", scaled = TRUE, known_tissue_type="Immune system",custom_marker_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx",name="sctype_classification")

p4=DimPlot(sub,group.by = "sctype_classification",label=T,label.size = 4,reduction="reclus_umap",repel = T) +NoLegend()
ggsave(p4,filename="plots/anno.png",width=15,height=10,units="in",dpi=300)
saveRDS(sub,"Seurat/Tcellsubset_soupx_annotated.RDS")

#### Project back to integrated data (Load tcell annotations from tcell subset)####

meta=scrna@meta.data
meta.t =sub@meta.data %>% select("celltype") %>% rename("possible_anno"="celltype")
meta.b= sub.b@meta.data %>% select(possible_anno)
meta.mono= sub.mono@meta.data %>% select("sctype_classification") %>% rename("possible_anno"="sctype_classification")
meta.mono2= sub.mono2@meta.data %>% select("possible_anno") 

all_anno=rbind(meta.t,meta.b,meta.mono,meta.mono2)
all_anno =all_anno %>% rownames_to_column("id")
meta= meta %>% rownames_to_column("id") %>% select(-possible_anno) %>% left_join(.,all_anno,by="id") %>% column_to_rownames("id")

scrna@meta.data=meta
scrna$celltype=scrna@meta.data$possible_anno
DimPlot(scrna,reduction = "umap.rpca",group.by = "celltype",label=T)

saveRDS(scrna,"Seurat/Final_integration_soupx_annotated.RDS")

