library(Seurat)
library(tidyverse)
library(HGNChelper)
library(openxlsx)
library(readxl)
library(dplyr)
library(RColorBrewer)
library(scCustomize)
library(ggsci)
library(scplotter)
library(readxl)
library(qs)
library(ggsignif)
library(patchwork)
library(ggpubr)
library(qs)
library(UCell)
library(data.table)
library(nichenetr)

dir.create("plots/paper_figures/")
setwd("/lindell_lab/Data/Sepsis_scRNAseq_TCR/")
tcells <- readRDS("Seurat/Tcellsubset_soupx_annotated.RDS")
scrna=qread("Final_integration_soupx_annotated.qs")

cpallette=c("#FFFF99","#1F78B4","#00CCCC","#9367BC", "#33A02C","#666666","#990000","#FF7E0E","#FDBF6F","#B2DF8A","#FB9A99","#A6CEE3","#33FF99","#CAB2D6","#D7C1B1")
names(cpallette)=sort(unique(scrna$celltype))

##### 6B #####
Idents(scrna)="condition"
sub=subset(scrna,idents=c("Healthy","Sepsis"))

#### Downsample Sepsis data for umap ####
Idents(scrna)="condition"

Healthy_cells <- WhichCells(scrna, idents = "Healthy Control")
Sepsis_cells <- WhichCells(scrna, idents= "MODS")

# Total target number of cells for sepsis
target_sepsis <- length(Healthy_cells)

# Proportion of cell types in Sepsis
prop <- table(scrna$celltype[Sepsis_cells])/length(Sepsis_cells)

# Calculate target number of cells per cell type
target_counts <- round(prop * target_sepsis)

downsampled_Sepsis <- unlist(lapply(names(target_counts), function(cell_type) {
  cell_ids <- WhichCells(scrna, idents = "MODS")
  cell_ids_filtered <- cell_ids[scrna$celltype[cell_ids] == cell_type]
  sample(cell_ids_filtered, size = target_counts[cell_type], replace = FALSE)
}))

# Combine downsampled Healthy cells with all Sepsis cells
downsampled_cells <- c(downsampled_Sepsis, Healthy_cells)

# Subset the Seurat object
scrna_downsampled <- subset(scrna, cells = downsampled_cells)

pdf("plots/paper_figures/Final_umap_downsampled.pdf",width=10,height=5)
DimPlot(scrna_downsampled,reduction="umap.rpca",group.by = "celltype",split.by = "condition",cols =  cpallette,label=T,repel=T,label.size=2)+ theme(legend.text=element_text(size=5))
dev.off()

#### 6C ####
sub=subset(scrna,idents =sort(unique(scrna$celltype))[c(1:3,5:6,10:11,13)] ,invert=T) #non T-cells
sub_t=subset(scrna,idents=sort(unique(scrna$celltype))[c(1:3,5:6,10:11,13)],invert=F) #tcells
sub_t$celltype=factor(sub_t$celltype,levels=unique(sub_t$celltype)[c(2,1,3,5,4,8,6,7)])
sub$celltype=factor(sub$celltype,levels=c("Naive B cells","Immature B cell","Memory B cell","Plasma B cell","NK cells","Classical Monocytes","Non-classical monocytes")) #order them

p1=sub@meta.data %>% 
  group_by(condition,celltype) %>%
  summarise(n = n()) %>%
  mutate(pct = n / sum(n)) %>% 
  ggplot(.,aes(x=condition,y=pct,fill=celltype)) +  
  geom_bar(position="fill", stat="identity")  + 
  scale_y_continuous(labels = scales::percent_format()) + labs(y="Proportion of CD3- Cells")+
  scale_fill_manual(values = cpallette)+ theme_bw() + theme(legend.position = 'bottom',legend.text=element_text(size=4)) 

p2=sub_t@meta.data %>% 
  group_by(condition,celltype) %>%
  summarise(n = n()) %>%
  mutate(pct = n / sum(n)) %>% 
  ggplot(.,aes(x=condition,y=pct,fill=celltype)) +  
  geom_bar(position="fill", stat="identity")  + 
  scale_y_continuous(labels = scales::percent_format()) + labs(y="Proportion of CD3+ Cells")+
  scale_fill_manual(values = cpallette)+ theme_bw() + theme(legend.position = 'bottom',legend.text=element_text(size=4))+
  guides(fill=guide_legend(nrow=3,byrow=TRUE))

pdf("plots/paper_figures/5C_barplot.pdf",width=8,height=10)
p1+p2
dev.off()

#### UCell scoring 6D and 6E ####
#Calculate module score
library(UCell)
my_protective=c("ZDHHC17", "FAS", "GK", "ICAM3","MME", "PDE4B", "PIK3CD","PTEN", "RAF1", "TLR1","PPP1R12A", "MAPK14", "SOS2",
                "TXN", "ASAH1", "ATG3", "BCAT1","BCL2L11", "BTK", "BTN2A2","CASP1", "CCL2", "CREB1","EP300", "GNAI3", "IL1A", "JAK2",
                "MAFB", "MAP3K1", "MAP3K3","PAK2", "PLEKHO1", "POU2F2","PRKAR1A", "PRKCB", "RHBDF2","SEMA6B", "SP1", "TLE4",
                "BMPR2", "CTNNB1", "INPP5D","ITGAV", "SLC12A7", "TBK1","VAMP5", "VRK2", "YKT6")
my_detrimental=c("ANXA3", "ARG1", "AZU1", "CAMP","CEACAM8", "CEP55", "CRISP2","CTSG", "DEFA4", "GADD45A",
                 "HMMR", "KIF15", "LCN2", "LTF","OLFM4", "ORM1", "PRC1", "SLPI","STOM", "AQP9", "BCL6", "KLHL2",
                 "PPL", "HTRA1", "TYK2", "SLC1A5","STX1A")
lym_protective=c("ARL14EP", "BPGM", "BTN3A2","BUB3", "CAMK4", "CASP8","CCNB1IP1", "CD247", "CD3E","CD3G", "DBT", "DDX6", "DYRK2",
                 "JAK1", "KLRB1", "MAP4K1","NCR3", "PIK3R1", "PLCG1","PPP2R5C", "SEMA4F", "SIDT1","SMAD4", "SMYD2", "TP53BP1",
                 "TRIB2", "ZAP70", "ZCCHC4","ZNF831")

scrna <- AddModuleScore_UCell(scrna, features=list(my_protective), name="Myeloid_protective_Score")
scrna <- SmoothKNN(scrna,signature.names="signature_1Myeloid_protective_Score")
scrna <- AddModuleScore_UCell(scrna, features=list(my_detrimental), name="Myeloid_detrimental_Score")
scrna <- SmoothKNN(scrna,signature.names="signature_1Myeloid_detrimental_Score")
scrna <- AddModuleScore_UCell(scrna, features=list(lym_protective), name="Lymphoid_protective_Score")
scrna <- SmoothKNN(scrna,signature.names="signature_1Lymphoid_protective_Score")

pdf(paste0("plots/paper_figures/lymphoid_prot_splitvln_fulldata_knn.pdf",sep=""),width=12,height=8)
VlnPlot(scrna,features="signature_1Lymphoid_protective_Score_kNN", group.by = "celltype",split.by = "condition",pt.size = 0,cols=brewer.pal(n=8,name = "Accent"),split.plot = T)
dev.off()
pdf(paste0("plots/paper_figures/reactome/lymphoid_prot_box_fulldata_knn.pdf",sep=""),width=12,height=10)
SCpubr::do_BoxPlot(sample = scrna,feature = "signature_1Lymphoid_protective_Score_kNN",group.by = "celltype_condition",use_test = F, split.by="condition",map_signif_level = T,test = "wilcox.text")+stat_compare_means(comparisons = comparisons)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


comparisons=list()
for(i in 1:length(unique(scrna$celltype))){
  comparisons[[i]]=c(paste0(unique(scrna$celltype),"_Healthy Control",sep="")[i],paste0(unique(scrna$celltype),"_MODS",sep="")[i])
}
comparisons


#### Module score KEGG  6G####
BiocManager::install("KEGGREST")
library(KEGGREST)

library("msigdbr")
c_gene_sets = msigdbr(species = "human", category = "C2")
c_gene_sets=c_gene_sets[c_gene_sets$gs_subcat=="CP:KEGG_LEGACY",]

mtor_genes=c_gene_sets$gene_symbol[c_gene_sets$gs_name=="KEGG_MTOR_SIGNALING_PATHWAY"]
oxp_genes=c_gene_sets$gene_symbol[c_gene_sets$gs_name=="KEGG_OXIDATIVE_PHOSPHORYLATION"]
gly_genes=c_gene_sets$gene_symbol[c_gene_sets$gs_name=="KEGG_GLYCOLYSIS_GLUCONEOGENESIS"]
fa_genes=c_gene_sets$gene_symbol[c_gene_sets$gs_name=="KEGG_FATTY_ACID_METABOLISM"]
purine_genes=c_gene_sets$gene_symbol[c_gene_sets$gs_name=="KEGG_PURINE_METABOLISM"]
pyruvate_genes=c_gene_sets$gene_symbol[c_gene_sets$gs_name=="KEGG_PYRUVATE_METABOLISM"]
genelist=list(mtor_genes=mtor_genes,oxp_genes=oxp_genes,gly_genes=gly_genes,fa_genes=fa_genes,purine_genes=purine_genes,pyruvate_genes=pyruvate_genes)

Idents(tcells)="celltype"
tcells_sub=subset(tcells,idents=c("CD8+ CM cells","CD8+ EM1 cells","CD8+ EMRA cells","Naive CD8+ T cells"))
Idents(tcells)="celltype_condition"
for(i in names(genelist)){
  tcells <- AddModuleScore_UCell(tcells, features=list(eval(parse(text = i))), name=gsub("genes","Ucell_Score",i))
  tcells <- SmoothKNN(tcells,signature.names=paste0("signature_1",gsub("genes","Ucell_Score",i),sep=""))
  
  g1=SCpubr::do_BoxPlot(sample = tcells_sub,feature = paste0("signature_1",gsub("genes","Ucell_Score",i),"_kNN",sep=""),group.by = "celltype_condition",use_test = F, split.by="condition",map_signif_level = T,test = "wilcox.text")+stat_compare_means(comparisons = comparisons)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(g1,filename=paste0("plots/paper_figures/kegg/",gsub("genes","Ucell_Score",i),"_box_knn.pdf",sep=""),width=12,height=10)
  
  g2=VlnPlot(tcells_sub,features =paste0("signature_1",gsub("genes","Ucell_Score",i),"_kNN",sep=""), group.by = "celltype",split.by = "condition",pt.size = 0,cols=brewer.pal(n=8,name = "Accent"),split.plot = T)
  ggsave(g2,filename=paste0("plots/paper_figures/kegg/",gsub("genes","Ucell_Score",i),"_vln_knn.pdf",sep=""),width=12,height=10)
  
}




#### Module scores from immune dictionary 6F####
full=full[full$cyt %in% c("OSM","TGF-beta-1","TNFa","IFNb","IFNg","IFNa1","TL1A","G-CSF","IL1a","IL1b","IL6","IL10","IL2","IL4","IL7","IL12","IL13","IL15","IL18","IL21","IL22","IL27","IL33","IL36a"),]
df=full[full$cell_type=="CD8+ T cell",]
df=df[df$p_val_adj<0.05,]
df=df[abs(df$avg_log2FC)>1,]
gl=c("TGF-beta-1","TNFa","IFNb","IFNg","IFNa1","TL1A","G-CSF","IL1a","IL1b","IL6","IL10","IL2","IL4","IL7","IL12","IL13","IL15","IL18","IL21","IL22","IL27","IL33","IL36a")

getlist=function(name){
  genes=df$gene[df$cyt == name]
  genes <- convert_mouse_to_human_symbols(genes)
  names(genes)=NULL
  genes=genes[!is.na(genes)]
}
genelist=list()
for(i in gl){
  genelist[[i]]=getlist(i)
}
tcells$celltype= factor(tcells$celltype,levels=unique(tcells$celltype)[c(2,1,3:5,8,7,6)])

for(i in names(genelist)){
  tcells <- AddModuleScore_UCell(tcells, features=list(genelist[[i]]), name=gsub("genes","Ucell_Score",i))
  tcells <- SmoothKNN(tcells,signature.names=paste0("signature_1",gsub("genes","Ucell_Score",i),sep=""))
  
  pdf(paste0("plots/paper_figures/imm_dic_perturb/",gsub("genes","Ucell_Score",i),"_splitvln_fulldata_knn.pdf",sep=""),width=12,height=8)
  VlnPlot(tcells,features =paste0("signature_1",gsub("genes","Ucell_Score",i),"_kNN",sep=""), group.by = "celltype",split.by = "condition",pt.size = 0,cols=brewer.pal(n=8,name = "Accent"),split.plot = T)
  dev.off()
  
  pdf(paste0("plots/paper_figures/imm_dic_perturb/",gsub("genes","Ucell_Score",i),"_box_fulldata_knn.pdf",sep=""),width=12,height=10)
  SCpubr::do_BoxPlot(sample = scrna,feature = paste0("signature_1",gsub("genes","Ucell_Score",i),"_kNN",sep=""),group.by = "celltype_condition",use_test = F, split.by="condition",map_signif_level = T,test = "wilcox.text")+stat_compare_means(comparisons = comparisons)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
}
saveRDS(tcells,"tcells_cytokine_scores.RDS")

#### Module score dotplot ####
gl=c("TNFa","IFNb","IFNg","IFNa1","TL1A","IL1b","IL6","IL10","IL2","IL4","IL7","IL12","IL22","IL27","IL33","IL36a","OSM")
ct=levels(tcells@meta.data[[group_col]])[c(3:6)]
getlist=function(name){
  genes=df$gene[df$cyt == name]
  genes <- convert_mouse_to_human_symbols(genes)
  names(genes)=NULL
  genes=genes[!is.na(genes)]
}
genelist=list()
for(i in gl){
  genelist[[i]]=getlist(i)
}
# Setup: metadata column names
group_col <- "celltype"      # column for cell type or cluster
condition_col <- "condition" # e.g. 'treated' vs 'baseline'
baseline_value <- "Healthy Control" # reference group

# Initialize results list
bubble_data <- list()

# Loop through gene sets
for (module in names(genelist)) {
  score_col <- paste0("signature_1",gsub("genes","Ucell_Score",module),sep="")
  
  # Loop through each cell type
  for (c in ct) {
    sub <- tcells@meta.data %>% filter(celltype == c)
    
    # Check if both groups are present
    conds <- unique(sub[[condition_col]])
    if (all(c(baseline_value, setdiff(conds, baseline_value)) %in% conds)) {
      
      baseline_scores <- sub %>% filter(.data[[condition_col]] == baseline_value) %>% pull(score_col)
      other_scores <- sub %>% filter(.data[[condition_col]] != baseline_value) %>% pull(score_col)
      
      # Log2 Fold Change
      median_baseline <- median(baseline_scores, na.rm = TRUE)
      median_other <- median(other_scores, na.rm = TRUE)
      log2fc <- log2((median_other + 1e-6) / (median_baseline + 1e-6))  # add small offset to avoid div by zero
      
      # Wilcoxon Test
      pval <- tryCatch(
        wilcox.test(other_scores, baseline_scores)$p.value,
        error = function(e) NA
      )
      
      bubble_data[[length(bubble_data) + 1]] <- tibble(
        Module = module,
        CellType = c,
        Log2FC = log2fc,
        Pval = pval
      )
    }
  }
}

# Combine results
df <- bind_rows(bubble_data) %>%
  mutate(NegLog10P = log10(Pval))
df$Cyt_group="unassigned"
df$Cyt_group =ifelse(df$Module %in% c("IL1b","IL33","IL36a"),"G1",df$Cyt_group)
df$Cyt_group =ifelse(df$Module %in% c("IL2","IL4","IL7"),"G2",df$Cyt_group)
df$Cyt_group =ifelse(df$Module %in% c("IL6","IL12","IL27","OSM"),"G3",df$Cyt_group)
df$Cyt_group =ifelse(df$Module %in% c("IL10","IL22"),"G4",df$Cyt_group)
df$Cyt_group =ifelse(df$Module %in% c("IFNa1","IFNb","IFNg"),"G5",df$Cyt_group)
df$Cyt_group =ifelse(df$Module %in% c("TNFa","TL1A"),"G6",df$Cyt_group)
df$CellType=factor(df$CellType,levels = unique(df$CellType)[c(1:4)])

my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
             "#E16462FF", "#FCA636FF", "#F0F921FF")
g3=ggballoonplot(
  df, x = "CellType", y = "Module",
  fill = "Log2FC", size = "NegLog10P",
  ggtheme = theme_classic()
)+facet_grid(Cyt_group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradientn(colors = my_cols,
                       values = scales::rescale(c(-0.5,-0.1,0,0.1, 0.5)),  # position of color breaks
                       limits = c(-0.5, 0.5),   )+
  labs(
    title = "Module Score Bubble Plot",
    x = "Cell Type",
    y = "Cytokines",
    size = "-log10(p-value)",
    color = "log2 Fold Change"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(g3,filename="plots/paper_figures/cytokine_celltype_dotplot_scaledfc.pdf",width=10,height = 10,dpi = 300)
