library("scRepertoire")
library(circlize)
library(scales)
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

setwd("/lindell_lab/Data/Sepsis_scRNAseq_TCR/")

files=list.files()
files=files[c(-1,-14)]

#Define sample names
samples=c("3005","4241","4901","5401","5791","6126","6499","8877","9480","HC-304","HC-581","HC-746")

#Define input files
vdj= paste0(files,"/outs/per_sample_outs/",files,"/vdj_t/filtered_contig_annotations.csv",sep="")
contig.list=list()
names(contig.list)=samples

#load data
for (i in 1:length(samples)){
  contig.list[[samples[i]]]=read.csv(paste0(files[i],"/outs/per_sample_outs/",files[i],"/vdj_t/filtered_contig_annotations.csv",sep=""))
}

for(i in seq_along(contig.list)) {
  contig.list[[i]]$barcode <- paste0(contig.list[[i]]$barcode, "_", i)
}

#combine contigs to clones
combined.TCR <- combineTCR(contig.list, 
                           samples = samples,
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = FALSE)

head(combined.TCR[[1]])

#Add experimental group information to the combined TCRData object
combined.TCR <- addVariable(combined.TCR, 
                            variable.name = "Group", 
                            variables = c("Sepsis","Sepsis","Sepsis","Non-sepsis","Non-sepsis","Sepsis","Sepsis","Non-sepsis","Sepsis","Healthy","Healthy","Healthy"))


#Export clonal information as csv file
exportClones(combined.TCR, 
             write.file = TRUE,
             dir = "results/",
             file.name = "clones_new.csv")



#return the total or relative numbers of unique clones per group and per sample
pdf("plots/num_clones_per_group.pdf",width=5,height=5)
clonalQuant(combined.TCR, cloneCall = "strict", group.by = "Group", scale = TRUE)
dev.off()

pdf("plots/num_clones_per_samp.pdf",width=5,height=5)
clonalQuant(combined.TCR,cloneCall="strict", chain = "both", scale = TRUE)             
dev.off()


#compare how proportions of certain clones compare between groups. 
png("plots/clonalcompare_aa.png",width=15,height=15,units = "in",res=300)
clonalCompare(combined.TCR,group.by = "Group", top.clones = 10, cloneCall="aa", graph = "alluvial")
dev.off()

#clonal space occupied by clonotypes of specific proportions (clone size -number of tcells making up the clone)
png("plots/clonal_homeostasis_aa_group.png",width=5,height=5,units="in",res=300)
clonalHomeostasis(combined.TCR, group.by = "Group",cloneCall = "aa")
dev.off()


#### Clonal diversity ####
c1=clonalDiversity(combined.TCR, cloneCall = "gene",x.axis = "sample") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
c2=clonalDiversity(combined.TCR, cloneCall = "aa",x.axis = "sample") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
c3=clonalDiversity(combined.TCR, cloneCall = "strict",x.axis = "sample") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
png("plots/clonaldiversity_allclonecalls.png",width = 15,height = 15,units = "in",res=300)
c1/c2/c3
dev.off()

#Export Diversity indices
t1=clonalDiversity(combined.TCR, cloneCall = "gene",x.axis = "sample",exportTable = T)
t2=clonalDiversity(combined.TCR, cloneCall = "aa",x.axis = "sample",exportTable = T)
t3=clonalDiversity(combined.TCR, cloneCall = "strict",x.axis = "sample",exportTable = T)
list=list(Gene=t1,AA=t2,gene_nt=t3)
writexl::write_xlsx(list,"results/clonal_diversity_index.xlsx")

t1=clonalDiversity(combined.TCR, cloneCall = "gene",x.axis = "sample",group.by = "Group",exportTable = T)
t2=clonalDiversity(combined.TCR, cloneCall = "aa",x.axis = "sample",group.by = "Group",exportTable = T)
t3=clonalDiversity(combined.TCR, cloneCall = "strict",x.axis = "sample",group.by = "Group",exportTable = T)
list=list(Gene=t1,AA=t2,gene_nt=t3)
writexl::write_xlsx(list,"results/clonal_diversity_index_bygroup.xlsx")


#### Combine with seurat ####
scrna <- readRDS("../Analysis/Sepsis_scRNAseq/Seurat/Tcellsubset_soupx_annotated.RDS")
scrna@meta.data$sample=scrna@meta.data$orig.ident
scrna$celltype_condition=paste0(scrna$celltype,"_",scrna$final_mods,sep="")

#make sure rownames of seurat metadata and tcr data match
scrna=RenameCells(scrna, new.names = paste0(scrna@meta.data$orig.ident,"_",rownames(scrna@meta.data),sep=""))

#Combine Expression
scrna <- combineExpression(combined.TCR, scrna, cloneCall="gene", group.by = "sample", proportion = TRUE)
scrna@meta.data$cdr3s_pat <- paste0(scrna@meta.data$CTaa, "_",scrna@meta.data$sample,sep="")


#### Remove non-sepsis ####
#### Clonal diversity by celltype - subset ####
#Define sample names
need.samps=c("3005","4241","4901","6126","6499","9480","HC-304","HC-581","HC-746")
Idents(scrna)="orig.ident"
sub=subset(scrna,idents=need.samps)
sub$celltype_cond=paste0(sub$sctype_classification,"_",sub$final_mods,sep="")

pdf("plots/shannon_gene_bycelltype_sub.pdf",width=9,height=6)
ClonalDiversityPlot(sub, group_by = "sctype_classification", plot_type = "box")
dev.off()

pdf("plots/shannon_gene_bycelltype_cond_sub.pdf",width=9,height=6)
ClonalDiversityPlot(sub, group_by = "celltype_cond", plot_type = "box")
dev.off()

###export diversity indices
t2=clonalDiversity(sub, cloneCall = "aa",x.axis = "final_mods",group.by = "celltype_condition",exportTable = T)
writexl::write_xlsx(t2,"results/clonal_diversity_index_subset.xlsx")

###Chord diagram to 
circles <- getCirclize(scrna,group.by = "celltype")
grid.cols <- hue_pal()(length(unique(scrna$celltype)))
names(grid.cols) <- unique(scrna$celltype)

#Graphing the chord diagram
png("plots/chord.png",width=8,height=8,res=300,units="in")
chordDiagram(circles, self.link = 1, grid.col = grid.cols)
dev.off()

####TRex Annotation ####
library("Trex")
scrna_anno.a <- annotateDB(scrna,chains = "TRA")
scrna_anno.b <- annotateDB(scrna,chains = "TRB")

#### We got top 500 clones in the data ####
Idents(scrna)="condition"
topclones=clonalCompare(scrna ,top.clones = 500, samples=c("Healthy Control","MODS"),cloneCall="aa",graph = "alluvial",exportTable = T)
topclones_healthy=topclones[topclones$Sample=="Healthy Control",]
topclones_sepsis=topclones[topclones$Sample=="MODS",]

#Map top500 sepsis clones to annotated scrna looking at both chains A and B
sepsis.a=scrna_anno.a@meta.data[scrna_anno.a$CTaa %in% topclones_sepsis$clones,] %>% select(orig.ident,final_mods,rpca_clusters,celltype,celltype_condition:TRA_Database)
sepsis.b=scrna_anno.b@meta.data[scrna_anno.b$CTaa %in% topclones_sepsis$clones,] %>% select(orig.ident,final_mods,rpca_clusters,celltype,celltype_condition:TRB_Database)

#Extract out cells mapped to CMV
cmv.cells=sepsis.a %>% rownames_to_column('id') %>% filter(grepl('CMV', TRA_Epitope.species)) %>% select(id,celltype)
noncmv.cells=sepsis.a %>% rownames_to_column('id') %>% filter(!grepl('CMV', TRA_Epitope.species)) %>% drop_na() %>% select(id)
noncmv.cells2=sepsis.a %>% rownames_to_column('id') %>% filter(!grepl('CMV', TRA_Epitope.species))  %>% select(id)

#Subset only CD8
tcells=subset(scrna,idents=c("Naive CD8+ T cells","CD8+ EM1 cells","CD8+ CM cells","CD8+ EMRA cells"))
tcells$id=paste0(tcells$orig.ident,"_",rownames(tcells@meta.data),sep="")
tcells$cmv.anno2=ifelse(tcells$id %in% cmv.cells$id,"CMV",ifelse(tcells$id %in% noncmv.cells$id,"Non-CMV","Not Annotated"))

####plot Ucell####
pdf("plots/paper_figures/exhaustion_score_ucell.pdf",width=10,height=7)
SCpubr::do_BoxPlot(sample = tcells,group.by = "cmv.anno2",feature = "signature_1Exhaustion_Ucell",use_test = TRUE,comparisons = list(c("CMV", "Non-CMV")),map_signif_level = F)
dev.off()

pdf("plots/paper_figures/glycolysis_score_ucell.pdf",width=10,height=7)
SCpubr::do_BoxPlot(sample = tcells,group.by = "cmv.anno2",feature = "signature_1Glycolysis_Ucell",use_test = TRUE,comparisons = list(c("CMV", "Non-CMV")),map_signif_level = F)
dev.off()

pdf("plots/paper_figures/oxph_score_ucell.pdf",width=10,height=7)
SCpubr::do_BoxPlot(sample = tcells,group.by = "cmv.anno2",feature = "signature_1Oxidative_Phosphorylation_Ucell",use_test = TRUE,comparisons = list(c("CMV", "Non-CMV")),map_signif_level = F)
dev.off()


