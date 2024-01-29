library(Seurat)
library(infercnv)
library(SeuratDisk)
library(Seurat)
library(rhdf5)
library(anndata)
library(stringr)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(sctransform)
library(ggsci)
library(umap)
library(limma)
library(viridis)
library(patchwork)
library(SingleR)
library(xlsx)
library(harmony)
library(NMF)
library(irGSEA)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(scales)
library(DoubletFinder)

UPS.integrated.tumor <- readRDS("/cluster/home/yflu/UPS/zl.integrated.tumor.b.monocle.rds")
DimPlot(UPS.integrated.tumor,reduction = "umap", label = T, group.by = "celltype")
seuratcluster <- UPS.integrated.tumor@meta.data$celltype
seuratcluster <- as.character(seuratcluster)
seuratcluster = c("CAF" = "CAF", "DKK1+ UPS" = "DKK1 UPS", "CRIP1+ UPS" = "CRIP1 UPS", "NME4+ UPS"= "DLK1 UPS", "MSC" = "MSC")[ as.character(seuratcluster)]
seuratcluster <- factor(seuratcluster , levels = c("CAF", "DKK1 UPS", "CRIP1 UPS", "DLK1 UPS", "MSC"))
UPS.integrated.tumor@meta.data$celltype <- seuratcluster
counts <- GetAssayData(UPS.integrated.tumor, slot = 'counts')
anno <- data.frame(UPS.integrated.tumor@meta.data$celltype)
rownames(anno) <- rownames(UPS.integrated.tumor@meta.data)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts,
                                    annotations_file=anno,
                                    delim="\t",
                                    gene_order_file="/cluster/home/yflu/RT/inferCNV/hg38_gencode_v27.txt",
                                    ref_group_names=c("CAF","MSC")) 
infercnv_obj_1 = infercnv::run(infercnv_obj,
                               cutoff = 0.1,
                               out_dir = "/cluster/home/yflu/UPS/inferCNV/output/", 
                               cluster_by_groups = T,
                               HMM = FALSE,
                               denoise = TRUE,
                               num_threads = 8)
infercnv_obj_2 = infercnv::run(infercnv_obj,
                               cutoff = 0.1,
                               out_dir = "/cluster/home/yflu/UPS/inferCNV/output3/", 
                               cluster_by_groups = T,
                               k_obs_groups = 3,
                               cluster_references = F,
                               HMM = T,
                               analysis_mode = c('samples'),
                               denoise = TRUE,
                               num_threads = 8,
                               inspect_subclusters = F)
infercnv_obj_1 <- readRDS("/cluster/home/yflu/RT/inferCNV/run.final.infercnv_obj")

plot_cnv(infercnv_obj_1,
         out_dir = "/cluster/home/yflu/RT/inferCNV/TEST",
         title = "inferCNV",
         obs_title = "Observations (Cells)",
         ref_title = "References (Cells)",
         cluster_by_groups = T,
         cluster_references = TRUE,
)
expr <- infercnv_obj_2@expr.data
rownames(expr)
match('DLK1',rownames(expr))
rownames(expr)[4947]
expr[4947,]
UPS.integrated.tumor@meta.data$DLK1_CNV <- expr[4947,]
VlnPlot(UPS.integrated.tumor,features = "DLK1_CNV")

infercnv_obj = readRDS("/cluster/home/yflu/UPS/inferCNV/output3/run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
normal_loc <- infercnv_obj@reference_grouped_cell_indices
test_loc <- infercnv_obj@observation_grouped_cell_indices

anno.df=data.frame(
  CB=c(colnames(expr)[normal_loc$CAF],colnames(expr)[normal_loc$MSC]),
  class=c(rep("CAF",length(normal_loc$CAF)),rep("MSC",length(normal_loc$MSC)))
)
anno.df_1=data.frame(
  CB=c(colnames(expr)[test_loc$`CRIP1 UPS`],colnames(expr)[test_loc$`DKK1 UPS`],colnames(expr)[test_loc$`DLK1 UPS`]),
  class=c(rep("CRIP1 UPS",length(test_loc$`CRIP1 UPS`)),rep("DKK1 UPS",length(test_loc$`DKK1 UPS`)),rep("DLK1 UPS",length(test_loc$`DLK1 UPS`)))
)
anno.df <- rbind(anno.df,anno.df_1)

expr2=expr-1
expr2=expr2 ^ 2
CNV_score=as.data.frame(colMeans(expr2))
colnames(CNV_score)="CNV_score"
CNV_score$CB=rownames(CNV_score)
rownames(anno.df) <- anno.df$CB
CNV_score=CNV_score%>%inner_join(anno.df,by="CB")
color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:5] 
names(color_v)=as.character(1:5)
CNV_score%>%ggplot(aes(class,CNV_score))+geom_violin(aes(fill=class))
  theme_bw()
ggsave("CNV level.pdf",width = 10,height = 2,units = "cm")
