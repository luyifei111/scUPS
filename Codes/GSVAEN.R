library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(fgsea)
library(presto)


tumornames <- Idents(zl.integrated.tumor)
tumornames <- unique(tumornames)
expr <- AverageExpression(zl.integrated.tumor, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)

genesets.hall <- msigdbr(species = "Homo sapiens", category = "H") 
genesets.hall <- subset(genesets.hall, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets.hall <- split(genesets.hall$gene_symbol, genesets.hall$gs_name)


gsva.res.hall <- gsva(expr, genesets.hall, method="ssgsea") 
saveRDS(gsva.res.hall, "gsva.res.hall.rds")
gsva.df.hall <- data.frame(Genesets=rownames(gsva.res.hall), gsva.res.hall, check.names = F)
write.csv(gsva.df.hall, "gsva_res.hall.csv", row.names = F)
pheatmap::pheatmap(gsva.res.hall, show_colnames = T, scale = "row")

genesets.kegg <- msigdbr(species = "Homo sapiens", category = "C2") 
genesets.kegg <- subset(genesets.kegg, gs_subcat=="CP:KEGG", select = c("gs_name", "gene_symbol")) %>% as.data.frame()
genesets.kegg <- split(genesets.kegg$gene_symbol, genesets.kegg$gs_name)

gsva.res.kegg <- gsva(expr, genesets.kegg, method="ssgsea") 
saveRDS(gsva.res.kegg, "gsva.res.kegg.rds")
gsva.df.kegg <- data.frame(Genesets=rownames(gsva.res.kegg), gsva.res.kegg, check.names = F)
write.csv(gsva.df.kegg, "gsva_res.kegg.csv", row.names = F)
pheatmap::pheatmap(gsva.res.kegg, show_colnames = T, scale = "row")

genesets.go <- msigdbr(species = "Homo sapiens", category = "C5") 
genesets.go <- subset(genesets.go, gs_subcat=="GO:BP", select = c("gs_name", "gene_symbol")) %>% as.data.frame()
genesets.go <- split(genesets.go$gene_symbol, genesets.go$gs_name)

gsva.res.go <- gsva(expr, genesets.go, method="ssgsea") 
saveRDS(gsva.res.go, "gsva.res.go.rds")
gsva.df.go <- data.frame(Genesets=rownames(gsva.res.go), gsva.res.go, check.names = F)
write.csv(gsva.df.go, "gsva_res.go.csv", row.names = F)
p.go <- pheatmap::pheatmap(gsva.res.go, show_colnames = T, scale = "row")

gsva.df.go3 <- gsva.df.go
gsva.df.go3 <- gsva.df.go[order(gsva.df.go[,2], decreasing = T),]
gsva.df.go1 <- gsva.df.go3[1:10,]
gsva.df.go3 <- gsva.df.go3[-c(1:10),]

for (i in c(3:6)) {
  gsva.df.go3 <- gsva.df.go3[order(gsva.df.go3[,i], decreasing = T),]
  gsva.df.go2 <- gsva.df.go3[1:10,]
  gsva.df.go1 <- rbind(gsva.df.go1,gsva.df.go2)
  gsva.df.go3 <- gsva.df.go3[-c(1:10),]
}

p.go <- pheatmap::pheatmap(gsva.df.go1[,-1], show_colnames = T, scale = "row")

expr.sc <- as.matrix(zl.integrated.tumor@assays$RNA@data)
expr.sc <- expr.sc[rowSums(expr.sc)>0,]  #选取非零基因

gsva.res.kegg.sc <- gsva(expr.sc, genesets.kegg, method="ssgsea") 
saveRDS(gsva.res.kegg.sc, "gsva.res.kegg.sc.rds")
gsva.df.kegg.sc <- data.frame(Genesets=rownames(gsva.res.kegg.sc), gsva.res.kegg.sc, check.names = F)
write.csv(gsva.df.kegg.sc, "gsva_res.kegg.sc.csv", row.names = F)
pheatmap::pheatmap(gsva.res.kegg.sc, show_colnames = F, scale = "row")

zl.integrated.tumor@meta.data$group <- recode(Idents(zl.integrated.tumor),
                                              "DKK1+ UPS" = "DKK1 Trajectory",
                                              "CRIP1+ UPS" = "DKK1 Trajectory",
                                              "MSC" = "MSC",
                                              "CAF" = "DLK1 Trajectory",
                                              "DLK1+ UPS" = "DLK1 Trajectory",)
zl.integrated.tumor@meta.data$cell_type <- Idents(zl.integrated.tumor)
zl.integrated.tumor.1 <- subset(zl.integrated.tumor,idents = c("DKK1+ UPS","CRIP1+ UPS","CAF","DLK1+ UPS"))

deg <- FindMarkers(zl.integrated.tumor.1, ident.1 = "DKK1+ UPS",ident.2 = "NME4+ UPS", assay = "RNA")
deg.1 <- FindMarkers(zl.integrated.tumor.1, ident.1 = "NME4+ UPS",ident.2 = "DKK1+ UPS", assay = "RNA")
deg.2 <- FindMarkers(zl.integrated.tumor.1, ident.1 = "CRIP1+ UPS",ident.2 = "CAF", assay = "RNA")
deg.3 <- FindMarkers(zl.integrated.tumor.1, ident.1 = "CAF",ident.2 = "CRIP1+ UPS", assay = "RNA")
write.csv(deg, file = "DEG.DKK1.csv", row.names = F)
write.csv(deg.1, file = "DEG.DLK1.csv", row.names = F)
write.csv(deg.2, file = "DEG.CRIP1.csv", row.names = F)
write.csv(deg.3, file = "DEG.CAF.csv", row.names = F)

auc <- wilcoxauc(zl.integrated.tumor.1)

genesets.kegg1 = msigdbr(species = "Homo sapiens", category = "C2") 
genesets.kegg1 <- subset(genesets.kegg1, gs_subcat=="CP:KEGG", select = c("gs_name", "gene_symbol"))

deg <- deg[order(deg$avg_log2FC, decreasing = T),]
genelist <- structure(deg$avg_log2FC, names = rownames(deg))
res <- GSEA(genelist, TERM2GENE = genesets.kegg1,eps = 0)


write.csv(res, file = "kegg.csv", row.names = F)

genesets.go1 = msigdbr(species = "Homo sapiens", category = "C5")
genesets.go2 <- genesets.go1 %>% split(x = .$gene_symbol, f = .$gs_name)
genesets.go1 <- subset(genesets.go1, gs_subcat=="GO:BP", select = c("gs_name", "gene_symbol"))

deg <- deg[order(deg$avg_log2FC, decreasing = T),]
genelist <- structure(deg$avg_log2FC, names = rownames(deg))
res1 <- GSEA(genelist, TERM2GENE = genesets.go1,eps = 0)

write.csv(res1, file = "go.DKK1.csv", row.names = F)

deg.1 <- deg.1[order(deg.1$avg_log2FC, decreasing = T),]
genelist2 <- structure(deg.1$avg_log2FC, names = rownames(deg.1))
res2 <- GSEA(genelist2, TERM2GENE = genesets.go1,eps = 0)

write.csv(res2, file = "go.DLK1.csv", row.names = F)

deg.2 <- deg.2[order(deg.2$avg_log2FC, decreasing = T),]
genelist3 <- structure(deg.2$avg_log2FC, names = rownames(deg.2))
res3 <- GSEA(genelist3, TERM2GENE = genesets.go1,eps = 0)

write.csv(res3, file = "go.CRIP1.csv", row.names = F)

deg.3 <- deg.3[order(deg.3$avg_log2FC, decreasing = T),]
genelist4 <- structure(deg.3$avg_log2FC, names = rownames(deg.3))
res4 <- GSEA(genelist4, TERM2GENE = genesets.go1,eps = 0)

write.csv(res4, file = "go.CAF.csv", row.names = F)


cluster.genes<- auc%>%
  dplyr::filter(group == "DLK1+ UPS") %>%
  arrange(desc(auc)) %>%
  dplyr::select(feature, auc)

ranks<- deframe(cluster.genes)

plotEnrichment(genesets.go2[["GOBP_RIBOSOME_ASSEMBLY"]],
               ranks) + labs(title="GOBP_RIBOSOME_ASSEMBLY")

fgseaRes <- fgsea(pathways = genesets.go2["GOBP_RIBOSOME_ASSEMBLY"], 
                  stats    = ranks)

deg1 <- FindMarkers(zl.integrated.tumor.1, ident.1 = "DLK1 Trajectory", assay = "RNA",
                   group.by = "group",only.pos = T)

deg1 <- deg1[order(deg1$avg_log2FC, decreasing = T),]
genelist1 <- structure(deg1$avg_log2FC, names = rownames(deg1))
res_1 <- GSEA(genelist1, TERM2GENE = genesets.kegg1,eps = 0)


write.csv(res_1, file = "kegg1.csv", row.names = F)

deg1 <- deg1[order(deg1$avg_log2FC, decreasing = T),]
genelist1 <- structure(deg1$avg_log2FC, names = rownames(deg1))
res1_1 <- GSEA(genelist1, TERM2GENE = genesets.go1,eps = 0)

write.csv(res1_1, file = "go1_1.csv", row.names = F)

cluster.genes <- auc%>%
  dplyr::filter(group == "DLK1 Trajectory") %>%
  arrange(desc(auc)) %>%
  dplyr::select(feature, auc)

ranks<- deframe(cluster.genes)

plotEnrichment(genesets.go2[["GOBP_POSITIVE_REGULATION_OF_MAPK_CASCADE"]],
               ranks) + labs(title="GOBP_POSITIVE_REGULATION_OF_MAPK_CASCADE")
openxlsx::write.xlsx(fgseaRes_1, "GSVAres.xlsx")
