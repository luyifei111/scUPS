macrophage.p <- read.xlsx("mmc5.xlsx","Table S6")
colnames(macrophage.p) <- macrophage.p[1,]
macrophage.p <- macrophage.p[-1,]
macrophage.p.M1 <- substr(macrophage.p$M1,1,(nchar(macrophage.p$M1)-1))
#macrophage.p.M1 <- paste(macrophage.p.M1,"+",sep = "")
macrophage.p.M1 <- macrophage.p.M1[1:16]
macrophage.p.M2 <- substr(macrophage.p$M2,1,(nchar(macrophage.p$M2)-1))
#macrophage.p.M2 <- paste(macrophage.p.M2,"+",sep = "")
macrophage.p.Angiogenesis <- substr(macrophage.p$Angiogenesis,1,(nchar(macrophage.p$Angiogenesis)-1))
macrophage.p.Angiogenesis <- macrophage.p.Angiogenesis[1:25]
macrophage.p.Phagocytosis <- substr(macrophage.p$Phagocytosis,1,(nchar(macrophage.p$Phagocytosis)-1))
macrophage.p.Phagocytosis <- macrophage.p.Phagocytosis[1:4]
macro.markers <- list()
macro.markers$M1 <- macrophage.p.M1
macro.markers$M2 <- macrophage.p.M2
macro.markers1 <- list()
macro.markers1$Angiogenesis <- macrophage.p.Angiogenesis
macro.markers1$Phagocytosis <- macrophage.p.Phagocytosis
zl.integrated.M <- irGSEA.score(object = zl.integrated.M, assay = "RNA", slot = "data", 
                             seeds = 123, ncores = 1,msigdb=F, 
                             custom = T, geneset = macro.markers, method = c("AUCell", "UCell", "singscore",
                                                                             "ssgsea"), 
                             kcdf = 'Gaussian')
zl.integrated.M <- irGSEA.score(object = zl.integrated.M, assay = "RNA", slot = "data", 
                                seeds = 123, ncores = 1,msigdb=F, 
                                custom = T, geneset = macro.markers1, method = c("AUCell", "UCell", "singscore",
                                                                                "ssgsea"), 
                                kcdf = 'Gaussian')
scatterplot <- irGSEA.density.scatterplot(object = zl.integrated.M,
                                          method = c("AUCell", "UCell", "singscore",
                                                     "ssgsea")[1],
                                          show.geneset = c("M1","M2"),
                                          reduction = "umap")
scatterplot

scatterplot1 <- irGSEA.density.scatterplot(object = zl.integrated.M,
                                          method = c("AUCell", "UCell", "singscore",
                                                     "ssgsea")[1],
                                          show.geneset = c("Phagocytosis"),
                                          reduction = "umap")
scatterplot1
