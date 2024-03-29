options(stringsAsFactors = F)
setwd('E:/UPS/Integrated/nmf')
library(Seurat)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
F1 <- readRDS("E:/UPS/Integrated/nmf/zl.integrated.tumor.rds")
#directories under nmf are module,program,raw,allmodulegene_value,top30,top50_value
dir.create('E:/UPS/Integrated/nmf/top50_value')

###step1: files preparation#####
samplelist=c("zl1","zl2","zl9","zl10","zl13")
a=paste(samplelist,'.gene_group.csv',sep = '')
dir = a
dir = list.files(pattern = '.gene_group.csv')
a = list.files(pattern = '.gene_group.csv')
dir=paste(getwd(),'/',dir,sep = '')
dir
n=length(dir)
n
for (i in 1:length(dir)) {
  print(dir[i])
  names=substr(a[i],1,nchar(a[i])-15)
  print(names)
  tumorgenemodule=read.csv(dir[i])
  #tumorgenemodule <- t(tumorgenemodule)
  rownames(tumorgenemodule) <- tumorgenemodule[,1]
  tumorgenemodule <- tumorgenemodule[,-1]
  tumorgenename <- rownames(tumorgenemodule)
  tumorgenemodule <- cbind(tumorgenename,tumorgenemodule)
  colnames(tumorgenemodule)[2:11] <-c(1:10) 
  tumorgenemodule <- as.data.frame(tumorgenemodule)
  head(tumorgenemodule)
  tumorgenemodule$max<-apply(tumorgenemodule[,2:11],1,which.max)#改成11
  tumorgenemodule_v1<-tumorgenemodule[which(tumorgenemodule$max==1),][,c(1,2)]
  tumorgenemodule_v1<-tumorgenemodule_v1[order(tumorgenemodule_v1[,2],decreasing = T),][1:50,]
  allgene<-tumorgenemodule_v1
  head(allgene)
  
  #calculate gene top50 value
  for (module in 3:11) {
    aa<-tumorgenemodule[which(tumorgenemodule$max== (module-1)),][,c(1,module)]
    aa<-aa[order(aa[,2],decreasing = T),][1:50,]
    allgene<-cbind(allgene,aa)
  }
  head(allgene)
  
  for (module in 1:10) {
    colnames(allgene)[2*module]=paste(names,'_v',module,sep = '')
    colnames(allgene)[2*module-1]=''
  }
  head(allgene)
  
  write.csv(allgene, file=paste('E:/UPS/Integrated/nmf/top50_value/',names,'_gene_top50_value.csv',sep = ''))
  
  aa<-as.data.frame(allgene[,1])
  for (module in 2:10) {
    mm<-as.data.frame(allgene[,2*module-1])
    aa<-cbind(aa,mm)
  }
  
  for (module in 1:10) {
    colnames(aa)[module]=paste(names,'_v',module,sep = '')
  }
  head(aa)
  #top50gene of each module
  write.csv(aa,file=paste('E:/UPS/Integrated/nmf/top50_value/',names,'_module_top50.csv',sep = ''))
  
  
}

#单独导出#
head(tumorgenemodule)
dir.create('E:/UPS/Integrated/nmf/test/')

if (T){
  for (i in 1:length(dir)) {
    print(i)
    tumorgenemodule=read.csv(dir[i])
    #head(tumorgenemodule)
    tumorgenemodule$max<-apply(tumorgenemodule[,2:11],1,which.max)
    #head(tumorgenemodule)
    names=substr(a[i],1,nchar(a[i])-15)
    print(names)
    for (j in 2:11) {
      aa<-tumorgenemodule[which(tumorgenemodule$max== (j-1)),][,c(1,j)]
      aa<-aa[order(aa[,2],decreasing = T),]
      colnames(aa)[2]=paste(names,j-1,sep = '_v')
      name=paste("E:/UPS/Integrated/nmf/test/",names,'_v',j-1,'_value.csv',sep = '')
      write.csv(aa,name)
    }
  }
}

#####step2:calculate module score######
# make a combined gene file containing all sample modules
a=list.files('E:/UPS/Integrated/nmf/top50_value/',pattern = '_gene_top50_value.csv')
a
dir=paste('E:/UPS/Integrated/nmf/top50_value/',a,sep = '')
dir
n=length(dir)
dir[1]
gene_top50<-read.csv(file = dir[1])
gene_top50<-gene_top50[,-1]
head(gene_top50)
for (i in 2:n) {
  gene_new<-read.csv(file = dir[i])
  gene_new=gene_new[,-1]
  gene_top50<-cbind(gene_top50,gene_new)
}
colnames(gene_top50)
ncol(gene_top50)
head(gene_top50)
View(gene_top50)
write.csv(gene_top50,'allmodulegene_top50.csv')
gene_top50=gene_top50[1:50,]
dim(gene_top50)

##add tumor score
#alltumor=alltumor_harm0907#seurat objecy
colnames(F1@meta.data)
#alltumor@meta.data=alltumor@meta.data[,1:6]
ncol(gene_top50)
View(gene_top50)
genes1 <- gene_top50
#将HLA.A转成HLA-A
if(F){
  genes1=data.frame(1:50)
  for (i in 1:length(rownames(F1@meta.data))) {
    #test=gene_top50[,i]
    test = rownames(F1@meta.data)
    test=gsub('[-]','.',test)#将HLA.A转成HLA-A
    test=as.data.frame(test)
    colnames(test)=colnames(gene_top50)[i]
    genes1=cbind(genes1,test)
  }
}
# genes1 <- genes1[,-1]
# View(genes)
# genes=genes[,-ncol(genes)]
# colnames(genes)
# ncol(genes)
# write.csv(genes,'allmodulegene_top50.csv')

# options(stringsAsFactors = F)
# genes=read.csv('allmodulegene_top30new_16sample.csv')
# genes=genes[,-1]
# 
# 
# colnames(alltumor@meta.data)[1:19]
# alltumor@meta.data=alltumor@meta.data[,1:19]
# head(alltumor@meta.data)
# head(genes)
# head(alltumornew@meta.data)

#genes=gene_top50
#head(gene_top50[,738:740])

for (ii in 1:50){
  i <- 2*ii
  print(i)
  module_names <- colnames(genes1)[i]
  module <- list(genes1[,i-1])
  print(module_names)
  F1<-AddModuleScore(F1,module,ctrl = 100,name= module_names)
  
}

grep('AP',rownames(F1),value = T)
ncol(F1@meta.data)
colnames(F1@meta.data)
head(F1@meta.data)

all<-F1@meta.data
ncol(all)
colnames(all)
View(all)
all_back=all
all=all[,-c(1:11)]
head(all)
ncol(all)
write.csv(all,'alltumor_modulescore_top50score.csv')

# options(stringsAsFactors = F)
# all<-read.csv('alltumor_modulescore.csv')
# rownames(all)=all[,1]
# all=all[,-1]
colnames(all)=substr(colnames(all),1,(nchar(colnames(all))-1))#改名
colnames(all)
ncol(all)

matrix=cor(all)# correlation
write.csv(matrix,'nmf_matrix.csv')

p1 <- pheatmap(matrix,show_rownames=T,show_colnames=T,
               color=rev(colorRampPalette(brewer.pal(11,'RdBu'))(50)),
               breaks = seq(-1,1,0.04),
               clustering_distance_rows = 'canberra',clustering_distance_cols = 'canberra'
)
#clustering_distances:'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'
matrix_cluster <- matrix[p1$tree_row$order, p1$tree_col$order]
colnames(matrix_cluster)
#colnames(matrix_cluster)=substr(colnames(matrix_cluster),1,(nchar(colnames(matrix_cluster))-1))
colname=colnames(matrix_cluster)
colname
match('zl1_v5',colname)
moduleA=colname[1:6]
moduleB=colname[7:11]
moduleC=colname[12:13]
moduleD=colname[19:22]
moduleE=colname[24:26]
moduleF=colname[29:32]
moduleG=colname[33:48]

a=list.files('E:/UPS/Integrated/nmf/top50_value/',pattern = '_gene_top50_value.csv')
a
dir=paste('E:/UPS/Integrated/nmf/top50_value/',a,sep = '')
dir
n=length(dir)
gene_value<-read.csv(file = dir[1])
head(gene_value)
gene_value <- gene_value[,-1]
dir.create('E:/UPS/Integrated/nmf/module')
dir.create('E:/UPS/Integrated/nmf/program')
for (i in 2:n) {
  gene_new<-read.csv(file = dir[i])
  gene_new <- gene_new[,-1]
  gene_value<-cbind(gene_value,gene_new)
}
colnames(gene_value)

top50_value=gene_value[1:50,]
write.csv(top50_value,'genetop50_valuecombine.csv')

data.col=data.frame(c(1:50))

moduleA.data=data.frame(c(1:50))
for (i in moduleA) {
  print(i)
  col=match(i,colnames(top50_value))
  mm=top50_value[,(col-1):col]
  moduleA.data=cbind(moduleA.data,mm)
}
moduleA.data=moduleA.data[,-1]
head(moduleA.data)
write.csv(moduleA.data,'E:/UPS/Integrated/nmf/module/moduleA.csv')

moduleB.data=data.frame(c(1:50))
for (i in moduleB) {
  col=match(i,colnames(top50_value))
  mm=top50_value[,(col-1):col]
  moduleB.data=cbind(moduleB.data,mm)
}
moduleB.data=moduleB.data[,-1]
head(moduleB.data)
write.csv(moduleB.data,'E:/UPS/Integrated/nmf/module/moduleB.csv')

moduleC.data=data.frame(c(1:50))
for (i in moduleC) {
  col=match(i,colnames(top50_value))
  mm=top50_value[,(col-1):col]
  moduleC.data=cbind(moduleC.data,mm)
}
moduleC.data=moduleC.data[,-1]
head(moduleC.data)
write.csv(moduleC.data,'E:/UPS/Integrated/nmf/module/moduleC.csv')

moduleD.data=data.frame(c(1:50))
for (i in moduleD) {
  col=match(i,colnames(top50_value))
  mm=top50_value[,(col-1):col]
  moduleD.data=cbind(moduleD.data,mm)
}
moduleD.data=moduleD.data[,-1]
head(moduleD.data)
write.csv(moduleD.data,'E:/UPS/Integrated/nmf/module/moduleD.csv')

moduleE.data=data.frame(c(1:50))
for (i in moduleE) {
  col=match(i,colnames(top50_value))
  mm=top50_value[,(col-1):col]
  moduleE.data=cbind(moduleE.data,mm)
}
moduleE.data=moduleE.data[,-1]
head(moduleE.data)
write.csv(moduleE.data,'E:/UPS/Integrated/nmf/module/moduleE.csv')

moduleF.data=data.frame(c(1:50))
for (i in moduleF) {
  col=match(i,colnames(top50_value))
  mm=top50_value[,(col-1):col]
  moduleF.data=cbind(moduleF.data,mm)
}
moduleF.data=moduleF.data[,-1]
head(moduleF.data)
write.csv(moduleF.data,'E:/UPS/Integrated/nmf/module/moduleF.csv')

moduleG.data=data.frame(c(1:50))
for (i in moduleG) {
  col=match(i,colnames(top50_value))
  mm=top50_value[,(col-1):col]
  moduleG.data=cbind(moduleG.data,mm)
}
moduleG.data=moduleG.data[,-1]
head(moduleG.data)
write.csv(moduleG.data,'E:/UPS/Integrated/nmf/module/moduleG.csv')

moduleH.data=data.frame(c(1:50))
for (i in moduleH) {
  col=match(i,colnames(top50_value))
  mm=top50_value[,(col-1):col]
  moduleH.data=cbind(moduleH.data,mm)
}
moduleH.data=moduleH.data[,-1]
head(moduleH.data)
write.csv(moduleH.data,'E:/UPS/Integrated/nmf/module/moduleH.csv')

setwd('E:/UPS/Integrated/nmf/module')
a=list.files(pattern = 'module')
a
dir=paste('E:/UPS/Integrated/nmf/module/',a,sep = '')
dir
n=length(dir)
n
m=LETTERS[1:n]###set program names according to module
m
i = dir[1]
for (ii in 1:n) {
  i <- dir[ii]
  print(i)
  ###select modules according to heatmap###
  moduleall<-read.csv(file = i)
  moduleall<-moduleall[,-1]
  head(moduleall)
  m1<-moduleall[1:50,]# select top50 genes
  
  #merge top50 genes of all modules#
  m2<-m1[,1:2]
  colnames(m2)[1]<-'symbol'
  
  for (i in 2:(ncol(moduleall)/2)){
    y=2*i-1
    aa<-m1[,y:(y+1)]
    colnames(aa)[1]<-'symbol'
    m2<-merge(m2,aa,by='symbol',all=T)
  }
  m2[is.na(m2)]<-0
  m2$mean<-apply(m2[,2:ncol(m2)],1,mean)
  m2<-m2[order(m2$mean,decreasing = T),]
  write.csv(m2, file = paste('module',LETTERS[ii],'data.csv', sep = ''))
  
  colname=colnames(m2)[2:ncol(m2)]
  colname
  csvname<-paste("E:/UPS/Integrated/nmf/test/",colname[1],'_value.csv',sep = '')
  csvname
  m3<-read.csv(csvname)[1:50,-1]
  head(m3)
  colnames(m3)[1]<-'symbol'
  m3<-m3[m3$symbol%in%m2$symbol,]
  
  for (i in colname[2:(ncol(m2)-2)]) {
    csvname<-paste("E:/UPS/Integrated/nmf/test/",i,'_value.csv',sep = '')
    aa<-read.csv(csvname)[1:50,-1]
    colnames(aa)[1]<-'symbol'
    aa<-aa[aa$symbol%in%m2$symbol,]
    m3<-merge(m3,aa,by='symbol',all=T)
  }
  
  print(sum(is.na(m2)))
  print(sum(is.na(m3)))
  #make adjustments to m3 and sort m3 according to the values of mean
  m3[is.na(m3)]<-0
  m3$mean<-apply(m3[,2:ncol(m3)],1,mean)
  m3<-m3[order(m3$mean,decreasing = T),]
  
  c<-m[ii]
  print(c)
  name<-paste("E:/UPS/Integrated/nmf/program/program",c,".csv",sep = '')
  write.csv(m3,name)
}  
##############select top 50 expressing gene of each module#################
i =1
namingtag1 <- c('1')
#for (i in 1:300) {
  #namingtag <- namingtag1[i]
  #col = match(namingtag,colnames(gene_top50))
  #mm = gene_top50[,(col-1)]
  #genea <- setdiff(mm, genetag)
  #genetag=c(genetag, genea)
#}
#namingtag1 <- namingtag1[-1]

#moduleA
genetag <- 'gene1'
for (i in 1:length(moduleA)) {
  namingtag <- moduleA[i]
  col = match(namingtag,colnames(moduleA.data))
  mm = moduleA.data[,(col-1)]
  genea <- setdiff(mm, genetag)
  genetag=c(genetag, genea)
}

moduleA_gene <- genetag[-1]
# key method: merge(all=T, by= 'geng_name')
a_data <- data.frame(moduleA_gene)
for (l in 1:length(moduleA)) {
  col1 <- match(moduleA[l],colnames(moduleA.data))
  mmm <- moduleA.data[,col1-1]
  rowname <- setdiff(moduleA_gene, mmm)
  rowname <- cbind(rowname,0)
  colnames(rowname) <- c()
  a_data1 <- moduleA.data[,(col1-1):col1]
  colnames(a_data1)[1] = colnames(a_data)[1]
  colnames(rowname) <- colnames(a_data1)
  a_data1 <- rbind(a_data1,rowname)
  a_data <- merge(a_data, a_data1, by = 'moduleA_gene')
  print(length(rownames(a_data)))
  print(l)
}  
rownames(a_data) <- a_data$moduleA_gene
a_data <- a_data[,-1]

mean = 0.1
for (i in 1:length(moduleA)) {
  mean1 <- a_data[,i]
  mean1 <- as.numeric(mean1)
  mean <- mean + mean1
}
a_data$sum <- mean-0.1
a_data<-a_data[order(a_data$sum,decreasing = T),][1:50,]
write.csv(a_data, file = 'module_A_gene.csv')
#####################
#module B
genetag <- 'gene1'
for (i in 1:length(moduleB)) {
  namingtag <- moduleB[i]
  col = match(namingtag,colnames(moduleB.data))
  mm = moduleB.data[,(col-1)]
  genea <- setdiff(mm, genetag)
  genetag=c(genetag, genea)
}

moduleB_gene <- genetag[-1]
# key method: merge(all=T, by= 'geng_name')
a_data <- data.frame(moduleB_gene)
for (l in 1:length(moduleB)) {
  col1 <- match(moduleB[l],colnames(moduleB.data))
  mmm <- moduleB.data[,col1-1]
  rowname <- setdiff(moduleB_gene, mmm)
  rowname <- cbind(rowname,0)
  colnames(rowname) <- c()
  a_data1 <- moduleB.data[,(col1-1):col1]
  colnames(a_data1)[1] = colnames(a_data)[1]
  colnames(rowname) <- colnames(a_data1)
  a_data1 <- rbind(a_data1,rowname)
  a_data <- merge(a_data, a_data1, all = T, by = 'moduleB_gene')
}  
mean = 0.1
rownames(a_data) <- a_data$moduleB_gene
a_data <- a_data[,-1]

for (i in 1:length(moduleB)) {
  mean1 <- a_data[,i]
  mean1 <- as.numeric(mean1)
  mean <- mean + mean1
}
a_data$sum <- mean-0.1
a_data<-a_data[order(a_data$sum,decreasing = T),][1:50,]

write.csv(a_data, file = 'module_B_gene.csv')
###############################
#module C
genetag <- 'gene1'
for (i in 1:length(moduleC)) {
  namingtag <- moduleC[i]
  col = match(namingtag,colnames(moduleC.data))
  mm = moduleC.data[,(col-1)]
  genea <- setdiff(mm, genetag)
  genetag=c(genetag, genea)
}

moduleC_gene <- genetag[-1]
# key method: merge(all=T, by= 'geng_name')
a_data <- data.frame(moduleC_gene)
for (l in 1:length(moduleC)) {
  col1 <- match(moduleC[l],colnames(moduleC.data))
  mmm <- moduleC.data[,col1-1]
  rowname <- setdiff(moduleC_gene, mmm)
  rowname <- cbind(rowname,0)
  colnames(rowname) <- c()
  a_data1 <- moduleC.data[,(col1-1):col1]
  colnames(a_data1)[1] = colnames(a_data)[1]
  colnames(rowname) <- colnames(a_data1)
  a_data1 <- rbind(a_data1,rowname)
  a_data <- merge(a_data, a_data1, all = T, by = 'moduleC_gene')
}  
head(a_data)
rownames(a_data) <- a_data$moduleC_gene
a_data <- a_data[,-1]
a_data_back <- a_data

mean = 0.1
for (i in 1:length(moduleC)) {
  mean1 <- a_data[,i]
  mean1 <- as.numeric(mean1)
  mean <- mean + mean1
}
a_data$sum <- mean-0.1
a_data<-a_data[order(a_data$sum,decreasing = T),][1:50,]

write.csv(a_data, file = 'module_C_gene.csv')
################################
#module D
genetag <- 'gene1'
for (i in 1:length(moduleD)) {
  namingtag <- moduleD[i]
  col = match(namingtag,colnames(moduleD.data))
  mm = moduleD.data[,(col-1)]
  genea <- setdiff(mm, genetag)
  genetag=c(genetag, genea)
}

moduleD_gene <- genetag[-1]
# key method: merge(all=T, by= 'geng_name')
a_data <- data.frame(moduleD_gene)
for (l in 1:length(moduleD)) {
  col1 <- match(moduleD[l],colnames(moduleD.data))
  mmm <- moduleD.data[,col1-1]
  rowname <- setdiff(moduleD_gene, mmm)
  rowname <- cbind(rowname,0)
  colnames(rowname) <- c()
  a_data1 <- moduleD.data[,(col1-1):col1]
  colnames(a_data1)[1] = colnames(a_data)[1]
  colnames(rowname) <- colnames(a_data1)
  a_data1 <- rbind(a_data1,rowname)
  a_data <- merge(a_data, a_data1, all = T, by = 'moduleD_gene')
}  
head(a_data)
rownames(a_data) <- a_data$moduleD_gene
a_data <- a_data[,-1]
a_data_back <- a_data

mean = 0.1
for (i in 1:length(moduleD)) {
  mean1 <- a_data[,i]
  mean1 <- as.numeric(mean1)
  mean <- mean + mean1
}
a_data$sum <- mean-0.1
a_data<-a_data[order(a_data$sum,decreasing = T),][1:50,]

write.csv(a_data, file = 'module_D_gene.csv')
#################################
#module E
genetag <- 'gene1'
for (i in 1:length(moduleE)) {
  namingtag <- moduleE[i]
  col = match(namingtag,colnames(moduleE.data))
  mm = moduleE.data[,(col-1)]
  genea <- setdiff(mm, genetag)
  genetag=c(genetag, genea)
}

moduleE_gene <- genetag[-1]
# key method: merge(all=T, by= 'geng_name')
a_data <- data.frame(moduleE_gene)
for (l in 1:length(moduleE)) {
  col1 <- match(moduleE[l],colnames(moduleE.data))
  mmm <- moduleE.data[,col1-1]
  rowname <- setdiff(moduleE_gene, mmm)
  rowname <- cbind(rowname,0)
  colnames(rowname) <- c()
  a_data1 <- moduleE.data[,(col1-1):col1]
  colnames(a_data1)[1] = colnames(a_data)[1]
  colnames(rowname) <- colnames(a_data1)
  a_data1 <- rbind(a_data1,rowname)
  a_data <- merge(a_data, a_data1, all = T, by = 'moduleE_gene')
}  
mean = 0.1
rownames(a_data) <- a_data$moduleE_gene
a_data <- a_data[,-1]
a_data_back <- a_data

for (i in 1:length(moduleE)) {
  mean1 <- a_data[,i]
  mean1 <- as.numeric(mean1)
  mean <- mean + mean1
}
a_data$sum <- mean-0.1
a_data<-a_data[order(a_data$sum,decreasing = T),][1:50,]

write.csv(a_data, file = 'module_E_gene.csv')
##################################
#module F
genetag <- 'gene1'
for (i in 1:length(moduleF)) {
  namingtag <- moduleF[i]
  col = match(namingtag,colnames(moduleF.data))
  mm = moduleF.data[,(col-1)]
  genea <- setdiff(mm, genetag)
  genetag=c(genetag, genea)
}

moduleF_gene <- genetag[-1]
# key method: merge(all=T, by= 'geng_name')
a_data <- data.frame(moduleF_gene)
for (l in 1:length(moduleF)) {
  col1 <- match(moduleF[l],colnames(moduleF.data))
  mmm <- moduleF.data[,col1-1]
  rowname <- setdiff(moduleF_gene, mmm)
  rowname <- cbind(rowname,0)
  colnames(rowname) <- c()
  a_data1 <- moduleF.data[,(col1-1):col1]
  colnames(a_data1)[1] = colnames(a_data)[1]
  colnames(rowname) <- colnames(a_data1)
  a_data1 <- rbind(a_data1,rowname)
  a_data <- merge(a_data, a_data1, all = T, by = 'moduleF_gene')
}  
mean = 0.1
rownames(a_data) <- a_data$moduleF_gene
a_data <- a_data[,-1]
a_data_back <- a_data
for (i in 1:length(moduleF)) {
  mean1 <- a_data[,i]
  mean1 <- as.numeric(mean1)
  mean <- mean + mean1
}
a_data$sum <- mean-0.1
a_data<-a_data[order(a_data$sum,decreasing = T),][1:50,]

write.csv(a_data, file = 'module_F_gene.csv')

##################################
#module G
genetag <- 'gene1'
for (i in 1:length(moduleG)) {
  namingtag <- moduleG[i]
  col = match(namingtag,colnames(moduleG.data))
  mm = moduleG.data[,(col-1)]
  genea <- setdiff(mm, genetag)
  genetag=c(genetag, genea)
}

moduleG_gene <- genetag[-1]
# key method: merge(all=T, by= 'geng_name')
a_data <- data.frame(moduleG_gene)
for (l in 1:length(moduleG)) {
  col1 <- match(moduleG[l],colnames(moduleG.data))
  mmm <- moduleG.data[,col1-1]
  rowname <- setdiff(moduleG_gene, mmm)
  rowname <- cbind(rowname,0)
  colnames(rowname) <- c()
  a_data1 <- moduleG.data[,(col1-1):col1]
  colnames(a_data1)[1] = colnames(a_data)[1]
  colnames(rowname) <- colnames(a_data1)
  a_data1 <- rbind(a_data1,rowname)
  a_data <- merge(a_data, a_data1, all = T, by = 'moduleG_gene')
}  
mean = 0.1
rownames(a_data) <- a_data$moduleG_gene
a_data <- a_data[,-1]
a_data_back <- a_data
for (i in 1:length(moduleG)) {
  mean1 <- a_data[,i]
  mean1 <- as.numeric(mean1)
  mean <- mean + mean1
}
a_data$sum <- mean-0.1
a_data<-a_data[order(a_data$sum,decreasing = T),][1:50,]

write.csv(a_data, file = 'module_G_gene.csv')

##################################
#module H
genetag <- 'gene1'
for (i in 1:length(moduleH)) {
  namingtag <- moduleH[i]
  col = match(namingtag,colnames(moduleH.data))
  mm = moduleH.data[,(col-1)]
  genea <- setdiff(mm, genetag)
  genetag=c(genetag, genea)
}

moduleH_gene <- genetag[-1]
# key method: merge(all=T, by= 'geng_name')
a_data <- data.frame(moduleH_gene)
for (l in 1:length(moduleH)) {
  col1 <- match(moduleH[l],colnames(moduleH.data))
  mmm <- moduleH.data[,col1-1]
  rowname <- setdiff(moduleH_gene, mmm)
  rowname <- cbind(rowname,0)
  colnames(rowname) <- c()
  a_data1 <- moduleH.data[,(col1-1):col1]
  colnames(a_data1)[1] = colnames(a_data)[1]
  colnames(rowname) <- colnames(a_data1)
  a_data1 <- rbind(a_data1,rowname)
  a_data <- merge(a_data, a_data1, all = T, by = 'moduleH_gene')
}  
mean = 0.1
rownames(a_data) <- a_data$moduleH_gene
a_data <- a_data[,-1]
a_data_back <- a_data
for (i in 1:length(moduleH)) {
  mean1 <- a_data[,i]
  mean1 <- as.numeric(mean1)
  mean <- mean + mean1
}
a_data$sum <- mean-0.1
a_data<-a_data[order(a_data$sum,decreasing = T),][1:50,]

write.csv(a_data, file = 'module_H_gene.csv')

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
programA <- read.csv(file = 'moduleAdata.csv')
programB <- read.csv(file = 'moduleBdata.csv')
programC <- read.csv(file = 'moduleCdata.csv')
programD <- read.csv(file = 'moduleDdata.csv')
programE <- read.csv(file = 'moduleEdata.csv')
programF <- read.csv(file = 'moduleFdata.csv')
programG <- read.csv(file = 'moduleGdata.csv')
programH <- read.csv(file = 'moduleHdata.csv')
programA <- programA[,-1]
programB <- programB[,-1]
programC <- programC[,-1]
programD <- programD[,-1]
programE <- programE[,-1]
programF <- programF[,-1]
programG <- programG[,-1]
programH <- programH[,-1]
aa=programA$symbol[1:50]
bb=programB$symbol[1:50]
cc=programC$symbol[1:50]
dd=programD$symbol[1:50]
ee=programE$symbol[1:50]
ff=programF$symbol[1:50]
gg=programG$symbol[1:50]
hh=programH$symbol[1:50]

genelist_aa <- bitr(aa, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
genelist_bb <- bitr(bb, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
genelist_cc <- bitr(cc, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
genelist_dd <- bitr(dd, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
genelist_ee <- bitr(ee, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
genelist_ff <- bitr(ff, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
genelist_gg <- bitr(gg, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
genelist_hh <- bitr(hh, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")

head(genelist_aa)
str(genelist_aa)
go_aa <- enrichGO(genelist_aa$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                  qvalueCutoff = 0.01,keyType = 'ENTREZID',readable = T)
go_bb <- enrichGO(genelist_bb$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                  qvalueCutoff = 0.01,keyType = 'ENTREZID',readable = T)
go_cc <- enrichGO(genelist_cc$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                  qvalueCutoff = 0.01,keyType = 'ENTREZID',readable = T)
go_dd <- enrichGO(genelist_dd$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                  qvalueCutoff = 0.01,keyType = 'ENTREZID',readable = T)
go_ee <- enrichGO(genelist_ee$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                  qvalueCutoff = 0.01,keyType = 'ENTREZID',readable = T)
go_ff <- enrichGO(genelist_ff$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                  qvalueCutoff = 0.01,keyType = 'ENTREZID',readable = T)
go_gg <- enrichGO(genelist_gg$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                  qvalueCutoff = 0.01,keyType = 'ENTREZID',readable = T)
go_hh <- enrichGO(genelist_hh$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05,
                  qvalueCutoff = 0.01,keyType = 'ENTREZID',readable = T)
write.csv(genelist_aa,file = 'genelist_aa_go.csv')
write.csv(genelist_bb,file = 'genelist_bb_go.csv')
write.csv(genelist_cc,file = 'genelist_cc_go.csv')
write.csv(genelist_dd,file = 'genelist_dd_go.csv')
write.csv(genelist_ee,file = 'genelist_ee_go.csv')
write.csv(genelist_ff,file = 'genelist_ff_go.csv')
write.csv(genelist_gg,file = 'genelist_gg_go.csv')
write.csv(genelist_hh,file = 'genelist_hh_go.csv')
write.csv(go_aa@result,file = 'program_aa_go.csv')
write.csv(go_bb@result,file = 'program_bb_go.csv')
write.csv(go_cc@result,file = 'program_cc_go.csv')
write.csv(go_dd@result,file = 'program_dd_go.csv')
write.csv(go_ee@result,file = 'program_ee_go.csv')
write.csv(go_ff@result,file = 'program_ff_go.csv')
write.csv(go_gg@result,file = 'program_gg_go.csv')
write.csv(go_hh@result,file = 'program_hh_go.csv')

#programA
barplot(go_aa,showCategory =10,drop=T)
back=go_aa
result<-go_aa@result
result<-result[order(result$p.adjust),]
back@result=result
barplot(back,showCategory=20,drop=T)

#programB
barplot(go_bb,showCategory =10,drop=T)
back=go_bb
result<-go_bb@result
result<-result[order(result$p.adjust),]
back@result=result
barplot(back,showCategory=10,drop=T)


#programC
barplot(go_cc,showCategory =10,drop=T)
back=go_cc
result<-go_cc@result
result<-result[order(result$Count,decreasing = T),]
back@result=result
barplot(back,showCategory=10,drop=T)

#programD
barplot(go_dd,showCategory =10,drop=T)
back=go_dd
result<-go_dd@result
result<-result[order(result$p.adjust),]
back@result=result
barplot(back,showCategory=10,drop=T)


#programE
barplot(go_ee,showCategory =10,drop=T)
back=go_ee
result<-go_ee@result
result<-result[order(result$p.adjust),]
back@result=result
barplot(back,showCategory=10,drop=T)

#programF
barplot(go_ff,showCategory =10,drop=T)
back=go_ff
result<-go_ff@result
result<-result[order(result$p.adjust),]
back@result=result
barplot(back,showCategory=10,drop=T)

#programG
barplot(go_gg,showCategory =10,drop=T)
back=go_gg
result<-go_gg@result
result<-result[order(result$p.adjust),]
back@result=result
barplot(back,showCategory=10,drop=T)

#programH
barplot(go_hh,showCategory =10,drop=T)
back=go_hh
result<-go_hh@result
result<-result[order(result$p.adjust),]
back@result=result
barplot(back,showCategory=10,drop=T)

kegg <- enrichKEGG(genelist_ff$ENTREZID, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                   minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,use_internal_data = FALSE)



kegg_result<-kegg@result
kegg_result<-kegg_result[order(kegg_result$p.adjust),]
kegg@result=kegg_result
barplot(kegg,showCategory = 10)
write.csv(kegg@result,file = 'program_D_kegg.csv')


setwd('E:/UPS/Integrated/nmf')
library(pheatmap)
library(RColorBrewer)
programA<-read.csv('E:/UPS/Integrated/nmf/program/programA.csv')
programB<-read.csv('E:/UPS/Integrated/nmf/program/programB.csv')
programC<-read.csv('E:/UPS/Integrated/nmf/program/programC.csv')
programD<-read.csv('E:/UPS/Integrated/nmf/program/programD.csv')
programE<-read.csv('E:/UPS/Integrated/nmf/program/programE.csv')
programF<-read.csv('E:/UPS/Integrated/nmf/program/programF.csv')
programG<-read.csv('E:/UPS/Integrated/nmf/program/programG.csv')
programH<-read.csv('E:/UPS/Integrated/nmf/program/programH.csv')
a<-data.frame(symbol=programA$symbol[1:50],program=c(rep('programA',50)))
b<-data.frame(symbol=programB$symbol[1:50],program=c(rep('programB',50)))
c<-data.frame(symbol=programC$symbol[1:50],program=c(rep('programC',50)))
d<-data.frame(symbol=programD$symbol[1:50],program=c(rep('programD',50)))
e<-data.frame(symbol=programE$symbol[1:50],program=c(rep('programE',50)))
f<-data.frame(symbol=programF$symbol[1:50],program=c(rep('programF',50)))
g<-data.frame(symbol=programG$symbol[1:50],program=c(rep('programG',50)))
h<-data.frame(symbol=programH$symbol[1:50],program=c(rep('programH',50)))
allgene<-rbind(a,b,c,d,e,f,g)
allgene_back=allgene
write.csv(allgene,'allprogramgene_top50.csv')
seuratcluster <- F1@meta.data$seurat_clusters
seuratcluster <- as.character(seuratcluster)
seuratcluster = c("0" = "CAF", "1" = "DKK1+ UPS", "2" = "DKK1 UPS", "3"= "CRIP1 UPS", 
                  "4" = "DKK1 UPS", "5" = "DLK1 UPS","6" = "DLK1 UPS",
                  "7" = "DLK1 UPS", "8" = "CRIP1 UPS", "9" = "CAF", 
                  "10" = "MSC", "11" = "DKK1+ UPS", "12" = "DLK1 UPS",
                  "13" = "DLK1 UPS","14" = "CRIP1 UPS", "15" = "CRIP1 UPS",
                  "16" = "CRIP1 UPS","17" = "DKK1 UPS", "18" = "MSC"
                  )[ as.character(seuratcluster)]
seuratcluster <- factor(seuratcluster , levels = c("DKK1+ UPS", "CRIP1+ UPS", "DLK1+ UPS"))
F1@meta.data$cell_type <- seuratcluster
alltumor.data<-as.data.frame(GetAssayData(F1))
ac=data.frame(cluster=F1@meta.data$cell_type,sample=F1@meta.data$orig.ident)
rownames(ac) <- colnames(alltumor.data)
head(ac)
head(allgene)
allgene=allgene[!duplicated(allgene$symbol),]
ar=data.frame(row.names=allgene[,1],program=allgene[,2])
head(ar)
plot<-alltumor.data[allgene[,1],]
plot[1:5,1:5]
max(plot)
min(plot)
bk=seq(-1.5,1.5,0.06)
p=pheatmap(plot,show_rownames=T,
           scale = 'row',
           clustering_method = 'ward.D2',
           show_colnames=F,cluster_cols = TRUE,cluster_rows = FALSE,
           color=rev(colorRampPalette(brewer.pal(11,'RdBu'))(50)),
           annotation_col = ac,annotation_row=ar,
           breaks = bk,
           #clustering_distance_rows = 'binary'
)

#分群尝试
head(F1@meta.data)
colnames(F1@meta.data)
alltumor <- F1
alltumor@meta.data=alltumor@meta.data[,c(1:9,15,17,18,19)]

head(allgene)
allgene=split(allgene$symbol,allgene$program)
head(allgene)
for (j in 1:length(allgene)) {
  print(j)
  F1=AddModuleScore(F1,features = list(allgene[[j]]),name = names(allgene)[[j]])
}

#for (j in 1:length(allgene)) {
#  print(j)
#  alltumornew=AddModuleScore(alltumornew,features = list(allgene[[j]]),name = names(allgene)[[j]])
#}

Idents(F1)=F1$seurat_clusters
?FeaturePlot
FeaturePlot(F1,features = c(colnames(F1@meta.data)[63:69]),reduction = 'umap',
            order  = T,label = F,pt.size = 0.6,
            #min.cutoff = 'q10',
            cols = rev(colorRampPalette(c('red4','orangered3','orangered','#FFFF4D','#4169E1','#2A52BE','#000080'))(50)))
ggsave('E:/UPS/Integrated/nmf/figs/programB.png',width = 7.5,height = 6)

FeaturePlot(alltumornew,features = c(colnames(alltumornew@meta.data)[13]),reduction = 'umap',
            sort.cell = T,label = F,pt.size = 0.4,
            #min.cutoff = 'q5',
            cols = rev(colorRampPalette(c('red4','orangered3','orangered','#FFFF4D','#4169E1','#2A52BE','#000080'))(50)))
ggsave('../paper_figures/programE.tiff',width = 7,height = 6)
saveRDS(alltumornew,'../rds/alltumor_harm0902.rds')

DimPlot(F1,label = T)
head(alltumornew@meta.data)
colnames(alltumor@meta.data)
aa=unique(alltumor$sample)
aa
table(Idents(alltumor))

Idents(alltumor)=factor(alltumor$sample,levels = c(aa[c(1,3:4,7:13,6,5,2,16,15,14)]))
for (i in 62:69) {
  p=VlnPlot(F1,features = c(colnames(F1@meta.data)[i]),group.by="seurat_clusters")
  print(p)
}

#求每个群每个模块平均�?
alltumor.data<-as.data.frame(GetAssayData(F1))
Idents(F1)=F1@meta.data$seurat_clusters
allgene=allgene_back
celltype=levels(Idents(F1))
celltype
cluster<-subset(F1,idents = celltype[1])
cluster.data<-as.data.frame(GetAssayData(cluster))[allgene$symbol,]
cluster.data[1:5,1:5]
cluster.data_mean<-apply(cluster.data,1,mean)
all<-as.data.frame(cluster.data_mean)
for (i in 2:length(celltype)) {
  cluster<-subset(F1,idents = celltype[i])
  cluster.data<-as.data.frame(GetAssayData(cluster))[allgene$symbol,]
  cluster.data[1:5,1:5]
  cluster.data_mean<-apply(cluster.data,1,mean)
  all<-cbind(all,cluster.data_mean)
}
colnames(all)<-paste('cluster',celltype,sep = '')
rownames(all)<-rownames(cluster.data)
table(allgene$program)

plot=as.data.frame(matrix(0:14,nrow = 1))
for (i in unique(allgene$program)) {
  filt=allgene[which(allgene$program==i),]$symbol
  mm=all[filt,]
  mm_mean=apply(mm, 2, mean)
  plot=rbind(plot,mm_mean)
}
plot=plot[-1,]
rownames(plot)=unique(allgene$program)
colnames(plot)=colnames(all)

deg=FindAllMarkers(F1)


##remove module C

head(all)
col=match(moduleC,colnames(all))
all=all[,-col]