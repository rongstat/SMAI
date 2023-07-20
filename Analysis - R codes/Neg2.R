#######################################################################
###############  real data analysis: Task Neg2
###############     May 28, 2023
#######################################################################

library(HDF5Array)
library(zellkonverter)
library(uwot)
library(ggrepel)
library(Seurat)
library(ggplot2)
library(cluster)
library(uwot)
library(fossil)
library(batchelor)
library(harmony)
library(rliger)
library(reticulate)
library(SMAI)
########################################
######### load data and preprocessing
########################################

data0=readH5AD("Lung_atlas_public.h5ad",reader = "python")

counts = data0@assays@data@listData$counts
dim(counts)
colnames(counts) = rownames(data0@colData)
rownames(counts) = names(data0@rowRanges)

meta.data= data.frame(method = data0@colData$batch, cell_type=data0@colData$cell_type)
row.names(meta.data)=c(colnames(counts))
all.data <- CreateSeuratObject(counts = counts, project = "lung", min.cells = 3,
                               min.features = 10, meta.data = meta.data)
ct.X1 = meta.data$cell_type[which(meta.data$method=="1")]
ct.X2 = meta.data$cell_type[which(meta.data$method=="A3")]
table(ct.X1)
table(ct.X2)
ct.X1=factor(ct.X1)
ct.X2=factor(ct.X2)

data.X1 = as.matrix(counts[,which(meta.data$method=="1")])
data.X2 =  as.matrix(counts[,which(meta.data$method=="A3")])


################################################
########### alignment with different methods
################################################

#preprocessing and normalization
all.RNA = cbind(data.X1, data.X2)
dim(all.RNA)
meta.data=data.frame(cell_type = factor(c(ct.X1,ct.X2)),
                     method = factor(c(rep("data1",length(ct.X1)), rep("data2",length(ct.X2)))))
row.names(meta.data)=c(colnames(all.RNA))
all.data <- CreateSeuratObject(counts = all.RNA, project = "panc8", min.cells = 3,
                               min.features = 10, meta.data = meta.data)
all.data <- SplitObject(all.data, split.by = "method")
all.data <- lapply(X = all.data, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- NormalizeData(x)
})
features <- SelectIntegrationFeatures(object.list = all.data)
data.anchors <- FindIntegrationAnchors(object.list = all.data, anchor.features = features)
data.combined <- IntegrateData(anchorset = data.anchors)
dim(data.combined@assays$integrated@data)
data1=as.matrix(data.combined@assays$integrated@data[,1:length(ct.X1)])
data2=as.matrix(data.combined@assays$integrated@data[,(length(ct.X1)+1):(length(ct.X1)+length(ct.X2))])
data.X1.smai = as.matrix(all.data$data1@assays$RNA@data)
data.X2.smai = as.matrix(all.data$data2@assays$RNA@data)
data.X1.smai = data.X1.smai[match(features,rownames(data.X1.smai)),]
data.X2.smai = data.X2.smai[match(features,rownames(data.X2.smai)),]



############## SMAI
#identify maximal correspondence set
mnn.out = findMutualNN(t(data1),t(data2), k1=50)
prop.align.par = min(c(length(unique(mnn.out$first))/dim(data1)[2],
                       length(unique(mnn.out$second))/dim(data2)[2]))
prop.align.par

#start SMAI
set.seed(33)
out.smai <- align(data.X1.smai, data.X2.smai, sel1=unique(mnn.out$first),  sel2=unique(mnn.out$second),
                  t=5, knn=30, r.max=200, dir.map = "auto", denoise="scree",
                  prop.align=prop.align.par,  cutoff = 1.001, cutoff.t=1.5)
out.smai$p.value
#[1] 0.002272086

meta_data=data.frame(cell_type = factor(c(ct.X1,ct.X2)),
                     batch = factor(c(rep("X1",length(ct.X1)), rep("X2",length(ct.X2)))))
data.int = t(cbind(out.smai$data1.integrate,out.smai$data2.integrate))

out=matrix(ncol=6,nrow=6)
set.seed(10)
out0=assess(data.int, meta_data)
out[1,] = unlist(out0)

f.out = matrix(ncol=6,nrow=6)
f.out[1,1:3] = faithful(data.int[1:length(ct.X1),], t(data.X1.smai))
f.out[1,4:6] = faithful(data.int[-c(1:length(ct.X1)),], t(data.X2.smai))


############Seurat
all.RNA = cbind(data.X1,data.X2)
dim(all.RNA)
meta.data= data.frame(method = c(rep("data1", dim(data.X1)[2]), rep("data2", dim(data.X2)[2])))
row.names(meta.data)=c(colnames(data.X1), colnames(data.X2))
all.data <- CreateSeuratObject(counts = all.RNA, project = "align", min.cells = 3, min.features = 10,
                               meta.data = meta.data)
all.data <- SplitObject(all.data, split.by = "method")
all.data <- lapply(X = all.data, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- NormalizeData(x)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = all.data)
data.anchors <- FindIntegrationAnchors(object.list = all.data, anchor.features = features)
data.combined <- IntegrateData(anchorset = data.anchors)
dim(data.combined@assays$integrated@data)
data1.seurat=as.matrix(data.combined@assays$integrated@data[,1:length(ct.X1)])
data2.seurat=as.matrix(data.combined@assays$integrated@data[,(length(ct.X1)+1):(length(ct.X1)+length(ct.X2))])

data.int.seurat = t(cbind(data1.seurat,data2.seurat))

f.out[2,1:3] = faithful(data.int.seurat[1:length(ct.X1),], t(data.X1.smai))
f.out[2,4:6] = faithful(data.int.seurat[-c(1:length(ct.X1)),], t(data.X2.smai))


######################## Harmony
all.RNA = cbind(data.X1,data.X2)
meta.data= data.frame(method = c(rep("data1", dim(data.X1)[2]), rep("data2", dim(data.X2)[2])))
row.names(meta.data)=c(colnames(data.X1), colnames(data.X2))
#normalization and QC
all.data <- CreateSeuratObject(counts = all.RNA, project = "align", min.cells = 5, min.features = 50, meta.data = meta.data)
all.data <- NormalizeData(all.data)
all.data <- FindVariableFeatures(all.data, selection.method = "vst", nfeatures = 2000)
all.data <- ScaleData(all.data, verbose = FALSE)
all.data <- RunPCA(all.data, pc.genes = all.data@var.genes, npcs = 50, verbose = FALSE)
all.data <- RunHarmony(all.data, "method")
dim(all.data@reductions$harmony@cell.embeddings)
data1.har=t(as.matrix(all.data@reductions$harmony@cell.embeddings[1:length(ct.X1),]))
data2.har=t(as.matrix(all.data@reductions$harmony@cell.embeddings[(length(ct.X1)+1):(length(ct.X1)+length(ct.X2)),]))

data.int.har = t(cbind(data1.har,data2.har))

f.out[3,1:3] = faithful(data.int.har[1:length(ct.X1),], t(data.X1.smai))
f.out[3,4:6] = faithful(data.int.har[-c(1:length(ct.X1)),], t(data.X2.smai))


################liger
ifnb_liger <- createLiger(list(data1=(data.X1),data22=(data.X2)))
ifnb_liger <- rliger::normalize(ifnb_liger)
ifnb_liger <- selectGenes(ifnb_liger)
ifnb_liger <- scaleNotCenter(ifnb_liger)
ifnb_liger <- optimizeALS(ifnb_liger, k = 50)
ifnb_liger <- quantile_norm(ifnb_liger)
data1.liger=t(ifnb_liger@H$data1)
data2.liger=t(ifnb_liger@H$data2)

temp=factor(c(ct.X1,ct.X2))
names(temp)=names(ifnb_liger@clusters)
sum(names(temp)!=c(colnames(data.X1),colnames(data.X2)))
ifnb_liger@clusters = temp
cluster.results <- runWilcoxon(ifnb_liger, compare.method = "clusters")

#calcAlignment(ifnb_liger)
#[1] 0.15

data.int.liger = t(cbind(data1.liger,data2.liger))
f.out[4,1:3] = faithful(data.int.liger[1:length(ct.X1),], t(data.X1.smai))
f.out[4,4:6] = faithful(data.int.liger[-c(1:length(ct.X1)),], t(data.X2.smai))


#############fastMNN
library(batchelor)
out.fastmnn = fastMNN(data.X1.smai,data.X2.smai)
data1.fmnn=(as.matrix(assay(out.fastmnn)[,1:length(ct.X1)]))
data2.fmnn=(as.matrix(assay(out.fastmnn)[,(length(ct.X1)+1):(length(ct.X1)+length(ct.X2))]))
data.int.fmnn = t(cbind(data1.fmnn,data2.fmnn))
f.out[5,1:3] = faithful(data.int.fmnn[1:length(ct.X1),], t(data.X1.smai))
f.out[5,4:6] = faithful(data.int.fmnn[-c(1:length(ct.X1)),], t(data.X2.smai))


##############SCANORAMA
datasets <- list(t(data.X1.smai),t(data.X2.smai))
genes_list <- list((rownames(data.X1.smai)),(rownames(data.X2.smai)))
reticulate::use_python("/Users/rongma/opt/miniconda3/envs/scanorama/bin/python")
scanorama <- import('scanorama')
integrated.data <- scanorama$integrate(datasets, genes_list)
corrected.data <- scanorama$correct(datasets, genes_list, return_dense=TRUE)
integrated.corrected.data <- scanorama$correct(datasets, genes_list,
                                               return_dimred=TRUE, return_dense=TRUE)
data1.sca=t(integrated.corrected.data[[2]][[1]])
data2.sca=t(integrated.corrected.data[[2]][[2]])
data1.sca=data1.sca[match(rownames(data.X1.smai),integrated.corrected.data[[3]]),]
data2.sca=data2.sca[match(rownames(data.X1.smai),integrated.corrected.data[[3]]),]
data.int.sca = t(cbind(data1.sca,data2.sca))
f.out[6,1:3] = faithful(data.int.sca[1:length(ct.X1),], t(data.X1.smai))
f.out[6,4:6] = faithful(data.int.sca[-c(1:length(ct.X1)),], t(data.X2.smai))



############################
### Make plots
############################


#### gene analysis, integrated data
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


table(ct.X2)/length((ct.X2))


sim.gene.det.mat=c()

k=0
#repeat chunk below for various k
k=k+1
sample1.in=which(ct.X1 == levels(c(ct.X1,ct.X2))[k])
sample2.in=which(ct.X2 == levels(c(ct.X1,ct.X2))[k])


pval.X1=c()
pval.X2=c()
pval.smai=c()
pval.fmnn=c()
pval.sca=c()
pval.seu=c()

for(i in 1:dim(data.X1.smai)[1]){
  if(length(sample1.in)>0){
    pval.X1[i]=t.test(data.X1.smai[i,sample1.in], data.X1.smai[i,-sample1.in], alternative="greater")$p.value
  }
  if(length(sample2.in)>0){
    pval.X2[i]=t.test(data.X2.smai[i,sample2.in], data.X2.smai[i,-sample2.in], alternative="greater")$p.value
  }
  pval.fmnn[i] = t.test(c(data2.fmnn[i,sample2.in],data1.fmnn[i,sample1.in]),
                        c(data2.fmnn[i,-sample2.in],data1.fmnn[i,-sample1.in]), alternative="greater")$p.value
  pval.sca[i] = t.test(c(data2.sca[i,sample2.in],data1.sca[i,sample1.in]),
                       c(data2.sca[i,-sample2.in],data1.sca[i,-sample1.in]), alternative="greater")$p.value
  pval.seu[i] = t.test(c(data2.seurat[i,sample2.in],data1.seurat[i,sample1.in]),
                       c(data2.seurat[i,-sample2.in],data1.seurat[i,-sample1.in]), alternative="greater")$p.value
}
pval.liger=cluster.results[which(cluster.results$group==levels(factor(c(ct.X1,ct.X2)))[k]),]$pval
pval.liger=pval.liger[match(rownames(data.X1.smai),cluster.results[which(cluster.results$group==levels(factor(ct.X2))[k]),]$feature)]

gene.out = list()
alpha=0.01
if(length(sample1.in)>0){
  gene.out[[1]] = which(p.adjust(pval.X1, method="BH")<alpha)
}
if(length(sample2.in)>0){
  gene.out[[1]] = which(p.adjust(pval.X2, method="BH")<alpha)
}
gene.out[[2]] = which(p.adjust(pval.smai, method="BH")<alpha)
gene.out[[3]] = which(p.adjust(pval.fmnn, method="BH")<alpha)
gene.out[[4]] = which(p.adjust(pval.sca, method="BH")<alpha)
gene.out[[5]] = which(p.adjust(pval.seu, method="BH")<alpha)
gene.out[[6]] = which(p.adjust(pval.liger, method="BH")<alpha)

sim.gene.det = c()
for(i in 1:5){
  sim.gene.det[i] = jaccard(gene.out[[1]], gene.out[[i+1]])
}

sim.gene.det.mat = rbind(sim.gene.det.mat, sim.gene.det)
#repeat above chunk for various k

data.plot = data.frame(similarity=c(sim.gene.det.mat), method=rep(c("SMAI", "fastMNN", "Scanorama","Seurat","LIGER"),each=11))
data.plot = data.plot[-(1:11),]

pdf(file = paste0("markers_boxplot_Neg2.pdf"), width=3,height=3)
ggplot(data=data.plot,aes(x=method,y=similarity))+geom_boxplot(outlier.shape = NA)+ylim(0,1)+ylab("jaccard similarity")+
  geom_point(position = position_jitter(width = 0.2), color = "black") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14))
dev.off()

############ false discovery

fp <- function(a, b) {
  f.p = length(b)-length(intersect(a, b))
  t.p = length(b)
  return (f.p/t.p)
}

sim.gene.det.mat=c()

k=0
#repeat chunk below for various k
k=k+1
sample1.in=which(ct.X1 == levels(c(ct.X1,ct.X2))[k])
sample2.in=which(ct.X2 == levels(c(ct.X1,ct.X2))[k])


pval.X1=c()
pval.X2=c()
pval.smai=c()
pval.fmnn=c()
pval.sca=c()
pval.seu=c()

for(i in 1:dim(data.X1.smai)[1]){
  if(length(sample1.in)>0){
    pval.X1[i]=t.test(data.X1.smai[i,sample1.in], data.X1.smai[i,-sample1.in], alternative="greater")$p.value
  }
  if(length(sample2.in)>0){
    pval.X2[i]=t.test(data.X2.smai[i,sample2.in], data.X2.smai[i,-sample2.in], alternative="greater")$p.value
  }
  pval.fmnn[i] = t.test(c(data2.fmnn[i,sample2.in],data1.fmnn[i,sample1.in]),
                        c(data2.fmnn[i,-sample2.in],data1.fmnn[i,-sample1.in]), alternative="greater")$p.value
  pval.sca[i] = t.test(c(data2.sca[i,sample2.in],data1.sca[i,sample1.in]),
                       c(data2.sca[i,-sample2.in],data1.sca[i,-sample1.in]), alternative="greater")$p.value
  pval.seu[i] = t.test(c(data2.seurat[i,sample2.in],data1.seurat[i,sample1.in]),
                       c(data2.seurat[i,-sample2.in],data1.seurat[i,-sample1.in]), alternative="greater")$p.value
}
pval.liger=cluster.results[which(cluster.results$group==levels(factor(c(ct.X1,ct.X2)))[k]),]$pval
pval.liger=pval.liger[match(rownames(data.X1.smai),cluster.results[which(cluster.results$group==levels(factor(ct.X2))[k]),]$feature)]

gene.out = list()
alpha=0.01
gene.out[[1]] = union(which(p.adjust(pval.X1, method="BH")<alpha),which(p.adjust(pval.X2, method="BH")<alpha))
gene.out[[2]] = which(p.adjust(pval.smai, method="BH")<alpha)
gene.out[[3]] = which(p.adjust(pval.fmnn, method="BH")<alpha)
gene.out[[4]] = which(p.adjust(pval.sca, method="BH")<alpha)
gene.out[[5]] = which(p.adjust(pval.seu, method="BH")<alpha)
gene.out[[6]] = which(p.adjust(pval.liger, method="BH")<alpha)

sim.gene.det = c()
for(i in 1:5){
  sim.gene.det[i] = fp(gene.out[[1]], gene.out[[i+1]])
}

sim.gene.det.mat = rbind(sim.gene.det.mat, sim.gene.det)
#repeat above chunk

data.plot = data.frame(similarity=c(sim.gene.det.mat), method=rep(c("SMAI", "fastMNN", "Scanorama","Seurat","LIGER"),each=11))
data.plot = data.plot[-(1:11),]
pdf(file = paste0("FD_boxplot_Neg2.pdf"), width=3,height=3)
ggplot(data=data.plot,aes(x=method,y=similarity))+geom_boxplot(outlier.shape = NA)+ylim(-0.01,1.01)+ylab("False Positive Rate")+
  geom_point(position = position_jitter(width = 0.2), color = "black") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) #+geom_label_repel(aes(label = cell_type),
dev.off()

############ power

pow <- function(a, b) {
  t.d = length(intersect(a, b))
  t.p = length(a)
  return (t.d/t.p)
}


sim.gene.det.mat=c()
k=0
#repeat chunk below
k=k+1
sample1.in=which(ct.X1 == levels(c(ct.X1,ct.X2))[k])
sample2.in=which(ct.X2 == levels(c(ct.X1,ct.X2))[k])

pval.X1=c()
pval.X2=c()
pval.smai=c()
pval.fmnn=c()
pval.sca=c()
pval.seu=c()
pval.liger=c()
for(i in 1:dim(data.X1.smai)[1]){
  if(length(sample1.in)>0){
    pval.X1[i]=t.test(data.X1.smai[i,sample1.in], data.X1.smai[i,-sample1.in], alternative="greater")$p.value
  }
  if(length(sample2.in)>0){
    pval.X2[i]=t.test(data.X2.smai[i,sample2.in], data.X2.smai[i,-sample2.in], alternative="greater")$p.value
  }
  pval.fmnn[i] = t.test(c(data2.fmnn[i,sample2.in],data1.fmnn[i,sample1.in]),
                        c(data2.fmnn[i,-sample2.in],data1.fmnn[i,-sample1.in]), alternative="greater")$p.value
  pval.sca[i] = t.test(c(data2.sca[i,sample2.in],data1.sca[i,sample1.in]),
                       c(data2.sca[i,-sample2.in],data1.sca[i,-sample1.in]), alternative="greater")$p.value
  pval.seu[i] = t.test(c(data2.seurat[i,sample2.in],data1.seurat[i,sample1.in]),
                       c(data2.seurat[i,-sample2.in],data1.seurat[i,-sample1.in]), alternative="greater")$p.value
}
pval.liger=cluster.results[which(cluster.results$group==k),]$pval
pval.liger=pval.liger[match(rownames(data.X1.smai),cluster.results[which(cluster.results$group==levels(factor(ct.X2))[k]),]$feature)]

gene.out = list()
alpha=0.01
gene.out[[1]] = union(which(p.adjust(pval.X1, method="BH")<alpha),which(p.adjust(pval.X2, method="BH")<alpha))
gene.out[[2]] = which(p.adjust(pval.smai, method="BH")<alpha)
gene.out[[3]] = which(p.adjust(pval.fmnn, method="BH")<alpha)
gene.out[[4]] = which(p.adjust(pval.sca, method="BH")<alpha)
gene.out[[5]] = which(p.adjust(pval.seu, method="BH")<alpha)
gene.out[[6]] = which(p.adjust(pval.liger, method="BH")<alpha)

sim.gene.det = c()
for(i in 1:5){
  sim.gene.det[i] = pow(gene.out[[1]], gene.out[[i+1]])
}

sim.gene.det.mat = rbind(sim.gene.det.mat, sim.gene.det)

data.plot = data.frame(similarity=c(sim.gene.det.mat), method=rep(c("SMAI", "fastMNN", "Scanorama","Seurat","LIGER"),each=11))
data.plot = data.plot[-(1:11),]
pdf(file = paste0("power_boxplot_Neg2.pdf"), width=3,height=3)
ggplot(data=data.plot,aes(x=method,y=similarity))+geom_boxplot(outlier.shape = NA)+ylim(-0.01,1.01)+ylab("Power")+
  geom_point(position = position_jitter(width = 0.2), color = "black") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) #+geom_label_repel(aes(label = cell_type),
dev.off()


############ distortion bar plot
method.names = c("Seurat","Harmony", "LIGER", "fastMNN", "Scanorama")
data.plot = data.frame(similarity = (f.out[-1,1]+f.out[-1,4])/2, method = method.names)
pdf(file = paste0("Neg2_dist.pdf"), width=3,height=3)
ggplot(data=data.plot,aes(x=method,y=similarity))+geom_bar(stat="identity")+ylim(0,1)+ylab("correlation")+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14))
dev.off()
