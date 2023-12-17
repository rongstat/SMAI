library(data.table)
library(uwot)
library(Seurat)
library(ggrepel)
library(ggplot2)
library(cluster)
library(uwot)
library(fossil)
library(batchelor)
library(harmony)
library(rliger)
library(reticulate)
library(SMAI)
library(Matrix)


########################################
######### load data and preprocessing
########################################

meta.data = read.table("data/GSE178429_PBMCs_stim_scRNAseq_cellMeta.txt", header=T)
data = readMM("data/PBMCs_stim_scRNAseq_counts.txt")
dim(data)
gene.name = read.table("~/Dropbox/Reading/working/Jason/data/GSE178429_PBMCs_stim_scRNAseq_geneNames.txt")
row.names(data) = unlist(gene.name)

table(meta.data$Condition)

data.X1 = as.matrix(data[,which(meta.data$Condition=="Control_1h")])
data.X2 = as.matrix(data[,which(meta.data$Condition=="LPS_1h")])
dim(data.X1)
dim(data.X2)


colnames(data.X1) = meta.data$cellBarcode[which(meta.data$Condition=="Control_1h")]
colnames(data.X2) = meta.data$cellBarcode[which(meta.data$Condition=="LPS_1h")]
row.names(data.X1) = unlist(gene.name)
row.names(data.X2) = unlist(gene.name)
all.RNA = cbind(data.X1, data.X2)
dim(all.RNA)
meta.data=data.frame( method = factor(c(rep("data1",dim(data.X1)[2]), rep("data2",dim(data.X2)[2]))))
ct.X1 = rep("data1",dim(data.X1)[2])
ct.X2 = rep("data2",dim(data.X2)[2])
row.names(meta.data)=c(colnames(all.RNA))
all.data <- CreateSeuratObject(counts = all.RNA, project = "align", min.cells = 3,
                               min.features = 10, meta.data = meta.data)
all.data <- SplitObject(all.data, split.by = "method")
all.data <- lapply(X = all.data, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- NormalizeData(x)
})
features <- SelectIntegrationFeatures(object.list = all.data)
data.X1.smai = as.matrix(all.data$data1@assays$RNA$data)
data.X2.smai = as.matrix(all.data$data2@assays$RNA$data)
data.X1.smai = data.X1.smai[match(features,rownames(data.X1.smai)),]
data.X2.smai = data.X2.smai[match(features,rownames(data.X2.smai)),]

##########save pre-processed data

write.table(data.X1.smai, file="pos7_data1.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(data.X2.smai, file="pos7_data2.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


############## SMAI
#identify maximal correspondence set
data.anchors <- FindIntegrationAnchors(object.list = all.data, anchor.features = features)
data.combined <- IntegrateData(anchorset = data.anchors)
dim(data.combined@assays$integrated@data)
data1=as.matrix(data.combined@assays$integrated@data[,1:length(ct.X1)])
data2=as.matrix(data.combined@assays$integrated@data[,(length(ct.X1)+1):(length(ct.X1)+length(ct.X2))])
prop.palign = 0.5
#mutual nearest neighbor matching
k1=10
prop.mc=0
while(prop.mc<prop.palign/2+0.01){
  mnn.out = findMutualNN(t(data1),t(data2), k1=k1) 
  prop.mc = min(c(length(unique(mnn.out$first))/dim(data1)[2],
                  length(unique(mnn.out$second))/dim(data2)[2]))
  k1=k1+5
}

#start SMAI
set.seed(9)
out.smai <- align(data.X1.smai, data.X2.smai, sel1=unique(mnn.out$first),  sel2=unique(mnn.out$second),
                  t=3, knn=30, r.max=200, dir.map = "auto",  denoise="scree",
                  prop.align=prop.palign/2, cutoff = 1.001, cutoff.t=1.5)
out.smai$p.value
#[1] 0.4704584

meta_data=data.frame(batch = factor(c(rep("X1",length(ct.X1)), rep("X2",length(ct.X2)))))
data.int = t(cbind(out.smai$data1.integrate,out.smai$data2.integrate))

out=matrix(ncol=3,nrow=8)
set.seed(10)
out0=assess(data.int, meta_data)
out[1,] = unlist(out0)

f.out = matrix(ncol=6,nrow=8)
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
out0=assess(data.int.seurat, meta_data)
out[2,] = unlist(out0)


######################## Harmony
library(harmony)
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

out0=assess(data.int.har, meta_data)
out[3,] = unlist(out0)

################liger
library(rliger)
ifnb_liger <- createLiger(list(data1=(data.X1),data2=(data.X2)))
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


out0=assess(data.int.liger, meta_data)
out[4,] = unlist(out0)


#############fastMNN
library(batchelor) 
out.fastmnn = fastMNN(data.X1.smai,data.X2.smai)
data1.fmnn=(as.matrix(assay(out.fastmnn)[,1:length(ct.X1)]))
data2.fmnn=(as.matrix(assay(out.fastmnn)[,(length(ct.X1)+1):(length(ct.X1)+length(ct.X2))]))
data.int.fmnn = t(cbind(data1.fmnn,data2.fmnn))
f.out[5,1:3] = faithful(data.int.fmnn[1:length(ct.X1),], t(data.X1.smai))
f.out[5,4:6] = faithful(data.int.fmnn[-c(1:length(ct.X1)),], t(data.X2.smai))

out0=assess(data.int.fmnn, meta_data)
out[5,] = unlist(out0)

##############SCANORAMA
datasets <- list(t(data.X1.smai),t(data.X2.smai))
genes_list <- list((rownames(data.X1.smai)),(rownames(data.X2.smai)))
library(reticulate)
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

out0=assess(data.int.sca, meta_data)
out[6,] = unlist(out0)

save(out.smai, data.X1.smai, data.X2.smai, 
     data1.seurat, data2.seurat,
     data1.har, data2.har,
     data1.liger, data2.liger,
     data1.fmnn, data2.fmnn,
     data1.sca, data2.sca,cluster.results,
     file="pos7_align.RData")

load("pos7_align.RData")
################Pamona

#run pamona_run.py in terminal

data1.pamona = t(read.csv("~/Dropbox/Eric-Rong-James/Alignment/RCodes/revision/pamona/pos7_data1.csv", header=F))
data2.pamona = t(read.csv("~/Dropbox/Eric-Rong-James/Alignment/RCodes/revision/pamona/pos7_data2.csv", header=F))
data.int.pamona = t(cbind(data1.pamona,data2.pamona))
f.out[7,1:3] = faithful(data.int.pamona[1:length(ct.X1),], t(data.X1.smai))
f.out[7,4:6] = faithful(data.int.pamona[-c(1:length(ct.X1)),], t(data.X2.smai))
out0=assess(data.int.pamona, meta_data)
out[7,] = unlist(out0)

################scot

#run scot_run.py in terminal

data1.scot = t(read.csv("~/Dropbox/Eric-Rong-James/Alignment/RCodes/revision/scot/pos7_data1.csv", header=F))
data2.scot = t(read.csv("~/Dropbox/Eric-Rong-James/Alignment/RCodes/revision/scot/pos7_data2.csv", header=F))
f.out[8,4:6] = faithful(t(data2.scot), t(data.X2.smai))
f.out[8,1:3] = faithful(t(data1.scot), t(data.X1.smai))
data.int.scot = t(cbind(data1.scot,data2.scot))
out0=assess(data.int.scot, meta_data)
out[8,] = unlist(out0)

############################
### Make plots
############################

## db vs. pearson
method.names = c("SMAI", "Seurat","Harmony", "LIGER", "fastMNN", "Scanorama","pamona","SCOT")
data.plot = data.frame(corr = (f.out[,1]+f.out[,4])/2, ind.bat = log(f.out[,3]), method = method.names)
highlight_df = data.plot[1,]
pdf(file = "Figs/Pos7_db&corr.pdf",  width=3.5,height=3.5)
sp <- ggplot(data.plot, aes(x=ind.bat, y=corr)) +
  geom_point(size=2) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) +
  geom_point(data=highlight_df,
             aes(x=ind.bat, y=corr),
             color='red',
             size=3) +
  ylab("pearson correlation") + xlab("log(batch D-B index)")+geom_label_repel(aes(label = method),
                                                                              force         = 180,
                                                                              max.overlaps = 10,
                                                                              segment.color = 'grey50')
sp
dev.off()

# ch vs. spearman
data.plot = data.frame(corr = (f.out[,2]+f.out[,5])/2, ind.bat = f.out[,2], method = method.names)
highlight_df = data.plot[1,]
pdf(file = "Pos7_ch&sp.pdf",  width=3.5,height=3.5)
sp <- ggplot(data.plot, aes(x=ind.bat, y=corr)) +
  geom_point(size=2) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) +
  geom_point(data=highlight_df,
             aes(x=ind.bat, y=corr),
             color='red',
             size=3) +
  ylab("spearman's rho") + xlab("1/(batch C-H index)")+geom_label_repel(aes(label = method),
                                                                        force         = 70,
                                                                        max.overlaps = 30,
                                                                        segment.color = 'grey50')
sp
dev.off()

######### UMAP plot
data.int = t(cbind(out.smai$data1.integrate,out.smai$data2.integrate))
umap.out = umap(data.int, n_neighbors = 50, metric = "cosine", spread = 5)
data.plot = data.frame(UMAP1 = c(umap.out[,1]), UMAP2= c(umap.out[,2]),
                       cell_type = factor(meta_data$cell_type), batch = factor(meta_data$batch))
p1 <- ggplot(data.plot, aes(x=UMAP1, y=UMAP2, colour = cell_type)) + geom_point(size = 0.7)
p1 <- p1 +
  theme(axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14,face="bold"), #change legend title font size
        legend.text = element_text(size=14,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=5)))#9*7 under pdf
pdf(file = "Pos7_SMAI1.pdf", width=8,height=6)
p1
dev.off()

p1 <- ggplot(data.plot, aes(x=UMAP1, y=UMAP2, colour = batch, alpha=0.5)) + geom_point(size = 0.7)
p1 <- p1 +
  theme(axis.text=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=14,face="bold"), #change legend title font size
        legend.text = element_text(size=14,face="bold")) +
  guides(colour = guide_legend(override.aes = list(size=5)))#7*6 under pdf
pdf(file = "Pos7_SMAI2.pdf", width=7,height=6)
p1
dev.off()
