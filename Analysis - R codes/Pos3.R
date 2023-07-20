#######################################################################
###############  real data analysis: Task Pos3
###############     July 14, 2023
#######################################################################

library(Seurat)
library(SeuratData)
library(uwot)
library(ggrepel)
library(ggplot2)
library(cluster)
library(uwot)
library(fossil)
library(batchelor)
library(harmony)
library(rliger)
library(reticulate)

########################################
######### load data and preprocessing
########################################
# install dataset
InstallData("pbmcsca")

# load dataset
LoadData("pbmcsca")

# split the dataset into a list
data.list <- SplitObject(pbmcsca, split.by = "Method")

# normalize and identify variable features for each dataset independently
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = data.list)
data.X1 = as.matrix(data.list$`10x Chromium (v2)`@assays$RNA@counts)
ct.X1 = data.list$`10x Chromium (v2)`@meta.data$CellType
data.X2 = as.matrix(data.list$`10x Chromium (v3)`@assays$RNA@counts)
ct.X2 = data.list$`10x Chromium (v3)`@meta.data$CellType
table(ct.X1)
table(ct.X2)

########### feature selection
data.X1 = data.X1[match(features,rownames(data.X1)),]
data.X2 = data.X2[match(features,rownames(data.X2)),]
dim(data.X1)
dim(data.X2)



################################################
########### alignment with different methods
################################################


#preprocessing and normalization
all.RNA = cbind(data.X1, data.X2)
dim(all.RNA)
meta.data=data.frame(cell_type = factor(c(ct.X1,ct.X2)),
                     method = factor(c(rep("data1",length(ct.X1)), rep("data2",length(ct.X2)))))
row.names(meta.data)=c(colnames(all.RNA))
all.data <- CreateSeuratObject(counts = all.RNA, project = "align", min.cells = 1,
                               min.features = 1, meta.data = meta.data)
all.data <- SplitObject(all.data, split.by = "method")
all.data <- lapply(X = all.data, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- NormalizeData(x)
})
features <- SelectIntegrationFeatures(object.list = all.data)
data.X1.smai = as.matrix(all.data$data1@assays$RNA@data)
data.X2.smai = as.matrix(all.data$data2@assays$RNA@data)
data.X1.smai = data.X1.smai[match(features,rownames(data.X1.smai)),]
data.X2.smai = data.X2.smai[match(features,rownames(data.X2.smai)),]

############## SMAI
#identify maximal correspondence subsets
data.anchors <- FindIntegrationAnchors(object.list = all.data, anchor.features = features)
data.combined <- IntegrateData(anchorset = data.anchors)
dim(data.combined@assays$integrated@data)
data1=as.matrix(data.combined@assays$integrated@data[,1:length(ct.X1)])
data2=as.matrix(data.combined@assays$integrated@data[,(length(ct.X1)+1):(length(ct.X1)+length(ct.X2))])
mnn.out = findMutualNN(t(data1),t(data2), k1=15)
prop.align.par = min(c(length(unique(mnn.out$first))/dim(data1)[2],
                       length(unique(mnn.out$second))/dim(data2)[2]))
prop.align.par

#start SMAI
set.seed(8)
out.smai <- align(data.X1.smai, data.X2.smai, sel1=unique(mnn.out$first),  sel2=unique(mnn.out$second),
                  t=3, knn=30, r.max=200, dir.map = "auto",  denoise="scree",
                  prop.align=0.3,  cutoff = 1.001, cutoff.t=1.5)
out.smai$p.value
#[1] 0.7851726

meta_data=data.frame(cell_type = factor(c(ct.X1,ct.X2)),
                     batch = factor(c(rep("X1",length(ct.X1)), rep("X2",length(ct.X2)))))
data.int = t(cbind(out.smai$data1.integrate,out.smai$data2.integrate))

out=matrix(ncol=3,nrow=6)
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
all.data <- CreateSeuratObject(counts = all.RNA, project = "panc8", min.cells = 1, min.features = 1,
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
all.data <- CreateSeuratObject(counts = all.RNA, project = "panc8", min.cells = 1, min.features = 1, meta.data = meta.data)
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
data1.lig=t(ifnb_liger@H$data1)
data2.lig=t(ifnb_liger@H$data2)

temp=factor(c(ct.X1,ct.X2))
names(temp)=names(ifnb_liger@clusters)
sum(names(temp)!=c(colnames(data.X1),colnames(data.X2)))
ifnb_liger@clusters = temp
cluster.results <- runWilcoxon(ifnb_liger, compare.method = "clusters")


data.int.lig = t(cbind(data1.lig,data2.lig))
f.out[4,1:3] = faithful(data.int.lig[1:length(ct.X1),], t(data.X1.smai))
f.out[4,4:6] = faithful(data.int.lig[-c(1:length(ct.X1)),], t(data.X2.smai))

out0=assess(data.int.lig, meta_data)
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
colnames(out)=names(unlist(out0))



############################
### Make plots
############################

## db vs. kendall
method.names = c("SMAI", "Seurat","Harmony", "LIGER", "fastMNN", "Scanorama")
data.plot = data.frame(corr = (f.out[,1]+f.out[,4])/2, ind.bat = log(out[,3]), method = method.names)
highlight_df = data.plot[1,]
pdf(file = "Pos3_db&k.pdf", width=3.5,height=3.5)
ggplot(data.plot, aes(x=ind.bat, y=corr)) +
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
  ylab("Kendall's tau") + xlab("log(batch D-B index)")+geom_label_repel(aes(label = method),
                                                                        force         = 70,
                                                                        max.overlaps = 30,
                                                                        segment.color = 'grey50')
dev.off()




# ch vs. spearman
data.plot = data.frame(corr = (f.out[,2]+f.out[,5])/2, ind.bat = out[,2], method = method.names)
highlight_df = data.plot[1,]
pdf(file = "Pos3_ch&sp.pdf",  width=3.5,height=3.5)
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

##eigenvalue plots

out=eigen(data.X1.smai %*% t(data.X1.smai))
plot(log(out$values[1:50]), ylab="log(eigenvalues)")#4.5*5.5

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
pdf(file = "Pos3_SMAI1.pdf", width=8,height=6)
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
pdf(file = "Pos3_SMAI2.pdf", width=7,height=6)
p1
dev.off()



#### gene analysis
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


sim.gene.det.mat=c()
marker.g=c()
k=0
#repeat following chunk for various k
k=k+1
sample1.in=which(ct.X1 == levels(factor(ct.X2))[k])
sample2.in=which(ct.X2 == levels(factor(ct.X2))[k])

pval.X1=c()
pval.X2=c()
pval.smai=c()
pval.fmnn=c()
pval.sca=c()
pval.seu=c()

for(i in 1:dim(data.X1.smai)[1]){
  pval.X1[i]=t.test(data.X1.smai[i,sample1.in], data.X1.smai[i,-sample1.in], alternative="greater")$p.value
  pval.X2[i]=t.test(data.X2.smai[i,sample2.in], data.X2.smai[i,-sample2.in], alternative="greater")$p.value
  pval.smai[i] = t.test(c(out.smai$data2.integrate[i,sample2.in],out.smai$data1.integrate[i,sample1.in]),
                        c(out.smai$data2.integrate[i,-sample2.in],out.smai$data1.integrate[i,-sample1.in]), alternative="greater")$p.value
  pval.fmnn[i] = t.test(c(data2.fmnn[i,sample2.in],data1.fmnn[i,sample1.in]),
                        c(data2.fmnn[i,-sample2.in],data1.fmnn[i,-sample1.in]), alternative="greater")$p.value
  pval.sca[i] = t.test(c(data2.sca[i,sample2.in],data1.sca[i,sample1.in]),
                       c(data2.sca[i,-sample2.in],data1.sca[i,-sample1.in]), alternative="greater")$p.value
  pval.seu[i] = t.test(c(data2.seurat[i,sample2.in],data1.seurat[i,sample1.in]),
                       c(data2.seurat[i,-sample2.in],data1.seurat[i,-sample1.in]), alternative="greater")$p.value

}

pval.liger=cluster.results[which(cluster.results$group==levels(factor(ct.X2))[k]),]$pval
pval.liger=pval.liger[match(rownames(data.X1.smai),cluster.results[which(cluster.results$group==levels(factor(ct.X2))[k]),]$feature)]


gene.out = list()
alpha=0.01
gene.out[[1]] = union(which(p.adjust(pval.X1, method="BH")<alpha),which(p.adjust(pval.X2, method="BH")<alpha))
#gene.out[[1]] = intersect(which(p.adjust(pval.X1, method="BH")<alpha),which(p.adjust(pval.X2, method="BH")<alpha))
#marker.g=c(marker.g, gene.out[[1]])
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

data.plot = data.frame(similarity=c(sim.gene.det.mat), method=rep(c("SMAI", "fastMNN", "Scanorama","Seurat","LIGER"),each=8))

pdf(file = paste0("markers_boxplot_Pos3.pdf"), width=3,height=3)
ggplot(data=data.plot,aes(x=method,y=similarity))+geom_boxplot(outlier.shape = NA)+ylim(0,1)+ylab("jaccard similarity")+
  geom_point(position = position_jitter(width = 0.2), color = "black") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) #+geom_label_repel(aes(label = cell_type),
dev.off()


#### additional more detailed gene analysis

i=1696
k=7
#ori
data.plot=data.frame(expr=log(1+c(data.X1.smai[i,])),
                     cell_type = as.character(ct.X1))
data.plot$cell_type[which(data.plot$cell_type!= levels(factor(ct.X2))[k])]="Others"


pdf(file = paste0("Pos3_",levels(factor(ct.X2))[k],"_",rownames(data.X1.smai)[i],"_ori1.pdf"), width=2,height=3.4)
ggplot(data=data.plot,aes(x=cell_type,y=expr))+geom_boxplot()+ylab("log(expression level)")+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) #+geom_label_repel(aes(label = cell_type),
dev.off()

data.plot=data.frame(expr=log(1+c(data.X2.smai[i,])),
                     cell_type = as.character(ct.X2))
data.plot$cell_type[which(data.plot$cell_type!= levels(factor(ct.X2))[k])]="Others"

pdf(file = paste0("Pos3_",levels(factor(ct.X2))[k],"_",rownames(data.X1.smai)[i],"_ori2.pdf"), width=2,height=3.4)
ggplot(data=data.plot,aes(x=cell_type,y=expr))+geom_boxplot()+ylab("log(expression level)")+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) #+geom_label_repel(aes(label = cell_type),
dev.off()


#smai
data.plot=data.frame(expr=log(3+c(out.smai$data1.integrate[i,],out.smai$data2.integrate[i,])),
                     cell_type = as.character(c(ct.X1,ct.X2)) )
data.plot$cell_type[which(data.plot$cell_type!= levels(factor(ct.X2))[k])]="Others"

pdf(file = paste0("Pos3_",levels(factor(ct.X2))[k],"_",rownames(data.X1.smai)[i],"_smai.pdf"), width=2,height=3.4)
ggplot(data=data.plot,aes(x=cell_type,y=expr))+geom_boxplot()+ylab("log(expression level)")+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) #+geom_label_repel(aes(label = cell_type),
dev.off()

#seura
data.plot=data.frame(expr=log(2+c(data1.seurat[i,],data2.seurat[i,])),
                     cell_type =  as.character(c(ct.X1,ct.X2)) )
data.plot$cell_type[which(data.plot$cell_type!= levels(factor(ct.X2))[k])]="Others"

pdf(file = paste0("task4_",levels(factor(ct.X2))[k],"_",rownames(data.X1.smai)[i],"_seurat.pdf"), width=2,height=3.4)
ggplot(data=data.plot,aes(x=cell_type,y=expr))+geom_boxplot()+ylab("log(expression level)")+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) #+geom_label_repel(aes(label = cell_type),
dev.off()

#fmnn
data.plot=data.frame(expr=log(2+c(data1.fmnn[i,],data2.fmnn[i,])),
                     cell_type =  as.character(c(ct.X1,ct.X2)) )
data.plot$cell_type[which(data.plot$cell_type!= levels(factor(ct.X2))[k])]="Others"

pdf(file = paste0("Pos3_",levels(factor(ct.X2))[k],"_",rownames(data.X1.smai)[i],"_fmnn.pdf"), width=2,height=3.4)
ggplot(data=data.plot,aes(x=cell_type,y=expr))+geom_boxplot()+ylab("log(expression level)")+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) #+geom_label_repel(aes(label = cell_type),
dev.off()

#scanorama
data.plot=data.frame(expr=log(1+c(data1.sca[i,],data2.sca[i,])),
                     cell_type =  as.character(c(ct.X1,ct.X2)) )
data.plot$cell_type[which(data.plot$cell_type!= levels(factor(ct.X2))[k])]="Others"

pdf(file = paste0("Pos3_",levels(factor(ct.X2))[k],"_",rownames(data.X1.smai)[i],"_sca.pdf"), width=2,height=3.4)
ggplot(data=data.plot,aes(x=cell_type,y=expr))+geom_boxplot()+ylab("log(expression level)")+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) #+geom_label_repel(aes(label = cell_type),
dev.off()

#### display batch effect

highlight = names(out.smai$gamma)[abs(out.smai$gamma)>quantile(abs(out.smai$gamma),0.995)]
data.plot=data.frame(gamma=out.smai$gamma, label=names(out.smai$gamma))

library(tidyverse)
library(cowplot)
data.plot <- data.plot %>%
  mutate(
    label = ifelse(label %in% highlight, label, "")
  )

data.plot=data.plot[order(out.smai$gamma),]
plt <- ggplot(data.plot, aes(x=1:dim(data.X1.smai)[1], y = gamma)) +
  geom_point(
    size = 1.5,
    alpha = 0.8,
    aes(colour= factor((label=="")))# It's nice to add some transparency because there may be overlap.
  ) + ylab("gamma")+xlab("index")+
  geom_text_repel(
    aes(label = label),
    force         = 20,
    max.overlaps = 10,
    segment.color = 'grey50')
pdf(file = "Pos3_gam.pdf", width=5,height=4)
plt + theme(legend.position = "none")
dev.off()

##### display genes

i=which(rownames(data.X1.smai)=="SLC31A2")

data.plot=data.frame(gene=data.X1.smai[i,], cell_type=meta_data$cell_type[1:length(ct.X1)])
sp <- ggplot(data.plot, aes(x=cell_type, y=gene)) +
  geom_boxplot(outlier.size = 0.5) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) + ylim(0,8.5)+
  ylab(names(out.smai$gamma)[i])
pdf(file = paste0("Pos3_markers_",rownames(data.X1.smai)[i],".pdf"), width=3,height=4)
sp
dev.off()


data.plot=data.frame(gene=data.X2.smai[i,], cell_type=ct.X2)
sp <- ggplot(data.plot, aes(x=cell_type, y=gene)) +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) +ylim(0,8.5)+
  geom_boxplot(outlier.size = 0.5) +
  ylab(names(out.smai$gamma)[i])
pdf(file = paste0("Pos3_markers_",rownames(data.X1.smai)[i],"2.pdf"), width=3,height=4)
sp
dev.off()


#### ROC

id.marker=rep(0,dim(data.X1.smai)[1])
unique(marker.g) #set alpha=0.001 and use intersection rather than union for selecting p-values
id.marker[unique(marker.g)]=1

library(pROC)
roc_score=roc(id.marker, abs(out.smai$gamma)) #AUC score
pdf(file = "Pos3_roc1.pdf", width=4,height=4)
plot(roc_score ,main ="ROC curve - gamma")#4*4
dev.off()
roc_score$auc
#Area under the curve: 0.88


library(pROC)
temp.mat=abs(out.smai$feature.rot)
diag(temp.mat)=0
roc_score=roc(id.marker, rowSums(temp.mat)) #AUC score
pdf(file = "Pos3_roc2.pdf", width=4,height=4)
plot(roc_score ,main ="ROC curve - R")#4*4
dev.off()
roc_score$auc
#Area under the curve: 0.88




