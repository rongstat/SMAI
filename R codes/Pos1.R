#######################################################################
###############  data analysis: Task Pos1
###############     July 13, 2023
#######################################################################


library(Seurat)
library(SeuratData)
library(uwot)
library(ggrepel)
library(ggplot2)
library(cluster)
library(fossil)
library(batchelor)
library(harmony)
library(rliger)
library(reticulate)
library(SMAI)
########################################
######### load data and preprocessing
########################################

# install dataset
InstallData("panc8")
# load dataset
LoadData("panc8")
# split the dataset into a list
data.list <- SplitObject(panc8, split.by = "tech")

# normalize and identify variable features for each dataset independently
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = data.list)
table(data.frame(panc8@meta.data$celltype,panc8@meta.data$tech))
data.X1 = as.matrix(data.list$smartseq2@assays$RNA@counts)
ct.X1 = data.list$smartseq2@meta.data$celltype
data.X2 = as.matrix(data.list$celseq2@assays$RNA@counts)
ct.X2 = data.list$celseq2@meta.data$celltype

# select cell types of moderate sizes
keep.ind = as.logical(match(ct.X1,names(table(ct.X1))[which(table(ct.X1)>20)]))
keep.ind[is.na(keep.ind)]=FALSE
data.X1=data.X1[,keep.ind]
ct.X1=ct.X1[keep.ind]
dim(data.X1)
table(ct.X1)
table(ct.X1)/length(ct.X1)

keep.ind = as.logical(match(ct.X2,names(table(ct.X2))[which(table(ct.X2)>20)]))
keep.ind[is.na(keep.ind)]=FALSE
data.X2=data.X2[,keep.ind]
ct.X2=ct.X2[keep.ind]
dim(data.X2)
table(ct.X2)
table(ct.X2)/length(ct.X2)

data.X1 = data.X1[match(features,rownames(data.X1)),]
dim(data.X1)
table(ct.X1)/length(ct.X1)
data.X2 = data.X2[match(features,rownames(data.X2)),]
dim(data.X2)
table(ct.X2)
table(ct.X2)/length(ct.X2)

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
#identify maximal correspondence subsets
mnn.out = findMutualNN(t(data1),t(data2), k1=20)
prop.align.par = min(c(length(unique(mnn.out$first))/dim(data1)[2],
                       length(unique(mnn.out$second))/dim(data2)[2]))
prop.align.par

#start SMAI
set.seed(3)
out.smai <- align(data.X1.smai, data.X2.smai, dir.map = "auto", denoise="scree",
                  sel1=unique(mnn.out$first),  sel2=unique(mnn.out$second),
                  t=5, knn=30, r.max=200, cutoff = 1.001, cutoff.t=1.5,
                  prop.align = 0.3)

out.smai$p.value#p-value from SMAI-test
#[1] 0.4394436
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
all.data <- CreateSeuratObject(counts = all.RNA, project = "panc8", min.cells = 3, min.features = 10,
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
all.RNA = cbind(data.X1,data.X2)
meta.data= data.frame(method = c(rep("data1", dim(data.X1)[2]), rep("data2", dim(data.X2)[2])))
row.names(meta.data)=c(colnames(data.X1), colnames(data.X2))
all.data <- CreateSeuratObject(counts = all.RNA, project = "panc8", min.cells = 3, min.features = 10, meta.data = meta.data)
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
ifnb_liger <- createLiger(list(data1=(data.X1),data2=(data.X2)))
ifnb_liger <- rliger::normalize(ifnb_liger)
ifnb_liger <- selectGenes(ifnb_liger)
ifnb_liger <- scaleNotCenter(ifnb_liger)
ifnb_liger <- optimizeALS(ifnb_liger, k = 50)
ifnb_liger <- quantile_norm(ifnb_liger)
data1.lig=t(ifnb_liger@H$data1)
data2.lig=t(ifnb_liger@H$data2)

temp=factor(as.numeric(factor(c(ct.X1,ct.X2))))
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
pdf(file = "Pos1_db&k.pdf", width=3.5,height=3.5)
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
pdf(file = "Pos1_ch&sp.pdf",  width=3.5,height=3.5)
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
pdf(file = "Pos1_SMAI1.pdf", width=8,height=6)
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
pdf(file = "Pos1_SMAI2.pdf", width=7,height=6)
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
#repeat the chunk below for various k
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

pval.liger=cluster.results[which(cluster.results$group==k),]$pval
pval.liger=pval.liger[match(rownames(data.X1.smai),cluster.results[which(cluster.results$group==k),]$feature)]

gene.out = list()
alpha=0.01

 gene.out[[1]] = union(which(p.adjust(pval.X1, method="BH")<alpha),which(p.adjust(pval.X2, method="BH")<alpha))
   #gene.out[[1]] = intersect(which(p.adjust(pval.X1, method="BH")<alpha),which(p.adjust(pval.X2, method="BH")<alpha))
   #marker.g = c(marker.g, gene.out[[1]])
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

#repeat the above chunk for various k

data.plot = data.frame(similarity=c(sim.gene.det.mat), method=rep(c("SMAI", "fastMNN", "Scanorama","Seurat","LIGER"),each=8))

pdf(file = paste0("markers_boxplot_Pos1.pdf"), width=3,height=3)
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


##### display genes

i=which(rownames(data.X1.smai)=="REG1A")

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
pdf(file = paste0("Pos1-1_markers_",rownames(data.X1.smai)[i],".pdf"), width=3,height=4)
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
pdf(file = paste0("Pos1-2_markers_",rownames(data.X1.smai)[i],"2.pdf"), width=3,height=4)
sp
dev.off()

boxplot(data.X1.smai[i,],data.X2.smai[i,])


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
plt <- ggplot(data.plot, aes(x=1:1987, y = gamma)) +
  geom_point(
    size = 1.5,
    alpha = 0.8,
    aes(colour= factor((label=="")))# It's nice to add some transparency because there may be overlap.
  ) + ylab("gamma")+xlab("index")+
  geom_text_repel(
    aes(label = label),
    force         = 30,
    max.overlaps = 30,
    segment.color = 'grey50')
pdf(file = "Pos1_gam.pdf", width=4,height=6)
plt + theme(legend.position = "none")+coord_flip()
dev.off()

#### ROC

id.marker=rep(0,dim(data.X1.smai)[1])
unique(marker.g) #set alpha=0.001 and use intersection rather than union for selecting p-values
id.marker[unique(marker.g)]=1

library(pROC)
roc_score=roc(id.marker, abs(out.smai$gamma)) #AUC score
pdf(file = "Pos1_roc1.pdf", width=4,height=4)
plot(roc_score ,main ="ROC curve - gamma")#4*4
dev.off()
roc_score$auc
#Area under the curve: 0.7775


library(pROC)
temp.mat=abs(out.smai$feature.rot)
diag(temp.mat)=0
roc_score=roc(id.marker, rowSums(temp.mat)) #AUC score
pdf(file = "Pos1_roc2.pdf", width=4,height=4)
plot(roc_score ,main ="ROC curve - R")#4*4
dev.off()
roc_score$auc
#Area under the curve: 0.858


#### display rotation and genes


plot(diag(out.smai$feature.rot))
library(igraph)
rot.adj = out.smai$feature.rot
rownames(rot.adj)=rownames(data.X1.smai)
colnames(rot.adj)=rownames(data.X1.smai)
rot.adj[abs(rot.adj)<quantile(abs(rot.adj),0.9999)]=0
diag(rot.adj)=0
0.0001*dim(rot.adj)^2
rm.ind=intersect(which(rowSums(rot.adj!=0)==0), which(colSums(rot.adj!=0)==0))
rot.adj = rot.adj[-rm.ind,-rm.ind]

g1 = graph_from_adjacency_matrix(abs(rot.adj), mode="directed",weighted = TRUE)
pdf(file = "Pos1_rot.pdf", width=8,height=8)
plot(g1,vertex.color="lightblue", vertex.size=3,edge.width = E(g1)$weight*10,
     edge.arrow.size=0.4, vertex.label.cex=0.8, vertex.label.dist=0.2)
dev.off()

#select genes
a=c("HSPA1B")
b=c("HSPA1A")

data.plot=data.frame(gene1=data.X1.smai[which(rownames(data.X1.smai)==a),],
                     gene2=data.X1.smai[which(rownames(data.X1.smai)==b),])
data.plot=data.frame(gene1=out.smai$data1.integrate[which(rownames(data.X1.smai)==a),],
                     gene2=out.smai$data1.integrate[which(rownames(data.X1.smai)==b),])
sp <- ggplot(data.plot, aes(x=gene1, y=gene2)) +
  geom_point(size=1.5, alpha=0.6) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) +
  ylab(b) + xlab(a) +       xlim(0,6)+ylim(0,6)+
  stat_smooth(method = "lm",
              formula = y ~ 0+x,
              geom = "smooth") +annotate("text",x=2, y=6,label=(paste0("slope==",coef(lm(data.plot$gene2~0+data.plot$gene1))[1])),parse=TRUE)
sp


data.plot=data.frame(gene1=data.X2.smai[which(rownames(data.X2.smai)==a),],
                     gene2=data.X2.smai[which(rownames(data.X2.smai)==b),])
sp <- ggplot(data.plot, aes(x=gene1, y=gene2)) +#xlim(0,4)+ylim(0,4)+
  geom_point(size=1.5, alpha=0.6) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14)) +
  ylab(b) + xlab(a) + xlim(0,8)+ylim(0,8)+
  stat_smooth(method = "lm",
              formula = y ~ 0+x,
              geom = "smooth")+
  annotate("text",x=1, y=8,label=(paste0("slope==",coef(lm(data.plot$gene2~0+data.plot$gene1))[1])),parse=TRUE)
sp







