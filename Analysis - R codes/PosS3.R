#######################################################################
###############  data analysis: Task PosS2
###############     July 18, 2023
#######################################################################


library(caret)
library(Seurat)
library(uwot)
library(ggrepel)
library(ggplot2)
library(batchelor)
library(harmony)
library(rliger)
library(reticulate)
library(SMAI)
reticulate::use_python("/Users/rongma/opt/miniconda3/envs/scanorama/bin/python")
scanorama <- import('scanorama')


RNA = read.table('scRNA_count/scRNA_count.txt')
spatial = read.table('/Users/rongma/Dropbox/Eric-Rong-James/Alignment/scRNA_count/Spatial_count.txt',
                     header = T)

dim(RNA)
dim(spatial)

RNA.sub = t(RNA[match(colnames(spatial),rownames(RNA)),])
dim(RNA.sub)

k=42
temp=abs(cov(spatial))
diag(temp)=0
gene.sel=order(rowMeans(temp),decreasing=T)[1:(k/2)]
out1 = matrix(ncol=3, nrow=42)
out2 = matrix(ncol=3, nrow=42)
out3 = matrix(ncol=3, nrow=42)
out4 = matrix(ncol=3, nrow=42)
out5 = matrix(ncol=3, nrow=42)
out6 = matrix(ncol=3, nrow=42)
for(i in gene.sel){
  spatial.sub=spatial[,-i]

  # ############Seurat
  all.RNA = t(rbind(spatial.sub,
                    RNA.sub[,-i]))
  dim(all.RNA)
  meta.data= data.frame(method = c(rep("spatial", dim(spatial.sub)[1]),
                                   rep("RNA", dim(RNA.sub)[1])))
  row.names(meta.data)=colnames(all.RNA)
  all.data <- CreateSeuratObject(counts = all.RNA, project = "imputation", min.cells = 0, min.features = 0,
                                 meta.data = meta.data)
  all.data <- SplitObject(all.data, split.by = "method")
  all.data <- lapply(X = all.data, FUN = function(x) {
    x <- NormalizeData(x,verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  })
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = all.data,verbose = FALSE)
  data.X1 = as.matrix(all.data$spatial@assays$RNA@data)
  data.X2 = as.matrix(all.data$RNA@assays$RNA@data)
  data.anchors <- FindIntegrationAnchors(object.list = all.data, anchor.features = features,verbose = FALSE)
  data.combined <- IntegrateData(anchorset = data.anchors, verbose = FALSE)
  dim(data.combined@assays$integrated@data)
  data1=as.matrix(data.combined@assays$integrated@data[,1:dim(data.X1)[2]])
  data2=as.matrix(data.combined@assays$integrated@data[,(dim(data.X1)[2]+1):(dim(data.X1)[2]+dim(data.X2)[2])])

  knnmodel = knnreg(x=t(data2), y=RNA.sub[,i])
  pred_y = predict(knnmodel, newdata=data.frame(t(data1)))
  out2[i,1] = cor(pred_y, spatial[,i])
  out2[i,2] = cor(pred_y, spatial[,i], method="kendall")
  out2[i,3] = cor(pred_y, spatial[,i], method="spearman")
  #
  #
  # ############## SMAI
  mnn.out = findMutualNN(t(data1),t(data2), k1=50)
  set.seed(20)
  out.smai <- align(data.X1, data.X2, denoise = "scree", dir.map = "2to1", r.max=40, outlier.cut=1,
                   sel1=unique(mnn.out$first),  sel2=unique(mnn.out$second),
                    prop.align=0.6, test=FALSE, t=1)
  data1=out.smai$data1.integrate
  data2=out.smai$data2.integrate
  knnmodel = knnreg(t(data2), RNA.sub[,i])
  pred_y = predict(knnmodel, data.frame(t(data1)))
  out1[i,1] = cor(pred_y, spatial[,i])
  out1[i,2] = cor(pred_y, spatial[,i], method="kendall")
  out1[i,3] = cor(pred_y, spatial[,i], method="spearman")


  # ######################## Harmony

  all.RNA = t(rbind(spatial.sub,
                    RNA.sub[,-i]))
  dim(all.RNA)
  meta.data= data.frame(method = c(rep("spatial", dim(spatial.sub)[1]),
                                   rep("RNA", dim(RNA.sub)[1])))
  row.names(meta.data)=colnames(all.RNA)
  all.data <- CreateSeuratObject(counts = all.RNA, project = "spatial", min.cells = 0, min.features = 0,
                                 meta.data = meta.data)
  all.data  <- NormalizeData(all.data,verbose = FALSE)
  all.data  <- FindVariableFeatures(all.data , selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  all.data <- ScaleData(all.data, verbose = FALSE)
  all.data <- RunPCA(all.data, pc.genes = all.data@var.genes, npcs = 50, verbose = FALSE)
  all.data <- RunHarmony(all.data, "method", verbose = FALSE)
  dim(all.data@reductions$harmony@cell.embeddings)
  data1=t(as.matrix(all.data@reductions$harmony@cell.embeddings[1:dim(spatial.sub)[1],]))
  data2=t(as.matrix(all.data@reductions$harmony@cell.embeddings[(dim(spatial.sub)[1]+1):(dim(spatial.sub)[1]+dim(RNA.sub)[1]),]))

  knnmodel = knnreg(t(data2), RNA.sub[,i])
  pred_y = predict(knnmodel, data.frame(t(data1)))
  out3[i,1] = cor(pred_y, spatial[,i])
  out3[i,2] = cor(pred_y, spatial[,i], method="kendall")
  out3[i,3] = cor(pred_y, spatial[,i], method="spearman")


  # ################liger

  rm(data_liger)
  X = as.matrix(t(spatial.sub))
  colnames(X) = row.names(spatial)
  Y=t(RNA.sub[,-i])
  data_liger <- createLiger(raw.data=list(data1=X,
                                          data2=Y))
  data_liger <- rliger::normalize(data_liger)
  data_liger <- selectGenes(data_liger)
  data_liger <- scaleNotCenter(data_liger)
  data_liger <- optimizeALS(data_liger, k = 15)
  data_liger <- quantile_norm(data_liger)
  data1=t(data_liger@H$data1)
  data2=t(data_liger@H$data2)

  knnmodel = knnreg(t(data2), RNA.sub[match(colnames(data2), rownames(RNA.sub)),i])
  pred_y = predict(knnmodel, data.frame(t(data1)))
  out4[i,1] = cor(pred_y, spatial[match(colnames(data1), rownames(spatial)),i])
  out4[i,2] = cor(pred_y, spatial[match(colnames(data1), rownames(spatial)),i], method="kendall")
  out4[i,3] = cor(pred_y, spatial[match(colnames(data1), rownames(spatial)),i], method="spearman")

  ##############fastMNN

  out.fastmnn = fastMNN(data.X1,data.X2)
  data1=(as.matrix(assay(out.fastmnn)[,1:dim(spatial.sub)[1]]))
  data2=(as.matrix(assay(out.fastmnn)[,(dim(spatial.sub)[1]+1):(dim(spatial.sub)[1]+dim(RNA.sub)[1])]))
  knnmodel = knnreg(t(data2), RNA.sub[,i])
  pred_y = predict(knnmodel, data.frame(t(data1)))
  out5[i,1] = cor(pred_y, spatial[match(colnames(data1), rownames(spatial)),i])
  out5[i,2] = cor(pred_y, spatial[match(colnames(data1), rownames(spatial)),i], method="kendall")
  out5[i,3] = cor(pred_y, spatial[match(colnames(data1), rownames(spatial)),i], method="spearman")



  # ##############SCANORAMA
  datasets <- list(t(data.X1),t(data.X2))
  genes_list <- list((rownames(data.X1)),(rownames(data.X2)))
  library(reticulate)
  reticulate::use_python("/Users/rongma/opt/miniconda3/envs/scanorama/bin/python")
  scanorama <- import('scanorama')
  integrated.data <- scanorama$integrate(datasets, genes_list)
  corrected.data <- scanorama$correct(datasets, genes_list, return_dense=TRUE)
  integrated.corrected.data <- scanorama$correct(datasets, genes_list,
                                                 return_dimred=TRUE, return_dense=TRUE)
  data1=t(integrated.corrected.data[[2]][[1]])
  data2=t(integrated.corrected.data[[2]][[2]])

  knnmodel = knnreg(t(data2), RNA.sub[,i])
  pred_y = predict(knnmodel, data.frame(t(data1)))
  out6[i,1] = cor(pred_y, spatial[,i])
  out6[i,2] = cor(pred_y, spatial[,i], method="kendall")
  out6[i,3] = cor(pred_y, spatial[,i], method="spearman")

  print(i)
  print(apply(cbind(out1[,3],out2[,3],out3[,3],out4[,3],out5[,3],out6[,3]),2,median, na.rm=T))
  print(apply(cbind(out1[,3],out2[,3],out3[,3],out4[,3],out5[,3],out6[,3]),2,mean, na.rm=T))

}


temp=abs(cov(spatial))
diag(temp)=0
gene.sel=order(rowMeans(temp),decreasing=T)[1:(k/2)]
data.plot=data.frame(tau = (c(out1[gene.sel,2],out2[gene.sel,2],out3[gene.sel,2],out4[gene.sel,2],out5[gene.sel,2],out6[gene.sel,2])),
                     method = rep(c("SMAI", "Seurat","Harmony", "LIGER", "fastMNN", "Scanorama"),
                                  each = k/2))
p <- ggplot(data.plot, aes(x=reorder(method, tau, FUN=median, na.rm=T), y=tau)) +
  geom_boxplot() + ylab("kendall's tau") + xlab("method")+
  #geom_point(position = position_jitter(width = 0.2), color = "black") +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        legend.key.size = unit(1, 'cm'),
        text=element_text(size=20),
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=14))
pdf(file = "PosS3_kend_top50_cor.pdf", width=3,height=3)
p
dev.off()

########## alignment test
i=1
spatial.sub=spatial[,-i]

#preprocessing
all.RNA = t(rbind(spatial.sub,
                  RNA.sub[,-i]))
dim(all.RNA)
meta.data= data.frame(method = c(rep("spatial", dim(spatial.sub)[1]),
                                 rep("RNA", dim(RNA.sub)[1])))
row.names(meta.data)=colnames(all.RNA)
all.data <- CreateSeuratObject(counts = all.RNA, project = "imputation", min.cells = 0, min.features = 0,
                               meta.data = meta.data)
all.data <- SplitObject(all.data, split.by = "method")
all.data <- lapply(X = all.data, FUN = function(x) {
  x <- NormalizeData(x,verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = all.data,verbose = FALSE)
data.X1 = as.matrix(all.data$spatial@assays$RNA@data)
data.X2 = as.matrix(all.data$RNA@assays$RNA@data)
data.anchors <- FindIntegrationAnchors(object.list = all.data, anchor.features = features,verbose = FALSE)
data.combined <- IntegrateData(anchorset = data.anchors, verbose = FALSE)
dim(data.combined@assays$integrated@data)
data1=as.matrix(data.combined@assays$integrated@data[,1:dim(data.X1)[2]])
data2=as.matrix(data.combined@assays$integrated@data[,(dim(data.X1)[2]+1):(dim(data.X1)[2]+dim(data.X2)[2])])

#SMAI test for partial alignability
mnn.out = findMutualNN(t(data1),t(data2), k1=30)
set.seed(20)
out.smai <- align(data.X1, data.X2, denoise = "scree", dir.map = "2to1", r.max=40, outlier.cut=1,
                  sel1=unique(mnn.out$first),  sel2=unique(mnn.out$second),
                  prop.align=0.5, t=3, cutoff.t=2.5)
out.smai$p.value
#[1] 0.6359227

