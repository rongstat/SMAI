
#' Spectral Manifold Alignment and Inference (SMAI)
#'
#' This function enables principled and interpretable alignability testing and structure-preserving integration of single-cell data.
#'
#' @param data1 data matrix (pxn1).
#' @param data2 data matrix (pxn2).
#' @param sel1 index set to select a subset of data1 for SMAI-align
#' @param sel2 index set to select a subset of data2 for SMAI-align
#' @param r.max intetger (default 200); the maximum number of singular values selected for desnoising and SMAI-alignment.
#' @param t integer (default 5); maximum number of iterations for SMAI-align.
#' @param denoise whether data should be denoised before searching for the spectral alignment; "scree" uses scree plot to denoise; "screeNOT" uses the ScreeNOT package to denoise; "none" does not perform denoising to the data.
#' @param test logical (default TRUE); whether SMAI-test is included in the workflow.
#' @param dir.map directionality of alignment; can be data1 aligned to data2 ("1to2"), data2 aligned to data1 ("2to1"), or internally determined acording to alignment error ("auto").
#' @param outlier.cut integer (default 20); the number of outliers to be removed from the reference data at each iteration.
#' @param knn integer (default 30); the number of nearest neighbors to be sampled from to do SMAI-test. A smaller value leads to more stringent selection of the test set.
#' @param prop.align numeric between 0 and 1 (default 0.3; and must not exceed 0.5 if test=TRUE, to ensure); proportion of samples used for SMAI-alignment. Smaller value improves speed with a possible cost of precision.
#' @param prop.inf numeric between 0 and 1-prop.align (default 0.5); proportion of samples used for SMAI-test. Smaller value improves speed with a possible cost of precision.
#' @param cutoff numeric larger than 1 (default 1.001); the threshold value for selection of number of informative singular values in SMAI-alignment. Effective only when denoise="scree."
#' @param cutoff.t numeric larger than 1 (default 1.5); the threshold value for data-driven selection of number of informative singular values in SMAI-test.
#' @return data1.integrate: integrated data matrices (pxn1);
#' @return data2.integrate: integrated data matrices (pxn2);
#' @return feature.rot: rotation matrix;
#' @return gamma: mean-shift vector;
#' @return beta: scaling factor;
#' @return p.value: p-value from alignability testing.
#' @export
align<-function(data1, data2, sel1=NULL, sel2=NULL, r.max=200, t=5, test=TRUE, dir.map = c("1to2","2to1","auto"),
                denoise = c("scree","screeNOT","none"), outlier.cut=20, knn=30,
                prop.align=0.3, prop.inf = 0.5, cutoff = 1.001, cutoff.t = 1.5){
  
  if(r.max>min(c(dim(data1), dim(data2)))){r.max = min(c(dim(data1), dim(data2)))}

  if(prop.align+prop.inf>1){ print("Error: prop.align + prop.inf > 1!") }else{
    
  n1 = dim(data1)[2]
  n2 = dim(data2)[2]

  ###denoising
  if(denoise == "scree"){
    svd.data1 = svds(data1, k=r.max)
    svd.data2 = svds(data2, k=r.max)
    if(sum(svd.data1$d[1:(length(svd.data1$d)-1)]/svd.data1$d[2:length(svd.data1$d)]>cutoff)==0){
      print("cutoff too large!")
    }
    if(sum(svd.data2$d[1:(length(svd.data2$d)-1)]/svd.data2$d[2:length(svd.data2$d)]>cutoff)==0){
      print("cutoff too large!")
    }
    r1=max(which(svd.data1$d[1:(length(svd.data1$d)-1)]/svd.data1$d[2:length(svd.data1$d)]>cutoff))
    r2=max(which(svd.data2$d[1:(length(svd.data2$d)-1)]/svd.data2$d[2:length(svd.data2$d)]>cutoff))
    data1.d = Tcrossprod(mat.mult(svd.data1$u[,1:r1],diag(svd.data1$d[1:r1])), svd.data1$v[,1:r1])
    data2.d = Tcrossprod(mat.mult(svd.data2$u[,1:r2],diag(svd.data2$d[1:r2])), svd.data2$v[,1:r2])
    r0 = min(c(r1,r2))
  }
  
  if(denoise == "none"){
    data1.d = data1
    data2.d = data2
    r0 = r.max
  }

  if(denoise == "screeNOT"){
    out1 = adaptiveHardThresholding(data1, r.max, strategy = "i")
    out2 = adaptiveHardThresholding(data2, r.max, strategy = "i")
    data1.d = out1$Xest
    data2.d = out2$Xest
    r0 = min(c(out1$r,out2$r))
  }


  ###select samples used for analysis
  if(is.null(sel1) & is.null(sel2)){
    id1.align = sample(n1, n1*prop.align)
    id2.align = sample(n2, n2*prop.align)
  }else{
    if(n1*prop.align>length(sel1)){print("prop.align too large!"); break}
    if(n2*prop.align>length(sel2)){print("prop.align too large!"); break}
    id1.align = sample(sel1, n1*prop.align)
    id2.align = sample(sel2, n2*prop.align)
  }



  ##############################
  ##   spectral alignment
  ##############################

  data1.align = data1.d[,id1.align]
  data2.align = data2.d[,id2.align]

  #sample matching
  if(length(id1.align)>=length(id2.align)){
    sub.id= sample(length(id1.align),length(id2.align))
    data.1.sub = data1.align[,sub.id]
    data.2.sub = data2.align
  }else{
    sub.id= sample(length(id2.align),length(id1.align))
    data.1.sub = data1.align
    data.2.sub = data2.align[,sub.id]
  }

  data.1.reg2 = data.1.sub-apply(data.1.sub, 1, mean)
  data.2.reg2 = data.2.sub-apply(data.2.sub, 1, mean)
  proc.dist = c()
  proc.dist.pair = c()
  M.set2 = 1:length(sub.id)
  M.set1 = 1:length(sub.id)
  for(tt in 1:t){

    if(tt>1){
      data.1.reg2 = Crossprod(rot.mat2, data.1.sub[,M.set1]-apply(data.1.sub[,M.set1], 1, mean))
      data.2.reg2 = Crossprod(t(rot.mat2), data.2.sub[,M.set2]-apply(data.2.sub[,M.set2], 1, mean))
    }

    #data.1.sub & data.2.sub: subset of data1.align or data2.align with matching sample sizes
    #data.1.reg2: centred data.1.sub[,M.set1] with feature space rotated
    prod.mat=Crossprod(data.1.reg2, data.2.sub[,M.set2]-apply(data.2.sub[,M.set2], 1, mean))

    #svd.prod.mat = svds(prod.mat, k=min(c(r0,dim(prod.mat))))
    svd.prod.mat = svds(prod.mat, k=min(c(r0, dim(prod.mat))))
    rot.mat = Tcrossprod(svd.prod.mat$u, svd.prod.mat$v)

    #data.1.reg: data.1.sub[,M.set1] with sample space aligned to data.2.sub[,M.set2]
    data.1.reg = mat.mult(data.1.sub[,M.set1], rot.mat)

    #data.2.reg: data.2.sub[,M.set2] with sample space aligned to data.1.sub[,M.set1]
    data.2.reg = mat.mult(data.2.sub[,M.set2], t(rot.mat))

    print(paste0("iteration ", tt,": sample space aligned!"))


    #feature space rotation
    prod.mat2= Tcrossprod(data.1.reg-apply(data.1.reg, 1, mean), data.2.sub[,M.set2]-apply(data.2.sub[,M.set2], 1, mean))
    svd.prod.mat2 = svds(prod.mat2, k=min(c(r0,dim(prod.mat2))))
    rot.mat2 = Tcrossprod(svd.prod.mat2$u,svd.prod.mat2$v)
    #data.1.reg3: centred data.1.sub[,M.set1] with sample and feature space rotated
    data.1.reg3 = Crossprod(rot.mat2, data.1.reg-apply(data.1.reg, 1, mean))

    #data.2.reg3: centred data.2.sub[,M.set2] with sample and feature space rotated
    data.2.reg3 = Crossprod(t(rot.mat2), data.2.reg-apply(data.2.reg, 1, mean))

    beta= sum(diag(Crossprod(data.1.reg3, data.2.sub[,M.set2]-apply(data.2.sub[,M.set2], 1, mean))))/sum(data.1.reg3^2)
    gamma = rowMeans(data.2.sub[,M.set2])
    Res = beta*data.1.reg3-(data.2.sub[,M.set2]-apply(data.2.sub[,M.set2], 1, mean))
    proc.dist[tt] = mean(Res^2)
    gamma.pair = rowMeans(data.1.sub[,M.set1])
    Res.pair = data.2.reg3-(data.1.sub[,M.set1]-apply(data.1.sub[,M.set1], 1, mean))*beta
    proc.dist.pair[tt] = mean(Res.pair^2)
    print(paste0("iteration ", tt,": feature space aligned!"))

    if(tt==1){
      if(((dir.map=="auto") & (proc.dist[tt]<proc.dist.pair[tt])) | (dir.map=="1to2")){
        data1.integrate = beta* Crossprod(rot.mat2, data1 - apply(data.1.sub[,M.set1], 1, mean)) + gamma %*% t(rep(1,n1))
        data2.integrate = data2
        feature.rot=  rot.mat2
        gamma.f = rowMeans(data.2.sub[,M.set2])
        beta.f = beta
        proc.dist.f=proc.dist
      }else{
        data1.integrate = data1
        data2.integrate = 1/beta* Crossprod(t(rot.mat2), data2 - apply(data.2.sub[,M.set2], 1, mean)) + gamma.pair %*% t(rep(1,n2))
        feature.rot=  t(rot.mat2)
        gamma.f = gamma.pair
        beta.f = 1/beta
        proc.dist.f=proc.dist.pair
      }
    }

    if(tt>1){

      if(((dir.map=="auto") & (proc.dist[tt]<proc.dist.pair[tt])) | (dir.map=="1to2")){
        if(proc.dist[tt]<proc.dist[tt-1]){
          data1.integrate = beta* Crossprod(rot.mat2, data1 - apply(data.1.sub[,M.set1], 1, mean)) + gamma %*% t(rep(1,n1))
          data2.integrate = data2
          feature.rot=  rot.mat2
          gamma.f = gamma
          beta.f = beta
          proc.dist.f=proc.dist
          bar2=mean(sort(colSums(Res^2),decreasing =F)[1:(floor(dim(Res)[2]/2))])*outlier.cut
          if(sum(colSums(Res^2)>bar2)>0){
            rm.set2 = which(colSums(Res^2)>bar2)
            M.set2 = M.set2[-rm.set2]
          }else{
            M.set2 = M.set2
          }
          bar1=mean(sort(colSums(Res.pair^2),decreasing =F)[1:(floor(dim(Res.pair)[2]/2))])*outlier.cut
          if(sum(colSums(Res.pair^2)>bar1)>0){
            rm.set1 = which(colSums(Res.pair^2)>bar1)
            M.set1 = M.set1[-rm.set1]
          }else{
            M.set1 = M.set1
          }

          if(tt==t){print("Data1 aligned to Data2!")}

        }else{

          print("Optimal alignment identified! Data1 aligned to Data2!")
          break
        }

      }

      if(((dir.map=="auto") & (proc.dist[tt]>proc.dist.pair[tt])) | (dir.map=="2to1")){
        if(proc.dist.pair[tt]<proc.dist.pair[tt-1]){
          data1.integrate = data1
          data2.integrate = 1/beta* Crossprod(t(rot.mat2), data2 - apply(data.2.sub[,M.set2], 1, mean)) + gamma.pair %*% t(rep(1,n2))
          feature.rot=  t(rot.mat2)
          gamma.f = gamma.pair
          beta.f = 1/beta
          proc.dist.f=proc.dist.pair
          bar2=mean(sort(colSums(Res^2),decreasing =F)[1:(floor(dim(Res)[2]/2))])*outlier.cut
          if(sum(colSums(Res^2)>bar2)>0){
            rm.set2 = which(colSums(Res^2)>bar2)
            M.set2 = M.set2[-rm.set2]
          }else{
            M.set2 = M.set2
          }
          bar1=mean(sort(colSums(Res.pair^2),decreasing =F)[1:(floor(dim(Res.pair)[2]/2))])*outlier.cut
          if(sum(colSums(Res.pair^2)>bar1)>0){
            rm.set1 = which(colSums(Res.pair^2)>bar1)
            M.set1 = M.set1[-rm.set1]
          }else{
            M.set1 = M.set1
          }

          if(tt==t){print("Data2 aligned to Data1!")}

        }else{
          print("Optimal alignment identified! Data2 aligned to Data1!")
          break
        }

      }

    }

  }

  print("Spectral alignment done!")

  print("Begin SMAI-test...")
  ##############################
  #######  SMAI-test
  ##############################

  ### preparation
  align.test <- function(data1, data2, K = NULL, cutoff.t = 1.5){


    p=dim(data1)[1]
    n=dim(data1)[2]

    lam1 = svd(data1)$d^2/n
    lam2 = svd(data2)$d^2/n


    #estimate k
    if(is.null(K)){
      if(min(lam1)<0.001){
        r0 = min(c(min(which(lam1<0.001))-1,30))
      }else{
        r0 = min(c(n,p))
      }

      if(sum(lam1[1:(r0-1)]/lam1[2:r0]>cutoff.t)>0){
        r1=max(which(lam1[1:(r0-1)]/lam1[2:r0]>cutoff.t))
      }else{
        r1=r0
        print("Warning: smaller cutoff.t recommended!")
      }

      if(min(lam2)<0.001){
        r0 = min(c(min(which(lam2<0.001))-1,30))
      }else{
        r0 = min(c(n,p))
      }
      if(sum(lam2[1:(r0-1)]/lam2[2:r0]>cutoff.t)>0){
        r2=max(which(lam2[1:(r0-1)]/lam2[2:r0]>cutoff.t))
      }else{
        r2=r0
        print("Warning: smaller cutoff.t recommended!")
      }
      K=min(c(r1,r2))
    }

    if(K==1){
      print("cutoff.t too smal! Used K=2.")
      K=2
    }

    #lam1 = lam1 * sum(lam2)/sum(lam1)
    lam1 = lam1 * sum(lam2[1:K])/sum(lam1[1:K])
    #estimate alpha_k
    m.lb <- function(z, lam, k, gamma){
      return(  -(1-gamma)/z+sum(1/(lam[-(1:k)]-z))/length(lam)    )
    }

    m.p <- function(z, lam, k, gamma){
      return(  -(1-gamma)/z^2+sum(1/(lam[-(1:k)]-z)^2)/length(lam)    )
    }

    alpha1 = c()
    alpha2 = c()
    psi1 = c()
    psi2 = c()
    for(i in 1:K){
      alpha1[i] = -1/m.lb(lam1[i], lam1, K, p/n)
      alpha2[i] = -1/m.lb(lam2[i], lam2, K, p/n)
      psi1[i] = 1/alpha1[i]^2/m.p(lam1[i], lam1, K, p/n)
      psi2[i] = 1/alpha2[i]^2/m.p(lam2[i], lam2, K, p/n)
    }


    #asymp covariance matrix
    asym.Sig1 = matrix(ncol=K,nrow=K)
    diag(asym.Sig1) = 2*alpha1^2*psi1

    asym.Sig2 = matrix(ncol=K,nrow=K)
    diag(asym.Sig2) = 2*alpha2^2*psi2

    t.stat = sqrt(n)*(lam1[1:K]-lam2[1:K])/sqrt(diag(asym.Sig1)+diag(asym.Sig2))
    p.value = 1-pchisq(sum(t.stat[1:K]^2),K)
    return(list(z.score = t.stat, p.value=p.value))
  }

  ####start SMAI-test

  if(test==TRUE){
    id1.inf = sample(setdiff(1:n1, id1.align), n1*prop.inf)
    id2.inf = sample(setdiff(1:n2, id2.align), n2*prop.inf)
    data1.inf = data1.d[,id1.inf]
    data2.inf = data2.d[,id2.inf]


    if(((dir.map=="auto") & (proc.dist[tt-1]<proc.dist.pair[tt-1])) | (dir.map=="1to2")){
      #find subset of data1.inf using data.1.sub[,M.set1]

      id.sel1=c()
      for(i in 1:length(M.set1)){
        #choose randomly from knn's, to reduce duplicates
        id.sel1[i]=sample(order(colSums((data.1.sub[,M.set1][,i]-data1.inf)^2),decreasing = F)[1:knn],1)
      }

      #find subset of data2.inf using f(data.1.sub[,M.set1])
      id.sel2=c()
      f.data.1.sub = beta* Crossprod(rot.mat2, data.1.sub[,M.set1] - apply(data.1.sub[,M.set1], 1, mean)) + gamma %*% t(rep(1,length(M.set1)))
      for(i in 1:length(M.set1)){
        #choose randomly from knn's, to reduce duplicates
        id.sel2[i]=sample(order(colSums((f.data.1.sub[,i]-data2.inf)^2),decreasing = F)[1:knn],1)
      }

      #alignability test
      data1.inf = data1[,id1.inf]
      data2.inf = data2[,id2.inf]


      p.value = align.test(t(data2.inf[,id.sel2]-apply(data2.inf[,id.sel2], 1, mean)),
                           t(data1.inf[,id.sel1]-apply(data1.inf[,id.sel1], 1, mean)),
                           cutoff.t = cutoff.t)$p.value
      print("Spectral inference done!")
    }else{

      id.sel2=c()
      for(i in 1:length(M.set2)){
        #choose randomly from knn's, to reduce duplicates
        id.sel2[i]=sample(order(colSums(abs(data.2.sub[,M.set2][,i]-data2.inf)^2),decreasing = F)[1:knn],1)
      }

      #find subset of data2.inf using f(data.1.sub[,M.set1])
      id.sel1=c()
      f.data.2.sub = 1/beta* Crossprod(t(rot.mat2), data.2.sub[,M.set2] - apply(data.2.sub[,M.set2], 1, mean)) + gamma.pair %*% t(rep(1,length(M.set2)))
      for(i in 1:length(M.set2)){
        #choose randomly from knn's, to reduce duplicates
        id.sel1[i]=sample(order(colSums(abs(f.data.2.sub[,i]-data1.inf)^2),decreasing = F)[1:knn],1)
      }

      #alignability test
      data1.inf = data1[,id1.inf]
      data2.inf = data2[,id2.inf]
      p.value = align.test(t(data2.inf[,id.sel2]-apply(data2.inf[,id.sel2], 1, mean)),
                           t(data1.inf[,id.sel1]-apply(data1.inf[,id.sel1], 1, mean)),
                           cutoff.t = cutoff.t)$p.value
      print("Spectral inference done!")


    }
  }else{
    p.value=NA
  }


  if(proc.dist[tt]>proc.dist.pair[tt]){
    return(list(data1.integrate = data1.integrate,
                data2.integrate = data2.integrate,
                feature.rot=  feature.rot,
                gamma = rowMeans(data1)-rowMeans(data2),
                beta = beta.f,
                proc.dist=proc.dist.f,
                p.value=p.value))
  }else{
    return(list(data1.integrate = data1.integrate,
                data2.integrate = data2.integrate,
                feature.rot=  feature.rot,
                gamma = rowMeans(data2)-rowMeans(data1),
                beta = beta.f,
                proc.dist=proc.dist.f,
                p.value=p.value))
  }

}

}


#' Batch correction evaluation metrics
#'
#' This function include a few metrics for evaluating the batch effect removal.
#'
#' @param data integrated data matrix (pxn).
#' @param meta_data label information.
#' @param straggler labels of samples with non-overlapping class labels.
#' @param cutoff cutoff parameter for PCA dimension reduction.
#' @return sihouette.batch: 1 - mean sihouette index; ch.batch: inverse CH index; db.batch: DB index.
#' @export
assess <- function(data, meta_data, straggler=NULL, cutoff=1.05){

  if(is.null(straggler)){
    overlap = 1:dim(meta_data)[1]
  }else{
    overlap = (1:dim(meta_data)[1])[-straggler]
  }


  if(dim(data)[2]>=30){
    #denoised integrated data
    svd.data = svds(data, k=30)
    if(sum(svd.data$d[1:(length(svd.data$d)-1)]/svd.data$d[2:length(svd.data$d)]>cutoff)>0){
      r=max(which(svd.data$d[1:(length(svd.data$d)-1)]/svd.data$d[2:length(svd.data$d)]>cutoff))
    }else{
      r=30
    }
    data.int.d = mat.mult(svd.data$u[,1:r],diag(svd.data$d[1:r]))
  }else{
    data.int.d =data
  }

  # #mean silhouette index
  dis.mat2=dist(data.int.d[overlap,])
  sihouette.batch = mean(silhouette(as.numeric(meta_data$batch)[overlap], dist = dis.mat2)[,3])


  #CH index
  ch.batch=calinhara(data.int.d[overlap,], as.numeric(meta_data$batch)[overlap])

  #DB index
  db.batch=index.DB(data.int.d[overlap,], cl=as.numeric(meta_data$batch)[overlap])$DB

  return(list(sihouette.batch=1-sihouette.batch,
              ch.batch=1/ch.batch,
              db.batch=db.batch))
}


#' Structure-preserving metrics
#'
#' This function include three metrics for evaluating structure preservation after integration.
#'
#' @param data.registered integrated data matrix (nxp).
#' @param data original data matrix (nxp).
#' @return Kendall's tau ("kendall"), Spearman's rho ("spearman"), and Pearson ("pearson").
#' @export
faithful <- function(data.registered, data){

  n=dim(data)[1]
  sub.id = sample(n,min(c(1000,n)))

  dist.reg = as.matrix(dist((data.registered[sub.id,])))
  dist.ori = as.matrix(dist((data[sub.id,])))
  n=dim(dist.reg)[1]
  faithful=matrix(nrow=3, ncol=n)

  for(i in 1:n){
    faithful[1,i]=cor(dist.reg[i,], dist.ori[i,], method="kendall")
    faithful[2,i]=cor(dist.reg[i,], dist.ori[i,], method="spearman")
    faithful[3,i]=cor(dist.reg[i,], dist.ori[i,], method="pearson")
  }
  out = apply(faithful,1,mean)
  return(c(kendall=out[1], spearman=out[2], pearson=out[3]))
}


