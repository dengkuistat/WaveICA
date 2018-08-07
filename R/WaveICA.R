#' @title WaveICA
#' @description Removing batch effects for metabolomics data.
#' @author Kui Deng
#' \email{dengkui_stat@163.com}
#' @param data Sample-by-matrix metabolomics data.
#' @param wf selecting wavelet functions, the default is "haar".
#' @param batch batch labels.
#' @param group denoting the biological group such as disease vs group.
#' This param is optional. The default is NULL.
#' @param K The maximal component that ICA decomposes.
#' @param t The threshold to consider a component associate with the batch,
#' should be between 0 and 1.
#' @param t2 The threshold to consider a component associate with the group,
#' should be between 0 and 1.
#' @param alpha The trade-off value between the independence of samples and those
#' of variables and should be between 0 and 1.
#' @return A list that contains the clean data.
#' @export
#' @examples
#' \dontrun{
#' ## load the demo data
#'data(Amide_data, package = "WaveICA")
#'
#'##create a folder for demo
#'dir.create("Demo for WaveICA")
#'setwd("Demo for WaveICA")
#'### Data management for original data
#' data_Amide_stat_merge_order<-Amide_data[order(Amide_data$injection.order),]
#' data.Amide.stat_order<-data_Amide_stat_merge_order[,-c(1:5)]
#' ### Generating group and batch
#' data_Amide_stat_merge_order$group[data_Amide_stat_merge_order$group=="QC"]<-2
#' group_zong_Amide<-as.numeric(data_Amide_stat_merge_order$group)
#' batch_zong_Amide<-data_Amide_stat_merge_order$batch
#' group_sample_Amide<-group_zong_Amide[group_zong_Amide!=2]
#' batch_sample_Amide<-batch_zong_Amide[group_zong_Amide!=2]
#' batch_qc_Amide<-batch_zong_Amide[group_zong_Amide==2]
#'
#' ##### Separation of QC samples and subject samples
#'data.Amide.stat.order.sample<-data.Amide.stat_order[group_zong_Amide!=2,]
#'data.Amide.stat.order.QC<-data.Amide.stat_order[group_zong_Amide==2,]
#'
#' ################# Performances of Original data
#'#### PCA score plots
#'bitmap("PCA_original(group).tiff",type="jpeg",res=500)
#'pca(data=data.Amide.stat_order,id=group_zong_Amide+1,plot=2,label=c("CE","CRC","QC"))
#'dev.off()
#'bitmap("PCA_original(batch).tiff",type="jpeg",res=500)
#'pca(data=data.Amide.stat_order,id=batch_zong_Amide,plot=2,label=c("Batch1","Batch2","Batch3","Batch4"))
#'dev.off()
#'#### Dustances between QCS
#'pc.original<-prcomp(data.Amide.stat_order,scale.=T)
#'pc.original_select<-pc.original$x[group_zong_Amide==2,1:3]
#'dist.original<-dist(pc.original_select,method ="euclidean")
#'dist.original<-as.matrix(dist.original)
#'sum(dist.original)/(85*85-85)
#'##### Distances between subject samples
#'pc.original<-prcomp(data.Amide.stat_order,scale.=T)
#'pc.original_select<-pc.original$x[group_zong_Amide!=2,1:3]
#'dist.original<-dist(pc.original_select,method ="euclidean")
#'dist.original<-as.matrix(dist.original)
#'sum(dist.original)/(644*644-644)
#'
#'## Variable selection
#'univariate_original_Amide<-var_select(data=data.Amide.stat.order.sample,label=group_sample_Amide,t.test=T,Wilcox=F,AUC=F,FDR=T,VIP=T,FC=T,comps=3)
#'var_select_original_Amide<-rownames(univariate_original_Amide[univariate_original_Amide$VIP>1&univariate_original_Amide$P.FDR<0.05,])
#'
#'
#'## Predictive accuracy
#'svm.model.Amide<-SVM_MODEL(X=data.Amide.stat.order.sample[,var_select_original_Amide],Y=group_sample_Amide,kernel="radial",kfold=5)
#'####  Predictive accuracy with the same number of variables
#'univariate_original_Amide_order<-univariate_original_Amide[order(univariate_original_Amide$VIP,decreasing=T),]
#'auc_zong_original_Amide<-c()
#'for (i in seq(50,1000,50)){
#'  cat(paste("###########################",i,"variable####################\n"))
#'  var_select_original<-rownames(univariate_original_Amide_order)[1:i]
#'  svm.model.original<-SVM_MODEL(X=data.Amide.stat.order.sample[,var_select_original],Y=group_sample_Amide,kernel="radial",kfold=5)
#'  AUC<-svm.model.original$auc
#'  auc_zong_original_Amide<-c(auc_zong_original_Amide,AUC)
#'}

#'#### Heatmap of correlation in QCS
#'library(ggfortify)
#'library(RColorBrewer)
#'myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
#'bitmap(file="heatmap_Original.tiff",type="jpeg",res=500)
#'rownames(data.Amide.stat.order.QC)<-paste("QC",1:85,sep="")
#'cor_original<-cor(t(data.Amide.stat.order.QC))
#'autoplot(cor_original,geom="tile")+scale_fill_gradientn(colors=myPalette(4),limits=c(0.5,1),name="r value")+labs(x="",y="")+
#'  theme(axis.text.x=element_text(angle=90))+
#'  scale_x_discrete(limits=paste("QC",c(1,c(1:17)*5),sep=""))+
#'  scale_y_discrete(limits=paste("QC",c(1,c(1:17)*5),sep=""))+
#'  theme(axis.text=element_text(face="bold"))+
#' guides(fill=F)
#'dev.off()
#'## Scatterplot Matrix
#'bitmap(file="QC_scatterplotmatrix_Original.tiff",type="jpeg",res=500)
#'pairs(~QC21+QC25+QC34+QC46+QC62+QC68,data=log(t(data.Amide.stat.order.QC)),pch=16,gap=0.5,font.labels = 2,cex=0.8)
#'dev.off()
#'############### run WaveICA
#'data_wave_reconstruct_Amide<-WaveICA(data=data.Amide.stat_order,wf="haar",batch=batch_zong_Amide,group=group_zong_Amide,K=20,t=0.05,t2=0.05,alpha=0)
#'data_Amide_zong_wave<-data_wave_reconstruct_Amide$data_wave

#'data_Amide_sample_wave<-data_Amide_zong_wave[group_zong_Amide!=2,]
#'data_Amide_qc_wave<-data_Amide_zong_wave[group_zong_Amide==2,]

#' ######## Performances of data after removing batch effects by WaveICA
#'### PCA score plot
#'bitmap(file="PCA_WaveICA(group).tiff",type="jpeg",res=500)
#'pca(data=data_Amide_zong_wave,id=group_zong_Amide+1,plot=2,label=c("CE","CRC","QC"))
#'dev.off()
#'bitmap(file="PCA_WaveICA(batch).tiff",type="jpeg",res=500)
#'pca(data=data_Amide_zong_wave,id=batch_zong_Amide,plot=2,label=c("Batch1","Batch2","Batch3","Batch4"))
#'dev.off()
#'#### Distances between QCS
#'pc.WaveICA<-prcomp(data_Amide_zong_wave,scale.=T)
#'pc.WaveICA_select<-pc.WaveICA$x[group_zong_Amide==2,1:3]
#'dist.WaveICA<-dist(pc.WaveICA_select,method ="euclidean")
#'dist.WaveICA<-as.matrix(dist.WaveICA)
#'sum(dist.WaveICA)/(85*85-85)
#'
#'
#'#### Distances between subject samples
#'pc.WaveICA<-prcomp(data_Amide_zong_wave,scale.=T)
#'pc.WaveICA_select<-pc.WaveICA$x[group_zong_Amide!=2,1:3]

#'dist.WaveICA<-dist(pc.WaveICA_select,method ="euclidean")
#'dist.WaveICA<-as.matrix(dist.WaveICA)
#'sum(dist.WaveICA)/(644*644-644)
#'
#'
#'### Variable selection
#'univariate_wave_Amide<-var_select(data=data_Amide_sample_wave,label=group_sample_Amide,t.test=T,Wilcox=F,AUC=F,FDR=T,VIP=T,FC=T,comps=3)
#'
#'var_select_wave_Amide<-rownames(univariate_wave_Amide[univariate_wave_Amide$VIP>1&univariate_wave_Amide$P.FDR<0.05,])
#'
#'
#'
#'## Predictive accuracy
#'svm.model.wave.Amide<-SVM_MODEL(X=data_Amide_sample_wave[,var_select_wave_Amide],Y=group_sample_Amide,kernel="radial",kfold=5)
#'
#'
#'### Predictive accuracy with the same number of variables
#'univariate_wave_Amide_order<-univariate_wave_Amide[order(univariate_wave_Amide$VIP,decreasing=T),]
#'
#'auc_zong_wave_Amide<-c()
#'for (i in seq(50,1000,50)){
#'  cat(paste("###########################",i,"variable####################\n"))
#' var_select_wave<-rownames(univariate_wave_Amide_order)[1:i]
#'  svm.model.wave<-SVM_MODEL(X=data_Amide_sample_wave[,var_select_wave],Y=group_sample_Amide,kernel="radial",kfold=5)
#'  AUC<-svm.model.wave$auc
#'  auc_zong_wave_Amide<-c(auc_zong_wave_Amide,AUC)
#'}
#'
#'## Pearson correlation coefficients
#'corr_wave_Amide<-correlation(data.Amide.stat.order.QC,data_Amide_qc_wave,method="pearson")
#'mean(corr_wave_Amide$V3)
#'mean(corr_wave_Amide$V4)
#'
#'
#'library(ggfortify)
#'library(RColorBrewer)
#'myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
#'
#'
#'#### Heatmap of correlation in QCS
#'bitmap(file="Heatmap_WaveICA.tiff",type="jpeg",res=500)
#'
#'rownames(data_Amide_qc_wave)<-paste("QC",1:85,sep="")
#'cor_WaveICA<-cor(t(data_Amide_qc_wave))
#'autoplot(cor_WaveICA,geom="tile")+scale_fill_gradientn(colors=myPalette(4),limits=c(0.5,1),name="r value")+labs(x="",y="")+
#'  theme(axis.text.x=element_text(angle=90))+
#' scale_x_discrete(limits=paste("QC",c(1,c(1:17)*5),sep=""))+
#'  scale_y_discrete(limits=paste("QC",c(1,c(1:17)*5),sep=""))+
#'  theme(axis.text=element_text(face="bold"))+
#'  guides(fill=F)
#'dev.off()
#'
#'
#'## Scatterplot Matrix
#'bitmap(file="QC_scatterplotmatrix_WaveICA.tiff",type="jpeg",res=500)
#'pairs(~QC21+QC25+QC34+QC46+QC62+QC68,data=log(t(data_Amide_qc_wave)),pch=16,gap=0.5,font.labels = 2,cex=0.8)
#'dev.off()
#'
#'
#'###### Comparing the predictive accuracy
#'data_auc<-data.frame(auc=c(auc_zong_original_Amide,auc_zong_wave_Amide),var=c(seq(50,1000,50),seq(50,1000,50)),
#'                     group=c(rep(1,length(auc_zong_original_Amide)),rep(2,length(auc_zong_wave_Amide))))

#'library(ggplot2)
#'bitmap(file="Comparing AUC.tiff",type="jpeg",res=500)
#'ggplot(data=data_auc,aes(x=var,y=auc,col=factor(group)))+geom_line(aes(linetype=factor(group)),size=1)+
#'  scale_colour_brewer(palette="Set1",labels=c("Original","WaveICA","ComBat","QC-RLSC","ICA"))+
#'  scale_linetype_manual(values=c(2,1,4,3,5),labels=c("Original","WaveICA","ComBat","QC-RLSC","ICA"))+
#'  labs(x="# features selected",y="AUC")+
#'  theme(legend.title = element_blank(),axis.title = element_text(size=rel(1.2),face="bold"),
#'        legend.text=element_text(size=rel(1.2),face="bold"),legend.key.width=unit(1.5,"cm"))
#'dev.off()

#' }


WaveICA<-function(data,wf="haar",batch,group=NULL,K=20,t=0.05,t2=0.05,alpha=0){
  ### Wavelet Decomposition
  library(waveslim)
  level<-floor(log(nrow(data),2))
  if (is.null(colnames(data))){
    stop("data must have colnames")
  }
  coef<-list()
  for (k in 1:(level+1)){
    coef[[k]] <-matrix(NA,nrow(data),ncol(data))
  }
  for (j in 1:ncol(data)){
    cat(paste("######Decomposition",j,"########\n"))
    data_temp<-data[,j]
    x_modwt<-modwt(data_temp,wf=wf,n.levels =level)
    for (k in 1:(level+1)){
      coef[[k]][,j]<-x_modwt[[k]]
    }
  }
  ##### ICA
  index<-level+1
  data_wave_ICA<-list()
  for (i in (1:index)){
    cat(paste("######### ICA",i,"#############\n"))
    data_coef<-coef[[i]]
    data_coef_ICA<-normFact(fact="stICA",X=t(data_coef),ref=batch,refType ="categorical",k=K,t=t,ref2=group,refType2="categorical",t2=t2,alpha)
    data_wave_ICA[[i]]<-t(data_coef_ICA$Xn)
  }
  ### Wavelet Reconstruction
  index<-ncol(data)
  index1<-length(data_wave_ICA)
  data_coef<-matrix(NA,nrow(data_wave_ICA[[1]]),index1)
  data_wave<-matrix(NA,nrow(data_wave_ICA[[1]]),ncol(data_wave_ICA[[1]]))
  for (i in 1:index){
    cat(paste("######Reconstruction",i,"########\n"))
    for (j in 1:index1){
      data_coef[,j]<-data_wave_ICA[[j]][,i]
    }
    data_temp<-data[,i]
    data_coef<-as.data.frame(data_coef)
    colnames(data_coef)<-c(paste("d",1:(index1-1),sep=""),paste("s",(index1-1),sep=""))
    y<-as.list(data_coef)
    attributes(y)$class<-"modwt"
    attributes(y)$wavelet<-wf
    attributes(y)$boundary<-"periodic"
    data_wave[,i]<-imodwt(y)+mean(data_temp)
  }
  rownames(data_wave)<-rownames(data)
  colnames(data_wave)<-colnames(data)
  return(list(data_wave=data_wave))
}



