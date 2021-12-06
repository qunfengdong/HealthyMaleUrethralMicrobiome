# PCA for IUMP samples with normalziations
# Author: Yue Xing (yxing4@luc.edu)

	library(ggbiplot)
	
# ALR
	
	d1=read.csv("hea_heatmap_bacteria_human.csv",row.names=1)
	d1.pca <- prcomp(t(d1), center = TRUE,scale. = TRUE)
	summary(d1.pca)

	cls=read.csv("hea_heatmap_bacteria_count.thres_0.data_trans_human.k_2.cluster.csv",row.names=1)
	rownames(cls)=cls[,1]
	cls=cls[rownames(t(d1)),]

	all(rownames(cls)==rownames(t(d1)))

	png("hea_heatmap_bacteria_human_pca.png",width=6000,height=2500,res=300)
	ggbiplot(d1.pca,ellipse=TRUE,var.axes=FALSE,groups=as.factor(cls[,2]))+theme_bw()

	dev.off()

# ALR with virus

	d1=read.csv("hea_count.heatmap_plus_virus_human.csv",row.names=1)
	d1.pca <- prcomp(t(d1), center = TRUE,scale. = TRUE)
	summary(d1.pca)

	cls=read.csv("hea_heatmap_bacteria_count.thres_0.data_trans_human.k_2.cluster.csv",row.names=1)
	rownames(cls)=cls[,1]
	cls=cls[rownames(t(d1)),]

	all(rownames(cls)==rownames(t(d1)))

	png("hea_heatmap_bacteria_virus_human_pca.png",width=6000,height=2500,res=300)
	ggbiplot(d1.pca,ellipse=TRUE,var.axes=FALSE,groups=as.factor(cls[,2]))+theme_bw()

	dev.off()

# CLR

	d1=read.csv("hea_heatmap_bacteria_clr.csv",row.names=1)
	d1.pca <- prcomp(t(d1), center = TRUE,scale. = TRUE)
	summary(d1.pca)
	
	cls=read.csv("hea_heatmap_bacteria_count.thres_0.data_trans_clr.k_2.cluster.csv",row.names=1)
	rownames(cls)=cls[,1]
	cls=cls[rownames(t(d1)),]

	all(rownames(cls)==rownames(t(d1)))

	png("hea_heatmap_bacteria_clr_pca.png",width=6000,height=2500,res=300)
	ggbiplot(d1.pca,ellipse=TRUE,var.axes=FALSE,groups=as.factor(cls[,2]))+theme_bw()

	dev.off()

