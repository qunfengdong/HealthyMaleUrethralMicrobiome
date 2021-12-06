# Differential abundance and information for UTs for ALR transformed IUMP data
# Author: Yue Xing (yxing4@luc.edu)

	source("scripts//functions_for_differential_abundance.R")

# UT1 vs UT2
	
	d1=read.csv("hea_heatmap_bacteria_human.csv",row.names=1)
	cls=read.csv("hea_heatmap_bacteria_count.thres_0.data_trans_human.k_2.cluster.csv",row.names=1)
	d2=read.csv("health_all.metaphlan3_species.ra.csv",row.names=1)
	all(rownames(d2)==rownames(d1))
	all(colnames(d2)==colnames(d1))
	
	tab0=wilcoxon_add_info(d1,d2,cls)
	
	write.csv(tab0,"hea_heatmap_bacteria_human.diff_abundance.csv",quote=FALSE,row.names=FALSE)
	
# Sub UTs of UT1 and UT2

	d1=read.csv("hea_heatmap_bacteria_human.csv",row.names=1)
	cls=read.csv("hea_heatmap_bacteria_count.thres_0.data_trans_human.k_2.cluster.csv",row.names=1)
	d2=read.csv("health_all.metaphlan3_species.ra.csv",row.names=1)
	all(rownames(d2)==rownames(d1))
	all(colnames(d2)==colnames(d1))

	cls1=cls[cls[,2]==1,]
	d11=d1[,colnames(d1)%in%cls1[,1]]
	d21=d2[,colnames(d2)%in%cls1[,1]]
	all(rownames(d21)==rownames(d11))
	all(colnames(d21)==colnames(d11))
	gp1=read.csv("hea_heatmap_bacteria_human_ut1.thres_NA.data_trans_none.k_2.cluster.csv",row.names=1)

	cls2=cls[cls[,2]==2,]
	d12=d1[,colnames(d1)%in%cls2[,1]]
	d22=d2[,colnames(d2)%in%cls2[,1]]
	all(rownames(d22)==rownames(d12))
	all(colnames(d22)==colnames(d12))
	gp2=read.csv("hea_heatmap_bacteria_human_ut2.thres_NA.data_trans_none.k_2.cluster.csv",row.names=1)

	tab1=wilcoxon_add_info(d11,d21,gp1)
	tab2=wilcoxon_add_info(d12,d22,gp2)
	write.csv(tab1,"hea_heatmap_bacteria_human.UT1.diff_abundance.csv",quote=FALSE,row.names=FALSE)
	write.csv(tab2,"hea_heatmap_bacteria_human.UT2.diff_abundance.csv",quote=FALSE,row.names=FALSE)
	
