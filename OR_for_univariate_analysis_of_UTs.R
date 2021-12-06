# Odds ratios for univariate analysis of UTs
# Author: Yue Xing (yxing4@luc.edu)

# Make tables

	require(MASS)

	a=read.csv("hea_heatmap_bacteria_human.csv", row.names=1)
	a=colnames(a)

	sexM=read.csv("meta.csv")
	rownames(sexM)=sexM[,1]
	sexM=sexM[,-1]
	sexM=sexM[,1:9]

	covs=read.csv("meta_9202020.csv")
	rownames(covs)=covs[,1]
	covs=covs[,4:ncol(covs)]
	
	sexM=sexM[rownames(sexM)%in%a,]
	covs=covs[rownames(covs)%in%a,]
	
	dim(covs)

	cnt=read.csv("hea_heatmap_bacteria_count.thres_0.data_trans_human.k_2.cluster.csv",row.names=1)

	rownames(cnt)=cnt[,1]
		
	for (i in 1:ncol(sexM)) {
		a=data.frame(sexM[,i])
		rownames(a)=rownames(sexM)
		colnames(a)=colnames(sexM)[i]
		print(colnames(a))
		
		dt=merge(cnt,a,by="row.names")
		rownames(dt)=dt[,1]
		dt=dt[,-1]
		write.csv(dt,paste0("urethraltype_logistic//hea_heatmap_bacteria_human.",colnames(a),".nocov.input.csv"),row.names=FALSE,quote=FALSE)
		
		
		dt2=merge(dt,covs,by="row.names")
		rownames(dt2)=dt2[,1]
		dt2=dt2[,-1]
		write.csv(dt2,paste0("urethraltype_logistic//hea_heatmap_bacteria_human.",colnames(a),".wcov.input.csv"),row.names=FALSE,quote=FALSE)
	}
	
# Calculate ORs and p-values
## Confusion matrix with correction

	source("scripts//calc_OR_urethraltype_universal_cof_matrix_correction_nocov.R")

	all_file=list.files(".","hea_heatmap_bacteria_human.*.nocov.input.csv")

	taball=data.frame()
	for (fl in all_file) {
		print(fl)
		tab0=get_stat(fl)
		taball=rbind(taball,tab0)
	}

	write.csv(taball,"urethraltype_nocov_conf_matrix.csv",quote=FALSE,row.names=FALSE)

## logistic regression with covariates

	source("scripts//calc_OR_one_y_add_covariate_ps.R")

	all_file=list.files(".","hea_heatmap_bacteria_human.*.wcov.input.csv")

	taball_tt=taball_e=data.frame()
	for (fl in all_file) {
		print(fl)
		tab0=calc_OR_cov(fl)
		taball_tt=rbind(taball_tt,tab0$tt)
		taball_e=rbind(taball_e,tab0$tabe)
	}

	taball_tt=taball_tt[-which(is.na(taball_tt[,1])),]
	taball_e=taball_e[-which(is.na(taball_e[,1])),]

	write.csv(taball_tt,"urethraltype_wcov_logistic.csv",quote=FALSE,row.names=FALSE)
	write.csv(taball_e,"urethraltype_wcov_logistic_cov_p.csv",quote=FALSE)

