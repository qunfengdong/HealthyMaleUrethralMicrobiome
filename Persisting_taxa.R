# Taxa that persist after reported vaginal sex exposure in two-time intervals
# Author: Yue Xing (yxing4@luc.edu)

# Prepare data

	cnt=read.csv("health_all.metaphlan3_species.ra.csv",row.names=1)
	cnt=t(cnt)

	covs=read.csv("meta_9202020.csv")
	rownames(covs)=covs[,1]
	covs=covs[,4:ncol(covs)]
	covs=covs[rownames(covs)%in%rownames(cnt),]

# 1) (V_60d- and (V_1yr+ or V_ever+)) vs (V_Ever-)

	sexM=read.csv("long1.csv")
	rownames(sexM)=sexM[,1]
	sexM=sexM[rownames(sexM)%in%rownames(cnt),]

	dt=merge(cnt,sexM,all.y=TRUE,by="row.names")
	print(dim(dt))
	write.csv(dt,paste0("long1.nocov.input.csv"),row.names=FALSE,quote=FALSE)
		
	rownames(dt)=dt[,1]
	dt=dt[,-1]
	dt2=merge(dt,covs,all.x=TRUE,by="row.names")
	print(dim(dt2))
	write.csv(dt2,paste0("long1.wcov.input.csv"),row.names=FALSE,quote=FALSE)

	# Remove the extra ID column at "DO". Change the sex behavior's column name to format in "V_60d".

# 2) (V_1yr- and V_ever+) vs (V_Ever-)

	sexM=read.csv("long2.csv")
	rownames(sexM)=sexM[,1]
	sexM=sexM[rownames(sexM)%in%rownames(cnt),]

	dt=merge(cnt,sexM,all.y=TRUE,by="row.names")
	print(dim(dt))
	write.csv(dt,paste0("long2.nocov.input.csv"),row.names=FALSE,quote=FALSE)
		
	rownames(dt)=dt[,1]
	dt=dt[,-1]
	dt2=merge(dt,covs,all.x=TRUE,by="row.names")
	print(dim(dt2))
	write.csv(dt2,paste0("long2.wcov.input.csv"),row.names=FALSE,quote=FALSE)

	#Remove the extra ID column at "DO". Change the sex behavior's column name to format in "V_60d".

# 3) V 60d-1yr+ vs V 1yr-

	sexM=read.csv("long3.csv")
	rownames(sexM)=sexM[,1]
	sexM=sexM[rownames(sexM)%in%rownames(cnt),]

	dt=merge(cnt,sexM,all.y=TRUE,by="row.names")
	print(dim(dt))
	write.csv(dt,paste0("long3.nocov.input.csv"),row.names=FALSE,quote=FALSE)
		
	rownames(dt)=dt[,1]
	dt=dt[,-1]
	dt2=merge(dt,covs,all.x=TRUE,by="row.names")
	print(dim(dt2))
	write.csv(dt2,paste0("long3.wcov.input.csv"),row.names=FALSE,quote=FALSE)

	#Remove the extra ID column at "DO". Change the sex behavior's column name to format in "V_60d".

# Calculate OR and p-values

	library(epiR)

	source("scripts//calc_OR_taxa_universal_cof_matrix_correction_nocov_stability.R")

	all_file=list.files(".","*.nocov.input.csv")

	for (fl in all_file) {
		print(fl)
		for (cutoff in c(0)) {
			tab0=calc_OR_conf_matrix(fl,cutoff)
			write.csv(tab0,gsub("input.csv","conf.OR.csv",fl),quote=FALSE,row.names=FALSE)
		}
	}

# Find significant ORs for pairwise comparisons

	fl="long1.nocov.conf.OR.csv"
	#fl="long2.nocov.conf.OR.csv"
	#fl="long3.nocov.conf.OR.csv"
	d1=read.csv(fl)
	idx=numeric()
	for (i in 1:nrow(d1)) {
		if (!is.na(d1[i,"OR.strata.wald.CI2.5"]) & !is.na(d1[i,"OR.strata.wald.CI97.5"]) & ((d1[i,"OR.strata.wald.CI2.5"]<1 & d1[i,"OR.strata.wald.CI97.5"]<1) | (d1[i,"OR.strata.wald.CI2.5"]>1 & d1[i,"OR.strata.wald.CI97.5"]>1))) {
			if ((d1[i,"OR.strata.wald"] >= d1[i,"OR.strata.wald.CI2.5"]) & (d1[i,"OR.strata.wald"] <= d1[i,"OR.strata.wald.CI97.5"])) {
				idx=c(idx,i)
			}
		}
	}
	d2=d1[idx,]
	fl2=gsub("OR.csv","sigOR.csv",fl)
	write.csv(d2,fl2,quote=FALSE,row.names=FALSE)


