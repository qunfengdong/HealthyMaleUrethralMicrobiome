# OR on taxa associated with sex behavior
# Author: Yue Xing (yxing4@luc.edu)

# Make input data

	cnt=read.csv("health_all.metaphlan3_species.ra.csv",row.names=1)
	cnt=t(cnt)

	sexM=read.csv("meta.csv")
	rownames(sexM)=sexM[,1]
	sexM=sexM[,-1]

	covs=read.csv("meta_9202020.csv")
	rownames(covs)=covs[,1]
	covs=covs[,4:ncol(covs)]
	
	sexM=sexM[rownames(sexM)%in%rownames(cnt),]
	covs=covs[rownames(covs)%in%rownames(cnt),]

	for (i in 1:ncol(sexM)) {
		a=data.frame(sexM[,i])
		rownames(a)=rownames(sexM)
		colnames(a)=colnames(sexM)[i]
		print(colnames(a))
		
		dt=merge(cnt,a,by="row.names")
		write.csv(dt,paste0("health_all.metaphlan3_species.ra.",colnames(a),".nocov.input.csv"),row.names=FALSE,quote=FALSE)
		
		rownames(dt)=dt[,1]
		dt=dt[,-1]
		dt2=merge(dt,covs,by="row.names")
		write.csv(dt2,paste0("health_all.metaphlan3_species.ra.",colnames(a),".wcov.input.csv"),row.names=FALSE,quote=FALSE)
	}

# Calculate OR and p-values
## Without covariates using confusion matrix with correction

	require(MASS)
	library(epiR)

	source("~/scripts/calc_OR_taxa_universal_cof_matrix_correction_nocov.R")

	all_file=list.files(".","health_all.metaphlan3_species.ra.*.nocov.input.csv")

	taball=data.frame()
	for (fl in all_file) {
		print(fl)
		for (cutoff in c(0)) {
			tab0=calc_OR_conf_matrix(fl,cutoff)
			taball=rbind(taball,tab0)
		}
	}

	write.csv(taball,"health_all.metaphlan3_species.ra.nocov.conf_table.corrected.OR.csv",quote=FALSE,row.names=FALSE)

## With covariates using logistic regression

	source("scripts//calc_OR_taxa_universal_cov.miss.levels.R")

	all_file=list.files(".","health_all.metaphlan3_species.ra.*.wcov.input.csv")

	taball_tt=taball_e=data.frame()
	for (fl in all_file) {
		print(fl)
		for (cutoff in c(0)) {
			tab0=calc_OR_cov(fl,cutoff)
			taball_tt=rbind(taball_tt,tab0$tt)
			taball_e=rbind(taball_e,tab0$tabe)
		}
	}

	which(is.na(taball_tt[,1]))
	which(is.na(taball_e[,1]))
	taball_tt=taball_tt[-which(is.na(taball_tt[,1])),]
	taball_e=taball_e[-which(is.na(taball_e[,1])),]

	write.csv(taball_tt,"health_all.metaphlan3_species.ra.wcov.OR.csv",quote=FALSE,row.names=FALSE)
	write.csv(taball_e,"health_all.metaphlan3_species.ra.wcov.cov_p.csv",quote=FALSE)

## Check if covariates are significant

	source("scripts//calc_OR_taxa_universal_correction_wcov.R")
	
	all_file=list.files(".","health_all.metaphlan3_species.ra.*.wcov.input.csv")
	
	taball_tt=data.frame()
	for (fl in all_file) {
		print(fl)
		for (cutoff in c(0)) {
			tab0=calc_OR_test_cov(fl,cutoff)
			taball_tt=rbind(taball_tt,tab0)
		}
	}
	
	taball_tt=taball_tt[-which(is.na(taball_tt[,1])),]
	
	colnames(taball_tt)=c("Intercept","Cutoff","a_Yes","Age","b_Yes","c_Yes",
		"d_Yes","e_Yes","Race_Asian","Race_Black_or_African_American",
		"Race_More_than_one_race","Race_White","STD_life_Yes",
		"Taxa","Type")
	
	# Adjust for multi-test corrections
	tabb=character()
	for (a in names(table(taball_tt$Type))) {
		tmp=taball_tt[taball_tt$Type==a,]
		for (i in 3:13) {
			tmp[,i]=p.adjust(tmp[,i], method="BH")
		}
		tabb=rbind(tabb,tmp)
	}
	
	write.csv(tabb,"health_all.metaphlan3_species.ra.wcov.cov_p.adj.csv",row.names=FALSE,quote=FALSE)

# Find significant ORs for pairwise comparisons
## Without covariates

	fl="health_all.metaphlan3_species.ra.nocov.conf_table.corrected.OR.csv"
	d1=read.csv(fl)
	idx=numeric()
	for (i in 1:nrow(d1)) {
		if (!is.na(d1[i,"OR.strata.wald.CI2.5"]) & !is.na(d1[i,"OR.strata.wald.CI97.5"]) & ((d1[i,"OR.strata.wald.CI2.5"]<1 & d1[i,"OR.strata.wald.CI97.5"]<1) | (d1[i,"OR.strata.wald.CI2.5"]>1 & d1[i,"OR.strata.wald.CI97.5"]>1))) {
			if ((d1[i,"OR.strata.wald"] >= d1[i,"OR.strata.wald.CI2.5"]) & (d1[i,"OR.strata.wald"] <= d1[i,"OR.strata.wald.CI97.5"])) {
				#if (d1[i,"chi2.strata.fisher.p"]<0.05) {
					idx=c(idx,i)
				#} else {
				#	print(d1[i,1:5])
				#}

			}
		}
	}
	d2=d1[idx,]
	fl2=gsub("OR.csv","sigOR.csv",fl)
	write.csv(d2,fl2,quote=FALSE,row.names=FALSE)

## With covariates

	fls=list.files(".","wcov.OR.csv")

	for (fl in fls) {
		d1=read.csv(fl)
		idx=numeric()
		for (i in 1:nrow(d1)) {
			if (!is.na(d1[i,"CI_2.5"]) & !is.na(d1[i,"CI_97.5"]) & ((d1[i,"CI_2.5"]<1 & d1[i,"CI_97.5"]<1) | (d1[i,"CI_2.5"]>1 & d1[i,"CI_97.5"]>1))) {
				if ((d1[i,"OR"] >= d1[i,"CI_2.5"]) & (d1[i,"OR"] <= d1[i,"CI_97.5"])) {
					#if (d1[i,"p_val"]<0.05) {
						idx=c(idx,i)
					#}

				}
			}
		}
		d2=d1[idx,]
		fl2=gsub("OR.csv","sigOR.csv",fl)
		write.csv(d2,fl2,quote=FALSE,row.names=FALSE)
	}

