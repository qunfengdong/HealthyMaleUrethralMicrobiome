# Author: Yue Xing (yxing4@luc.edu)

calc_OR_test_cov = function(nm,ct) {

	d1=read.csv(nm,row.names=1,na.strings=".")
	d2=numeric()
	
	if (ct==0) {
		# Y: 1 = bacteria, 0 = no bacteria
		for (i in 1:(ncol(d1)-10)) {
			line=d1[,i]
			for (j in 1:length(line)) {
				if (line[j]!=0) {line[j]=1}
			}
			d2=cbind(d2,line)
		}
	} else if (ct==0.1) {
		# Y: 1 = bacteria abundance > 0.1%, 0 = bacteria abundance <= 0.1%
		for (i in 1:(ncol(d1)-10)) {
			line=d1[,i]
			for (j in 1:length(line)) {
				if (line[j]>0.1) {
					line[j]=1
				} else {
					line[j]=0
				}
			}
			d2=cbind(d2,line)
		}
	}  else if (ct==0.5) {
		# Y: 1 = bacteria abundance > 0.5%, 0 = bacteria abundance <= 0.5%
		for (i in 1:(ncol(d1)-10)) {
			line=d1[,i]
			for (j in 1:length(line)) {
				if (line[j]>0.5) {
					line[j]=1
				} else {
					line[j]=0
				}
			}
			d2=cbind(d2,line)
		}
	} else if (ct==1) {
		# Y: 1 = bacteria abundance > 1%, 0 = bacteria abundance <= 1%
		for (i in 1:(ncol(d1)-10)) {
			line=d1[,i]
			for (j in 1:length(line)) {
				if (line[j]>1) {
					line[j]=1
				} else {
					line[j]=0
				}
			}
			d2=cbind(d2,line)
		}
	} else {
		print("nm error")
		quit(status=1)
	}

	d2=data.frame(d2)
	d2=cbind(d2,d1[,((ncol(d1)-9):ncol(d1))])
	colnames(d2)=colnames(d1)
	rownames(d2)=rownames(d1)

	d2$Race=factor(d2$Race)
	if ("Other" %in% names(table(d2$Race))) {
		d2$Race=relevel(d2$Race,ref="Other")
	}

	flag_all=c("3l_a","3l_b","3l_c","3l_d","3l_e",
		"360_a","360_b","360_c","360_d","360_e",
		"3l_abcde","360_abcde","3l_be","360_be",
		"3l_n","360_n")
	
	tab=data.frame(rep(NA,15))
	rownames(tab)=c("(Intercept)", "d2[, \"Age\"]", "d2[, \"Race\"]Asian", 
		"d2[, \"Race\"]Black or African American", "d2[, \"Race\"]More than one race", 
		"d2[, \"Race\"]White", "d2[, \"STD_life\"]Yes", "d2[, \"a\"]Yes", 
		"d2[, \"b\"]Yes", "d2[, \"c\"]Yes", "d2[, \"d\"]Yes", "d2[, \"e\"]Yes", 
		"taxa", "type", "cutoff")

	for (i in 1:(ncol(d2)-10)) {

		tmp=d2[,c(i,(ncol(d2)-9):ncol(d2))]

		taxa=gsub(".*__","",colnames(tmp)[1])
		tp=colnames(tmp)[2]

		y=factor(tmp[,1])
		y=relevel(y,ref="0")

		lmodel <- glm(y ~ d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"]+d2[,"a"]+d2[,"b"]+d2[,"c"]+d2[,"d"]+d2[,"e"], family = "binomial")

		dd=summary(lmodel)$coefficients
		dd=data.frame(dd[,"Pr(>|z|)"])
		dd[,1]=round(dd[,1],4)

		dd=rbind(dd,taxa,tp,ct)
		rownames(dd)[(nrow(dd)-2):nrow(dd)]=c("taxa","type","cutoff")

		tab=merge(tab,dd,by="row.names",all.x=TRUE)
		rownames(tab)=tab$Row.names
		tab=tab[,-1]
	}
	tab=t(tab)
	return(tab)
}


