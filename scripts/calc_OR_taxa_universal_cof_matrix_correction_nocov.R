# Author: Yue Xing (yxing4@luc.edu)

get_stat=function(testi,refi,tmpi) {

	n11i=nrow(tmpi[tmpi[,1]==1 & tmpi[,2]==testi,])
	n10i=nrow(tmpi[tmpi[,1]==1 & tmpi[,2]==refi,])
	n01i=nrow(tmpi[tmpi[,1]==0 & tmpi[,2]==testi,])
	n00i=nrow(tmpi[tmpi[,1]==0 & tmpi[,2]==refi,])

	ss=matrix(c(n11i,n10i,n01i,n00i),byrow=TRUE,ncol=2)
	ss=as.table(ss)

	if (0%in%ss) {
		ss=ss+0.5
	}

	info=data.frame(table(tmpi[,1],tmpi[,2]))
	info$inf=paste(info[,1],info[,2],sep="_")
	info=paste(info$inf,info$Freq,sep=":")
	info=paste(info,collapse=";")

	out=epi.2by2(ss)$massoc.detail
	outf=cbind(testi,refi,round(out$OR.strata.wald,4),round(out$OR.strata.score,4),
		round(out$chi2.strata.uncor[,4],4),round(out$chi2.strata.yates[,4],4),
		round(out$chi2.strata.fisher[,4],4),info)

	colnames(outf)=c("Test","Reference","OR.strata.wald","OR.strata.wald.CI2.5","OR.strata.wald.CI97.5",
		"OR.strata.score","OR.strata.score.CI2.5","OR.strata.score.CI97.5",
		"chi2.strata.uncor.p","chi2.strata.yates.p","chi2.strata.fisher.p",
		"taxa_behav")

	return(outf)
}

calc_OR_conf_matrix = function(nm,ct) {

	d1=read.csv(nm,row.names=1,na.strings=".")
	d2=numeric()

	if (ct==0) {
		# Y: 1 = bacteria, 0 = no bacteria
		for (i in 1:(ncol(d1)-1)) {
			line=d1[,i]
			for (j in 1:length(line)) {
				if (line[j]!=0) {line[j]=1}
			}
			d2=cbind(d2,line)
		}
	} else if (ct==0.1) {
		# Y: 1 = bacteria abundance > 0.1%, 0 = bacteria abundance <= 0.1%
		for (i in 1:(ncol(d1)-1)) {
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
		for (i in 1:(ncol(d1)-1)) {
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
		for (i in 1:(ncol(d1)-1)) {
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
		print("cut off error")
		quit(status=1)
	}

	d2=data.frame(d2)
	d2=cbind(d2,d1[,ncol(d1)])
	colnames(d2)=colnames(d1)
	rownames(d2)=rownames(d1)


	tab=data.frame()
	for (i in 1:(ncol(d2)-1)) {

		tmp=d2[,c(i,ncol(d2))]

		taxa=gsub(".*__","",colnames(tmp)[1])
		tp=colnames(tmp)[2]

		cc0=strsplit(tp,"_")[[1]]
	
		if (nchar(cc0[1])==1) {
			ref="0"
			test="1"

			outt=get_stat(test,ref,tmp)
			outt=cbind(taxa,tp,ct,outt)
			tab=rbind(tab,outt)

		} else if (nchar(cc0[1])==2) {
			cc1=strsplit(cc0[1],"")[[1]]
			lvs=c(paste0("n",cc1[1],cc1[2]),cc1[1],cc1[2],paste0("b",cc1[1],cc1[2]))

			test=lvs[2]
			ref=lvs[1]

			outt=get_stat(test,ref,tmp)
			outt=cbind(taxa,tp,ct,outt)
			tab=rbind(tab,outt)

			####

			test=lvs[3]
			ref=lvs[1]

			outt=get_stat(test,ref,tmp)
			outt=cbind(taxa,tp,ct,outt)
			tab=rbind(tab,outt)

			####

			test=lvs[4]
			ref=lvs[1]

			outt=get_stat(test,ref,tmp)
			outt=cbind(taxa,tp,ct,outt)
			tab=rbind(tab,outt)

			####

			test=lvs[3]
			ref=lvs[2]

			outt=get_stat(test,ref,tmp)
			outt=cbind(taxa,tp,ct,outt)
			tab=rbind(tab,outt)

			####

			test=lvs[4]
			ref=lvs[2]

			outt=get_stat(test,ref,tmp)
			outt=cbind(taxa,tp,ct,outt)
			tab=rbind(tab,outt)

			####

			test=lvs[4]
			ref=lvs[3]

			outt=get_stat(test,ref,tmp)
			outt=cbind(taxa,tp,ct,outt)
			tab=rbind(tab,outt)
		}
	}

	colnames(tab)[1:3]=c("Taxa","Type","Cutoff")
	return(tab)
}


