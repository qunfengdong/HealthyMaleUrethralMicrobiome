# Author: Yue Xing (yxing4@luc.edu)

get_stat=function(d0) {

	require(epiR)
	tmpi=read.csv(d0,row.names=1)

	n11i=nrow(tmpi[tmpi[,1]==1 & tmpi[,2]==1,])
	n10i=nrow(tmpi[tmpi[,1]==1 & tmpi[,2]==0,])
	n01i=nrow(tmpi[tmpi[,1]==2 & tmpi[,2]==1,])
	n00i=nrow(tmpi[tmpi[,1]==2 & tmpi[,2]==0,])

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

	d0=gsub("hea_heatmap_bacteria_human.","",d0)
	d0=gsub(".nocov.input.csv","",d0)

	outf=cbind(d0,round(out$OR.strata.wald,4),round(out$OR.strata.score,4),
		round(out$chi2.strata.uncor[,4],4),round(out$chi2.strata.yates[,4],4),
		round(out$chi2.strata.fisher[,4],4),info)

	colnames(outf)=c("Behavior","OR.strata.wald","OR.strata.wald.CI2.5","OR.strata.wald.CI97.5",
		"OR.strata.score","OR.strata.score.CI2.5","OR.strata.score.CI97.5",
		"chi2.strata.uncor.p","chi2.strata.yates.p","chi2.strata.fisher.p",
		"taxa_behav")

	return(outf)
}

