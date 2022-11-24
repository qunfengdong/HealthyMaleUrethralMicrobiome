# UTs vs metadata
# Author: Yue Xing (yxing4@luc.edu)

	d1=read.csv("meta_for_permanova.csv",na.strings=".",row.names=1)
	d2=read.csv("UTs.csv")
	d2=d2[!is.na(d2$UT),]
	colnames(d2)[1]="ID"
	
	d3=merge(d2,d1,all.x=TRUE,by.x="ID",by.y="row.names")
	tab=character()
	for (i in 3:ncol(d3)) {
		tmp=d3[,c(2,i)]
		if (length(table(tmp[,2]))>1) {
			ct2=chisq.test(tmp[,1],tmp[,2], simulate.p.value = TRUE)
			ct=chisq.test(tmp[,1],tmp[,2])
			tab=rbind(tab,c(colnames(tmp)[2],round(ct$p.value,4),round(ct2$p.value,4)))
		}
	}
	
	colnames(tab)=c("Name","p","simulated_p")
	tab=data.frame(tab)
	tab$p.adj=p.adjust(tab$p, method = "BH")
	tab$sim_p.adj=p.adjust(tab$simulated_p, method = "BH")
	
	write.csv(tab,"UT_vs_meta.csv",quote=FALSE,row.names=FALSE)