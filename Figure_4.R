# Figure 4
# Author: Yue Xing (yxing4@luc.edu)

	library(ggplot2)
	library(ggpubr)

	d1=read.csv("health_all.metaphlan3_species.rm.ra.rm_low.heatmap.csv",row.names=1)

	ut=read.csv("hea_heatmap_bacteria_count.thres_0.data_trans_human.k_2.cluster.csv",row.names=1)

	u1=ut[ut[,2]==1,]
	u2=ut[ut[,2]==2,]

	ut1=d1[,colnames(d1)%in%u1[,1]]
	ut2=d1[,colnames(d1)%in%u2[,1]]

	#

	a1=data.frame(rowMeans(ut1))
	a1[,2]=a1[,1]
	a1=a1[order(a1[,1],decreasing=TRUE),]
	b1=a1[1:20,]

	ut1=ut1[rownames(b1),]
	rownames(ut1)=sub("[[:lower:]]*_",". ",rownames(ut1))
	ut1m=data.frame(rowMeans(ut1))
	ut1m$Taxon=rownames(ut1m)
	colnames(ut1m)[1]="Mean_RA"

	ut1m$SE=NA
	for (i in 1:nrow(ut1m)) {
		ut1m[i,"SE"]=sd(ut1[i,])/sqrt(ncol(ut1))
	}
	ut1m$Taxon=factor(ut1m$Taxon,levels = rev(ut1m$Taxon))

	p1 = ggplot(ut1m) +
	geom_bar(aes(x=Taxon, y=Mean_RA), stat="identity", fill="skyblue", alpha=0.7) +
	coord_flip() +
	ggtitle("UT1") +
	geom_errorbar( aes(x=Taxon, ymin=Mean_RA-SE, ymax=Mean_RA+SE), width=0.4, colour="orange", alpha=0.9, size=1)

	#

	a1=data.frame(rowMeans(ut2))
	a1[,2]=a1[,1]
	a1=a1[order(a1[,1],decreasing=TRUE),]
	b1=a1[1:20,]

	ut2=ut2[rownames(b1),]
	rownames(ut2)=sub("[[:lower:]]*_",". ",rownames(ut2))
	ut2m=data.frame(rowMeans(ut2))
	ut2m$Taxon=rownames(ut2m)
	colnames(ut2m)[1]="Mean_RA"

	ut2m$SE=NA
	for (i in 1:nrow(ut2m)) {
		ut2m[i,"SE"]=sd(ut2[i,])/sqrt(ncol(ut2))
	}
	ut2m$Taxon=factor(ut2m$Taxon,levels = rev(ut2m$Taxon))

	p2 = ggplot(ut2m) +
	geom_bar(aes(x=Taxon, y=Mean_RA), stat="identity", fill="skyblue", alpha=0.7) +
	coord_flip() +
	ggtitle("UT2") +
	geom_errorbar( aes(x=Taxon, ymin=Mean_RA-SE, ymax=Mean_RA+SE), width=0.4, colour="orange", alpha=0.9, size=1)

	#

	png("E://IUMP//clr_080521//manuscript_v4//UTs.top20.v3.png",res=300,width=3000,height=3000)
	ggarrange(p1+ rremove("xlab"), p2, 
		labels = c("a)", "b)"),
		nrow = 2)

	dev.off()
