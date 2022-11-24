# Figure 4
# Author: Yue Xing (yxing4@luc.edu)

# 4a and 4b

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

	png("UTs.top20.png",res=300,width=3000,height=3000)
	ggarrange(p1+ rremove("xlab"), p2, 
		labels = c("a)", "b)"),
		nrow = 2)

	dev.off()

# 4c

	library(ggplot2)
	library(ggpubr)
	library(vegan)
	
	# Calculate alpha diversity
	nm="hea_heatmap_bacteria_UT1_ra"
	#nm="hea_heatmap_bacteria_UT2_ra"
	
	da <- read.csv(paste0(nm,'.csv'),row.names = 1)
	cntdata=t(da)
	richness <- t(estimateR(round(cntdata))[c("S.obs","S.chao1","S.ACE"),])
	shannon <- diversity(round(cntdata),index = "shannon")
	simpson <- diversity(round(cntdata),index = "simpson")
	evenness <- function(data) {diversity(data)/log(specnumber(data))}
	eveness <- evenness(round(cntdata))
	stats <- as.data.frame(list(richness,shannon,simpson,eveness))
	stats[is.na(stats)] <- 0
	colnames(stats) <- c("Obs","Chao1","ACE","Shannon","Simpson","Pielou")
	
	write.csv(stats,paste0(nm,".alpha_diversity.stats.csv"))

	# Plot alpha diversity
	ut2=read.csv("hea_heatmap_bacteria_UT1_ra.alpha_diversity.stats.csv",row.names=1)
	ut1=read.csv("hea_heatmap_bacteria_UT2_ra.alpha_diversity.stats.csv",row.names=1)

	ut2=data.frame(t(ut2))
	ut1=data.frame(t(ut1))

	ut1$M=rowMeans(ut1)
	ut1$SE=NA
	nc=(ncol(ut1)-2)
	for (i in 1:nrow(ut1)) {
		ut1[i,"SE"]=sd(ut1[i,1:nc])/sqrt(nc)
	}
	ut1=ut1[-nrow(ut1),]
	ut1$Name=rownames(ut1)

	ut2$M=rowMeans(ut2)
	ut2$SE=NA
	nc=(ncol(ut2)-2)
	for (i in 1:nrow(ut2)) {
		ut2[i,"SE"]=sd(ut2[i,1:nc])/sqrt(nc)
	}
	ut2=ut2[-nrow(ut2),]
	ut2$Name=rownames(ut2)

	ut1$UT="UT1"
	ut2$UT="UT2"

	ut1=ut1[,c("M","SE","Name","UT")]
	ut2=ut2[,c("M","SE","Name","UT")]

	uts=rbind(ut1,ut2)
	uts$Name[c(1,6)]=c("Observation","Observation")
	uts$Name=factor(uts$Name,levels = c("Observation","ACE","Chao1","Shannon","Simpson"))

	p1 = ggplot(uts, aes(x=Name, y=M, fill=UT, group=UT)) + 
		geom_bar(stat="identity", color="black", 
		position=position_dodge()) +
		geom_errorbar(aes(ymin=M-SE, ymax=M+SE), width=.4,
		size=1, colour="black", alpha=0.9, position = position_dodge(0.9)) +
		#coord_flip() +
		#ggtitle("MF") +
		labs(x = "Type", y="Alpha diversity") +
		theme_classic() +
		scale_fill_manual(values=c("deeppink1","steelblue1")) 

	png("a_diversity.png",width=3000,height=2000,res=300)
	plot(p1)
	dev.off()
