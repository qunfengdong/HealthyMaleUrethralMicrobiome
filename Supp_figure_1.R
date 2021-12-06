# Graph correlation between Q-PCR and G. vag WGS reads
# Author: Yue Xing (yxing4@luc.edu)

# Get the mp3 counts
	d1=read.csv("health_all.metaphlan3_species.ra.csv")
	a=colnames(d1)
	a=a[-1]

	d1=read.csv("Gardnerella_vaginalis_mp3_counts.csv")
	colnames(d1)[1]="ID"
	d1=d1[d1[,1]%in%a,]
	d1[,1]=gsub("\\..*","",d1[,1])

# Merge with ut type data of the final heatmap
	ut=read.csv("hea_heatmap_bacteria_count.thres_0.data_trans_human.k_2.cluster.csv",row.names=1)
	d1=merge(d1,ut,by="ID")

# Merge with qpcr data
	library(readxl)
	qpcr=read_excel("IUMP_GV Bacterial Loads_12-15-20.xlsx", sheet = "Sheet1")
	qpcr=data.frame(qpcr)
	d2=merge(d1,qpcr,by.x="ID",by.y="Participant.ID",all.x=TRUE)
	d2$ratio=d2$Gardnerella_vaginalis/d2$Human_reads

# Merge with alr data

	d0=read.csv("hea_heatmap_bacteria_human.csv",row.names=1)
	d0=t(data.frame(d0["Gardnerella_vaginalis",]))
	colnames(d0)[1]="Gvags_ALR"

	rownames(d2)=d2[,1]
	d2=merge(d2,d0,by="row.names",all.x=TRUE)

# Corr for all samples
	for (mt in c("pearson", "kendall", "spearman")) {
		a=cor.test(d2$ratio,d2$G..vag.load.1.ng.DNA, method=mt, exact=FALSE)
		#print(a)
		print(c(mt,a$estimate,a$p.value))
	}

	for (mt in c("pearson", "kendall", "spearman")) {
			a=cor.test(d2$Gvags_ALR,d2$G..vag.load.1.ng.DNA, method=mt, exact=FALSE)
			#print(a)
			print(c(mt,a$estimate,a$p.value))
		}

	for (mt in c("pearson", "kendall", "spearman")) {
			a=cor.test(d2$Gvags_ALR,d2$ratio, method=mt, exact=FALSE)
			#print(a)
			print(c(mt,a$estimate,a$p.value))
		}
	
### Plot
	
	png("E://IUMP//clr_080521//manuscript_v4//qpcr_gvags_corr_f.png",
	res=300,width=2000,height=2000)

	plot(d2$Gvags_ALR,log(d2$G..vag.load.1.ng.DNA+1),xlab="ALR transformed G. vags counts",ylab="Log qPCR loads",pch=19,col="blue")
	abline(lm(log(d2$G..vag.load.1.ng.DNA+1) ~ d2$Gvags_ALR),lwd=2)
	text(-17,13.5,"rho = 0.82 (p = 1.2e-23)")

	dev.off()
