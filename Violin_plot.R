# Violin plot
# Author: Yue Xing (yxing4@luc.edu)

library(ggplot2)

d1=read.csv("hea_heatmap_bacteria_human.csv",row.names=1)
d2=read.csv("Supplemental Table 5 Diff abundance for UT1 vs UT2 using ALR.csv")
ut=read.csv("UTs.csv")

d2=d2[d2$Wilcoxon<0.05,]
ta=d2[,1]
ta=gsub("\\.","",ta)
rownames(d1)=gsub("_"," ",rownames(d1))
d1=d1[ta,]
ut=ut[!is.na(ut[,2]),]
rownames(ut)=ut[,1]

for (i in 1:nrow(d1)) {
	tmp=t(d1[i,])
	tmp=merge(tmp,ut,all=TRUE,by="row.names")
	ta=colnames(tmp)[2]
	colnames(tmp)=c("ID","Abundance","ID2","UT")

	dp <- ggplot(tmp, aes(x=UT, y=Abundance, fill=UT)) + 
	  geom_violin(trim=FALSE)+
	  geom_boxplot(width=0.1, fill="white")+
	  labs(title=ta,x="UT", y = "ALR abundance") +
	  theme_classic() +
	  geom_jitter(shape=16, position=position_jitter(0.2))

	png(paste0(ta,"_violin_plot.png"),width=2000,height=2000,res=300)
	print(dp)
	dev.off()

}


