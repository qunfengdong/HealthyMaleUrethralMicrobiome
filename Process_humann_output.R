# Process humann output
# Author: Yue Xing (yxing4@luc.edu)

# Combine the same gene families of different tax, and Wilcox test on UT1 vs UT2

	d1=read.delim("healthy_genefamilies.ra.tsv")
	colnames(d1)=gsub("\\..*","",colnames(d1))
	d1[,1]=gsub("\\|.*","",d1[,1])

	library(dplyr)
	s=d1 %>%
	    group_by(X) %>%
	    summarise_all(funs(sum)) %>%
	    as.data.frame()
	
	write.table(s,"healthy_genefamilies.ra.combined.tsv",quote=FALSE,sep="\t",row.names=FALSE)

####

	d1=read.delim("healthy_genefamilies.ra.combined.tsv",row.names=1)

	ut=read.csv("UTs.csv")
	ut=ut[!(is.na(ut$UT)),]
	ut1=ut[ut$UT=="UT1",1]
	ut2=ut[ut$UT=="UT2",1]
	
	u1=d1[,colnames(d1)%in%ut1]
	u2=d1[,colnames(d1)%in%ut2]
	u1=as.matrix(u1)
	u2=as.matrix(u2)
	
	tab=character()
	for (i in 1:nrow(d1)) {
		tt=wilcox.test(u1[i,],u2[i,], paired = FALSE, alternative = "two.sided", exact=FALSE)
		if (!is.na(tt$p.value)) {
			if (tt$p.value<0.05) {
				tab=rbind(tab,c(rownames(d1)[i],tt$p.value))
			}
		}
	}
	write.table(tab,"healthy_genefamilies.ra.combined.ut1_vs_ut2.txt",quote=FALSE,sep="\t",row.names=FALSE)

# Extract the significantly diff abundant ones

	d1=read.delim("healthy_genefamilies.ra.combined.ut1_vs_ut2.txt")
	d1=d1[,1]
	
	d2=read.delim("healthy_genefamilies.ra.combined.tsv",row.names=1)
	colnames(d2)=gsub("\\..*","",colnames(d2))
	d2=d2[d1,]
	write.table(d2,"healthy_genefamilies.ra.combined.DA.tsv",quote=FALSE,sep="\t")

Add "# Gene Family	" at the beginning of healthy_genefamilies.ra.combined.DA.tsv.

## Map to GO terms

	humann_regroup_table -i healthy_genefamilies.ra.combined.DA.tsv -o healthy_genefamilies.ra.combined.DA.go.tsv -g uniref90_go -e 6
	humann_rename_table -i healthy_genefamilies.ra.combined.DA.go.tsv -o healthy_genefamilies.ra.combined.DA.go.go.tsv -n go

## Split to MF, CC, BP; find the top 20 abundant for plot

	s=read.delim("healthy_genefamilies.ra.combined.DA.go.go.tsv",row.names=1)

	s1=data.frame(rowSums(s))
	s1$Type=gsub(".*\\[","",rownames(s1))
	s1$Type=gsub("\\].*","",s1$Type)
	colnames(s1)[1]="Sum_RA"
	s1$Name=gsub(".*\\] ","",rownames(s1))
	s1=merge(s1,s,all.x=TRUE,by="row.names")
	
	s1mf=s1[s1$Type=="MF",]
	s1cc=s1[s1$Type=="CC",]
	s1bp=s1[s1$Type=="BP",]
	
	s1mf=s1mf[order(s1mf$Sum_RA,decreasing=TRUE),]
	s1cc=s1cc[order(s1cc$Sum_RA,decreasing=TRUE),]
	s1bp=s1bp[order(s1bp$Sum_RA,decreasing=TRUE),]
	
	write.table(s1mf,"healthy_genefamilies.ra.combined.DA.go.MF.tsv",quote=FALSE,sep="\t")
	write.table(s1cc,"healthy_genefamilies.ra.combined.DA.go.CC.tsv",quote=FALSE,sep="\t")
	write.table(s1bp,"healthy_genefamilies.ra.combined.DA.go.BP.tsv",quote=FALSE,sep="\t")
	
	rownames(s1mf)=s1mf[,4]
	s1mf=s1mf[,-c(1:4)]
	rownames(s1cc)=s1cc[,4]
	s1cc=s1cc[,-c(1:4)]
	rownames(s1bp)=s1bp[,4]
	s1bp=s1bp[,-c(1:4)]
	
	write.table(s1mf,"healthy_genefamilies.ra.combined.DA.go.MF.short.tsv",quote=FALSE,sep="\t")
	write.table(s1cc,"healthy_genefamilies.ra.combined.DA.go.CC.short.tsv",quote=FALSE,sep="\t")
	write.table(s1bp,"healthy_genefamilies.ra.combined.DA.go.BP.short.tsv",quote=FALSE,sep="\t")

####

	s1mf=s1mf[1:20,]
	s1cc=s1cc[1:20,]
	s1bp=s1bp[1:20,]

	trans_tab = function(dta) {
		tab=character()
		for (i in 1:nrow(dta)) {
			for (j in 1:ncol(dta)) {
				tab=rbind(tab,c(rownames(dta)[i],colnames(dta)[j],dta[i,j]))
			}
		}
		colnames(tab)=c("GO Term","Subject","RA")
		tab=data.frame(tab)
		tab$RA=as.numeric(as.character(tab$RA))
		#tab$log_RA=log(tab$RA+0.0001)
		#tab$Prev=ifelse(tab$RA==0,0,1)
		return(tab)
	}
	
	s2mf=trans_tab(s1mf)
	s2cc=trans_tab(s1cc)
	s2bp=trans_tab(s1bp)
	
	s2mf$GO.Term=factor(s2mf$GO.Term,levels = rev(rownames(s1mf)))
	s2cc$GO.Term=factor(s2cc$GO.Term,levels = rev(rownames(s1cc)))
	s2bp$GO.Term=factor(s2bp$GO.Term,levels = rev(rownames(s1bp)))

## Plot
	
	library(ggplot2)
	library(ggpubr)
	
	p1 <- ggplot(s2mf, aes(GO.Term,RA)) +
	geom_bar(position = 'dodge', stat = 'summary', fun = 'mean', fill="skyblue", alpha=0.7) +
	geom_errorbar(stat = 'summary', position = 'dodge', fun.data=mean_se, width=0.4, colour="orange", alpha=0.9, size=1) +
	#geom_point(aes(GO.Term,RA), shape = 19, position = position_dodge(width=1), col="yellow") +
	coord_flip() +
	ggtitle("MF") +
	labs(x = "Go terms", y="Relative abundance")
	
	p2 <- ggplot(s2cc, aes(GO.Term,RA)) +
	geom_bar(position = 'dodge', stat = 'summary', fun = 'mean', fill="skyblue", alpha=0.7) +
	geom_errorbar(stat = 'summary', position = 'dodge', fun.data=mean_se, width=0.4, colour="orange", alpha=0.9, size=1) +
	coord_flip() +
	ggtitle("CC") +
	labs(x = "Go terms", y="Relative abundance")
	
	p3 <- ggplot(s2bp, aes(GO.Term,RA)) +
	geom_bar(position = 'dodge', stat = 'summary', fun = 'mean', fill="skyblue", alpha=0.7) +
	geom_errorbar(stat = 'summary', position = 'dodge', fun.data=mean_se, width=0.4, colour="orange", alpha=0.9, size=1) +
	coord_flip() +
	ggtitle("BP") +
	labs(x = "Go terms", y="Relative abundance")
	
	png("healthy_genefamilies.ra.combined.DA.go.top20.png",width=10000,height=2000,res=300)
	plot(ggarrange(p3, p2+rremove("ylab"), p1+rremove("ylab"),
	      labels = c("a)", "b)", "c)"),
	      ncol = 3))
	dev.off()


