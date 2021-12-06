# SpiecEasi for microbial network for IUMP samples
# Author: Yue Xing (yxing4@luc.edu)

	library(SpiecEasi)
	
	nm="hea_heatmap_bacteria_count"
	d1=t(as.matrix(read.csv(paste0(nm,".csv"),row.names=1)))
	d1=d1+0.1
	
	dmb <- spiec.easi(d1, method='mb', lambda.min.ratio=1e-3,
	nlambda=20, pulsar.params=list(rep.num=50))

	dg <- spiec.easi(d1, method='glasso', lambda.min.ratio=1e-3,
	nlambda=20, pulsar.params=list(rep.num=50))
		
	####

	a=as.matrix(getRefit(dmb))
	colnames(a)=rownames(a)=colnames(d1)
	
	tab=character()
	for (i in 1:(nrow(a)-1)) {
		for (j in (i+1):ncol(a)) {
			if (a[i,j]>0) {
				tab=rbind(tab,c(rownames(a)[i],colnames(a)[j],a[i,j]))
			}
		}
	}
	colnames(tab)=c("Taxon1","Taxon2","edge")
	
	write.csv(tab,"hea_heatmap_bacteria_count.SpiecEasi.mb.csv",quote=FALSE,row.names=FALSE)
	
	####

	a=as.matrix(getRefit(dg))
	colnames(a)=rownames(a)=colnames(d1)
	
	tab=character()
	for (i in 1:(nrow(a)-1)) {
		for (j in (i+1):ncol(a)) {
			if (a[i,j]>0) {
				tab=rbind(tab,c(rownames(a)[i],colnames(a)[j],a[i,j]))
			}
		}
	}
	colnames(tab)=c("Taxon1","Taxon2","edge")
	
	write.csv(tab,"hea_heatmap_bacteria_count.SpiecEasi.glasso.csv",quote=FALSE,row.names=FALSE)
	
	# Get weights (beta for MB and inverted icov for glasso)

	a=sub("[[:lower:]]*_"," ",colnames(d1))
	
	b=getOptBeta(dmb)
	b=as.matrix(b)
	rownames(b)=colnames(b)=a
	
	tab=character()
	for (i in 1:(nrow(b)-1)) {
		for (j in (i+1):ncol(b)) {
			if (b[i,j]!=0) {
				tab=rbind(tab,c(rownames(b)[i],colnames(b)[j],b[i,j]))
			}
		}
	}
	colnames(tab)=c("Taxon1","Taxon2","Beta")

	d=getOptMerge(dmb)
	d=as.matrix(d)
	rownames(d)=colnames(d)=a
	tab=data.frame(tab)
	tab$conf.score=NA
	for (s in 1:nrow(tab)) {
		tab[s,4]=d[tab[s,1],tab[s,2]]
	}

	write.csv(tab,"hea_heatmap_bacteria_count.SpiecEasi.mb.beta.csv",quote=FALSE,row.names=FALSE)

	####
	
	b=getOptiCov(dg)
	b=as.matrix(b)
	rownames(b)=colnames(b)=a
	
	tab=character()
	for (i in 1:(nrow(b)-1)) {
		for (j in (i+1):ncol(b)) {
			if (b[i,j]!=0) {
				# inverted icov!
				tab=rbind(tab,c(rownames(b)[i],colnames(b)[j],1/b[i,j]))
			}
		}
	}
	colnames(tab)=c("Taxon1","Taxon2","i_icov")

	d=getOptMerge(dg)
	d=as.matrix(d)
	rownames(d)=colnames(d)=a
	tab=data.frame(tab)
	tab$conf.score=NA
	for (s in 1:nrow(tab)) {
		tab[s,4]=d[tab[s,1],tab[s,2]]
	}

	write.csv(tab,"hea_heatmap_bacteria_count.SpiecEasi.glasso.i_icov.csv",quote=FALSE,row.names=FALSE)

