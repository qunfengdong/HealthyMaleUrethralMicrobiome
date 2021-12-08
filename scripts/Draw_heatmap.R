# Author: Yue Xing (yxing4@luc.edu)

library(circlize)
library(ComplexHeatmap)
library(vegan)
library(cluster)
library(phyloseq)
require(clusterSim)

pam.clustering=function(x,k) {
     require(cluster)
     cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
     return(cluster)
}

euclid = function(m) {
	dist(m, method = "euclidean")
}

enterotype = function(data,nk,in_dname,cl_fun) {
	if (cl_fun=="euclid") {
		data.dist=euclid(t(data))
	}
	
	dd=as.matrix(data.dist,labels=TRUE)
	a=rownames(dd)

	nclusters=NULL

	for (k in 1:20) { 
		if (k==1) {
			nclusters[k]=NA 
		} else {
			data.cluster_temp=pam.clustering(data.dist, k)
			nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist,
			centrotypes = "medoids")
		}
	}

	png(paste0(in_dname,".nclusters.png"),width=2000,height=2000,res=300)
	print(plot(nclusters, type="h", xlab="k clusters", ylab="CH index"))
	dev.off()

	# 
	data.cluster=pam.clustering(data.dist, k=nk)

	b=data.frame(ID=a,clus=data.cluster)

	obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])

	ent=list("cluster" = b, "dis" = data.dist)
	return(ent)

}

preprocess = function(file,qc,cot,transf,in_nk,
	in_dir,ra) {
	
	db=read.csv(paste0(file,".csv"),row.names=1)

	qc=qc[rownames(qc)%in%colnames(db),]

	cb=cot[rownames(cot)%in%rownames(db),]
	cb=cb[,colnames(cb)%in%colnames(db)]
	#print(dim(cb)==dim(db))

	sumb=data.frame(apply(cb,2,sum))
	colnames(sumb)="Taxa"
	qc=merge(qc,sumb,by="row.names",all.x=TRUE)
	rownames(qc)=qc[,1]
	qc=qc[,-1]

	#

	qcb=qc[colnames(db),]
	db=as.matrix(db)

	#print(which(rownames(qcb)!=colnames(db)))

	if (transf=="log2") {
		db2=log2(db+1)
	} else if (transf=="log10") {
		db2=log10(db+1)
	} else if (transf=="loge") {
		db2=log(db+1)
	} else if (transf=="none") {
		db2=ra[rownames(db),colnames(db)]
	} else if (transf=="clr") {
		db2=ra[rownames(db),colnames(db)]
		ob=otu_table(db, taxa_are_rows=TRUE)
		psb=phyloseq(ob)
		psbc <- microbiome::transform(psb, 'clr')
		db=as.matrix(psbc)
	} else if (transf=="human") {
		db2=ra[rownames(db),colnames(db)]
		for (i in 1:ncol(db)) {
			a=colnames(db)[i]
			db[,i]=(db[,i]+0.1)/qcb[a,"Human_reads"]
		}
		db=log(db)
	}

	dname=paste0(in_dir,"//",file,".data_trans_",transf,".k_",in_nk)

	ent_b=enterotype(db,in_nk,dname,"euclid")

	write.csv(ent_b$cluster,paste0(dname,".cluster.csv"),quote=FALSE)

	nl <- list("db2" = db2, "qcb" = qcb, 
		"db"=db, 
		"dname"=dname,
		"b_clus"=ent_b$cluster, "b_dist"=ent_b$dis)

	assign("data.list", nl, envir = .GlobalEnv)
}

plot_heatmap = function(pref,d1,ta,cdr,cdc,cf,wd,ht,rdw,cdw) {
	png(paste0(pref,".png"),width=wd,height=ht,res=300)
	par(mar = c(2, 2, 0.5, 0.5))
	rownames(d1)=gsub("_"," ",rownames(d1))
	ht = Heatmap(d1, name = "Abundance",
		top_annotation = ta,
	    clustering_distance_rows = cdr,
	    clustering_distance_columns = cdc,
	    #column_names_rot = 45,
	    row_names_gp = gpar(fontsize = 10),
	    #row_names_side = "left",
	   	column_names_gp = gpar(fontsize = 10),
	   	row_title_gp = gpar(fontsize = 12),
	   	column_title_gp = gpar(fontsize = 12),
	    column_title = "Subjects", 
	    column_names_side = "top",
	    row_title = "Relative abundance",
	    row_dend_width = unit(rdw, "mm"),
	    column_dend_height = unit(cdw, "mm"),
	    row_dend_side = "right",
	    col = cf
	)
	ht=draw(ht)
	print(ht)

	dev.off()
	#print(column_order(ht))
}

add_bars_b = function(dl0,in_bar_file) {
	b2=read.csv(in_bar_file)
	b2=b2[b2[,1]%in%rownames(dl0$qcb),]
	rownames(b2)=b2[,1]
	b2=b2[rownames(dl0$qcb),]
	return(b2[,2])
}

add_bar_clus_b = function(dl0) {
	b0=dl0$b_clus
	b0=b0[b0[,1]%in%rownames(dl0$qcb),]
	rownames(b0)=b0[,1]
	b0=b0[rownames(dl0$qcb),]
	return(b0[,2])
}

plot_ht_all = function(dl,transf_bar) {

	db2=dl$db2
	qcb=dl$qcb

	col_fun = colorRamp2(c(0, 0.01, 0.1, 1, 100), 
		c("blue", "skyblue", "yellow", "orange", "red"))
	
	col_list=list(
		Urethraltype=c("1"="steelblue1","2"="deeppink1","3"="purple2",
			"4"="khaki1","5"="orange","6"="black","7"="grey"),
        Last_60_days=c("N"="grey","V"="blue","R"="green","VR"="red"),
        Last_year=c("N"="grey","V"="blue","R"="green","VR"="red"),
        Ever=c("N"="grey","V"="blue","R"="green","VR"="red"))

	bar1_b=add_bars_b(dl,
		"E://IUMP//sex_behavior_and _bacteria_new_list//new_list_11022020//bar_60d.csv")
	bar0_b=add_bar_clus_b(dl)

	bar2_b=add_bars_b(dl,
		"E://IUMP//sex_behavior_and _bacteria_new_list//new_list_11022020//bar_1yr.csv")
	
	bar3_b=add_bars_b(dl,
		"E://IUMP//sex_behavior_and _bacteria_new_list//new_list_11022020//bar_ever.csv")
	
	if (transf_bar=="log10") {
		breaks = c(1, 10, 100)

		abundance_fun_b=anno_barplot(log10(qcb$Taxa/qcb$Human_reads*1000+1), 
				axis_param = list(at = log10(breaks), labels =breaks),
				gp=gpar(fill=c("white")),
				height = unit(4, "cm"))

	} else if (transf_bar=="log2") {
		breaks = c(1, 10, 100)

		abundance_fun_b=anno_barplot(log2(qcb$Taxa/qcb$Human_reads*1000+1), 
				axis_param = list(at = log2(breaks), labels =breaks),
				gp=gpar(fill=c("white")),
				height = unit(4, "cm"))

	} else if (transf_bar=="log") {
		breaks = c(1, 10, 100)

		abundance_fun_b=anno_barplot(log(qcb$Taxa/qcb$Human_reads*1000+1), 
				axis_param = list(at = log(breaks), labels =breaks),
				gp=gpar(fill=c("white")),
				height = unit(4, "cm"))

	}

	column_taxa = HeatmapAnnotation(
		Urethraltype=bar0_b,
		Last_60_days=bar1_b,
		Last_year=bar2_b,
		Ever=bar3_b,
		Abundancex1000=abundance_fun_b,
		col = col_list)

	plot_heatmap(paste0(dl$dname,".bar_trans_",transf_bar),
		db2,column_taxa,euclid,euclid,
		col_fun,6000,5000,30,50)

}

plot_final = function(file0,transf0,nk0,transf_bar0,
	dir0,sample_clust,qc0,cot0,ra0) {
	dir.create(dir0)
	preprocess(file0,qc0,cot0,transf0,nk0,dir0,ra0)
	plot_ht_all(data.list,transf_bar0)

}
