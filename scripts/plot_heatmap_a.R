pam.clustering=function(x,k) {
     require(cluster)
     cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
     return(cluster)
}

euclid = function(m) {
	dist(m, method = "euclidean")
}

enterotype = function(data,nk,in_dname,cl_fun) {
	data.dist=euclid(t(data))
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

	#nclusters[is.na(nclusters)]=0

	png(paste0(in_dname,".nclusters.png"),width=2000,height=2000,res=300)
	print(plot(nclusters, type="h", xlab="k clusters", ylab="CH index"))
	dev.off()

	# 
	data.cluster=pam.clustering(data.dist, k=nk)

	b=data.frame(ID=a,clus=data.cluster)

	obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
	print(obs.silhouette)

	ent=list("cluster" = b, "dis" = data.dist)
	return(ent)

}

preprocess = function(file,qc,cot,in_thres,transf,in_nk,
	in_dir,in_cl_fun,renorm,ra) {
	db=read.csv(paste0(file,".csv"),row.names=1)

	qc=qc[rownames(qc)%in%colnames(db),]
	print(dim(qc))

	cb=cot[rownames(cot)%in%rownames(db),]
	cb=cb[,colnames(cb)%in%colnames(db)]
	print(dim(cb)==dim(db))

	sumb=data.frame(apply(cb,2,sum))
	colnames(sumb)="Taxa"
	qc=merge(qc,sumb,by="row.names",all.x=TRUE)
	rownames(qc)=qc[,1]
	qc=qc[,-1]
	qcb=qc[colnames(db),]

	db=as.matrix(db)

	print(which(rownames(qcb)!=colnames(db)))
	print(dim(db))

	if (transf=="clr") {
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

	dname=paste0(in_dir,"//",file,".thres_",in_thres,".data_trans_",transf,".k_",in_nk)

	ent_b=enterotype(db,in_nk,dname,in_cl_fun)

	write.csv(ent_b$cluster,paste0(dname,".cluster.csv"),quote=FALSE)

	nl <- list("db2" = db2, "qcb" = qcb, 
		"db"=db, 
		"dname"=dname,
		#"dname_b"=dname_b, "dname_v"=dname_v,
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

enterotype_clus_b = function(m) {
	data.list$b_dist
}

plot_ht_all = function(dl,in_cdr,transf_bar,col_f) {

	db2=dl$db2
	qcb=dl$qcb

	nn=0.1

	col_fun6 = colorRamp2(c(0, 0.01, 0.1, 1, 100), 
		c("blue", "skyblue", "yellow", "orange", "red"))

	col_list=list(
		Urethraltype=c("1"="steelblue1","2"="deeppink1","3"="purple2",
			"4"="khaki1","5"="orange","6"="black","7"="grey"),
		Urethraltype_ra=c("1"="steelblue1","2"="deeppink1","3"="purple2",
			"4"="khaki1","5"="orange","6"="black","7"="grey"),
		Urethraltype_clr=c("1"="steelblue1","2"="deeppink1","3"="purple2",
			"4"="khaki1","5"="orange","6"="black","7"="grey"),
		Urethraltype_human=c("1"="steelblue1","2"="deeppink1","3"="purple2",
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

	####
	bar_clr=add_bars_b(dl,
		"E://IUMP//clr_080521//hea_heatmap_bacteria_count.thres_0.data_trans_clr.k_2.cluster.csv")

	bar_hn=add_bars_b(dl,
		"E://IUMP//clr_080521//hea_heatmap_bacteria_count.thres_0.data_trans_human.k_2.cluster.csv")

	bar_ra=add_bars_b(dl,
		"E://IUMP//clr_080521//health_all.metaphlan3_species.rm.relative_abundance.rm_low.thres_0.data_trans_none.k_2.bacteria.cluster.csv")

	####
	
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
		Urethraltype_ra=bar_ra,
		Urethraltype_clr=bar_clr,
		Urethraltype_human=bar_hn,
		Last_60_days=bar1_b,
		Last_year=bar2_b,
		Ever=bar3_b,
		Abundancex1000=abundance_fun_b,
		col = col_list)

	if (col_f==6) {
		plot_heatmap(paste0(dl$dname,".bar_trans_",transf_bar),
			db2,column_taxa,in_cdr,enterotype_clus_b,
			col_fun6,6000,5500,30,50)
	}
}

plot_final = function(file0,thres0,transf0,nk0,transf_bar0,
	dir0,sample_clust,re_norm,qc0,cot0,in_col_f,ra0) {
	dir.create(dir0)
	preprocess(file0,qc0,cot0,
	thres0,transf0,nk0,dir0,sample_clust,re_norm,ra0)

	plot_ht_all(data.list,euclid,transf_bar0,in_col_f)
}
