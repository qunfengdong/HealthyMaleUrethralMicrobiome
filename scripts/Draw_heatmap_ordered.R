# Author: Yue Xing (yxing4@luc.edu)

# This code finds the sum of all taxa for an subj and add it to QC automatically.
# This script doesn't split bacteria or virus. Just input a table and produce a heatmap from that table.
# This script is used for relative abundance/counts.
# Filtering by a threshold still in use.
# Don't recommend doing re-normalization.
# No missing values allowed.

# Display RA on the heatmap. RA directly get from input ra table. Transformation still used for clustering for clr and human reads.
# Can use "NA" for skipping noise removal.
# If transf=="none", will display RA on heatmap.

library(taxize)
library(circlize)
library(ComplexHeatmap)
library(vegan)
library(cluster)
library(phyloseq)
require(clusterSim)

dist.JSD <- function(inMatrix,pseudocount=0.000001,...) {
	KLD <- function(x,y) sum(x *log(x/y))
	JSD <- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
	matrixColSize <- length(colnames(inMatrix))
	matrixRowSize <- length(rownames(inMatrix))
	colnames <- colnames(inMatrix)
	resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
	    
	inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))

	for(i in 1:matrixColSize) {
		for(j in 1:matrixColSize) { 
			resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
			as.vector(inMatrix[,j]))
		}
	}
	colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
	as.dist(resultsMatrix)->resultsMatrix
	attr(resultsMatrix, "method") <- "dist"
	return(resultsMatrix) 
}

weighted_jacc <- function(A) {
	sim.jac <- matrix(0, nrow=nrow(A), ncol=nrow(A))
	rownames(sim.jac) <- rownames(A)
	colnames(sim.jac) <- rownames(A)

	#weighted jaccard
	pairs <- t(combn(1:nrow(A), 2))
	for (i in 1:nrow(pairs)){
	  num <- sum(sapply(1:ncol(A), function(x)(min(A[pairs[i,1],x],A[pairs[i,2],x]))))
	  den <- sum(sapply(1:ncol(A), function(x)(max(A[pairs[i,1],x],A[pairs[i,2],x]))))
	  sim.jac[pairs[i,1],pairs[i,2]] <- num/den
	  sim.jac[pairs[i,2],pairs[i,1]] <- num/den  
	}
	sim.jac[which(is.na(sim.jac))] <- 0
	diag(sim.jac) <- 1

	dist.jac <- as.dist(1-sim.jac)
	return(dist.jac)
}

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
     require(cluster)
     cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
     return(cluster)
}

noise.removal_inhouse = function(d0,thres) {
	print(dim(d0))
	if (thres!="NA") {
		idx=numeric()
		for (i in 1:nrow(d0)) {
			s=(d0[i,]>thres)
			if (sum(s==TRUE)==0) {
				idx=c(idx,i)
			}
		}
		if (length(idx)!=0) {
			d0=d0[-idx,]
		}
	}

	print(dim(d0))
	return(d0)
}

## Define distance functions
braycurtis = function(m) {
	bc=vegdist(m, method="bray", binary=FALSE)
	as.dist(as.matrix(bc))
}

braycurtis2 = function(m) {
	bc=vegdist(m, method="bray", binary=TRUE)
	as.dist(as.matrix(bc))
}

euclid = function(m) {
	dist(m, method = "euclidean")
}

jaccard = function(m) {
	distance(m, method = "jaccard")
}


phylo = function(m) {
	uids <- get_uid(rownames(m))
	if (NA %in% uids) {
		print("Error: some taxa names do not match NCBI IDs. Stop.")
		print(m[which(is.na(uids))])
		quit(status=1)
	}
	uid_class <- classification(uids, db = "ncbi")
	uid_tree <- class2tree(uid_class, check = TRUE)
	uid_tree$distmat
	#plot(uid_tree)
}

enterotype = function(data,nk,in_dname,cl_fun) {
	if (cl_fun=="bc") {
		data.dist=braycurtis(t(data))
	} else if (cl_fun=="euclid") {
		data.dist=euclid(t(data))
	} else if (cl_fun=="jsd") {
		data.dist=dist.JSD(data)
	} else if (cl_fun=="bc_b") {
		data.dist=braycurtis2(t(data))
	} else if (cl_fun=="jacc") {
		data.dist=jaccard(t(data))
	} else if (cl_fun=="jacc_w") {
		data.dist=weighted_jacc(t(data))
	}
	
	
	#print(data.dist)
	#print("flag")
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
	print(db[1:5,1:5])

	# unclassified removed. 
	# Actually all are 0 for unclassified in IUMP healthy samples.
	#if ("unclassified" %in% rownames(d1)) {
	#	d1=d1[-which(rownames(d1)=="unclassified"),]
	#}

	#for (i in 1:ncol(db)) {
	#	db[,i]=db[,i]/sum(db[,i])*100
	#}

	#

	db=noise.removal_inhouse(db,in_thres)

	if (renorm=="norm") {
		for (i in 1:ncol(db)) {
			db[,i]=db[,i]/sum(db[,i])*100
		}

		db[is.na(db)]=0

		idxb=numeric()
		for (i in 1:ncol(db)) {
			if (sum(db[,i])==0) {idxb=c(idxb,i)}
		}
		if (length(idxb)!=0) {
			db=db[,-idxb]
		}

	} else if (renorm=="add_other") {
		flag00=0
	} else if (renorm=="none") {
		flag00=0
	}

	#

	#qc=read.csv(in_qc,row.names=1)
	qc=qc[rownames(qc)%in%colnames(db),]
	print(dim(qc))
	#rownames(qc)=gsub("\\..*","",rownames(qc))

	cb=cot[rownames(cot)%in%rownames(db),]
	cb=cb[,colnames(cb)%in%colnames(db)]
	print(dim(cb)==dim(db))

	sumb=data.frame(apply(cb,2,sum))
	#rownames(sumb)=gsub("\\..*","",rownames(sumb))
	colnames(sumb)="Taxa"
	qc=merge(qc,sumb,by="row.names",all.x=TRUE)
	rownames(qc)=qc[,1]
	qc=qc[,-1]

	#

	qcb=qc[colnames(db),]
	#print(head(qcb))

	##

	db=as.matrix(db)

	print(which(rownames(qcb)!=colnames(db)))
	print(dim(db))

	# db2 and dv2 for heatmap, db and dv for sample clustering
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

	dname=paste0(in_dir,"//",file,".thres_",in_thres,".data_trans_",transf,".k_",in_nk)
	#dname_b=paste0(dname,".bacteria")
	#dname_v=paste0(dname,".virus")

	#print(db)
	ent_b=enterotype(db,in_nk,dname,in_cl_fun)

	#dis_tmp <- as.data.frame(as.matrix(ent_b$dis))
	#write.csv(dis_tmp,paste0(dname,".ent_b.dis.csv"),quote=FALSE)

	write.csv(ent_b$cluster,paste0(dname,".cluster.csv"),quote=FALSE)

	nl <- list("db2" = db2, "qcb" = qcb, 
		"db"=db, 
		"dname"=dname,
		#"dname_b"=dname_b, "dname_v"=dname_v,
		"b_clus"=ent_b$cluster, "b_dist"=ent_b$dis)

	assign("data.list", nl, envir = .GlobalEnv)
	#return(nl)

}

plot_heatmap = function(pref,d1,ta,cdr,cf,wd,ht,rdw,cdw,li) {
	png(paste0(pref,".png"),width=wd,height=ht,res=300)
	par(mar = c(2, 2, 0.5, 0.5))
	rownames(d1)=gsub("_"," ",rownames(d1))
	ht = Heatmap(d1, name = "Abundance",
		top_annotation = ta,
	    clustering_distance_rows = cdr,
	    column_order = order(match(colnames(d1),li)),
	    #clustering_distance_columns = cdc,
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
	#draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
	dev.off()
}

add_bars_b = function(dl0,in_bar_file) {
	b2=read.csv(in_bar_file)
	#print(dim(b2))
	#print(dim(dl0$qc))
	b2=b2[b2[,1]%in%rownames(dl0$qcb),]
	rownames(b2)=b2[,1]
	b2=b2[rownames(dl0$qcb),]
	print(dim(b2))
	print(which(b2[,1]!=colnames(dl0$db2)))
	#print(which(b2[,1]!=colnames(dl0$dv2)))
	print(which(rownames(dl0$qcb)!=colnames(dl0$db2)))
	#print(which(rownames(dl0$qcv)!=colnames(dl0$dv2)))
	return(b2[,2])
}

add_bar_clus_b = function(dl0) {
	b0=dl0$b_clus
	#print(dim(b0))
	#print(dim(dl0$qcb))
	b0=b0[b0[,1]%in%rownames(dl0$qcb),]
	rownames(b0)=b0[,1]
	b0=b0[rownames(dl0$qcb),]
	print(dim(b0))
	print(which(b0[,1]!=colnames(dl0$db2)))
	#print(which(b0[,1]!=colnames(dl0$dv2)))
	print(which(rownames(dl0$qcb)!=colnames(dl0$db2)))
	#print(which(rownames(dl0$qcb)!=colnames(dl0$dv2)))
	return(b0[,2])
}

enterotype_clus_b = function(m) {
	#print(dim(as.matrix(data.list$b_dist)))
	data.list$b_dist
}

plot_ht_all = function(dl,in_cdr,transf_bar,col_f,in_li) {

	db2=dl$db2
	qcb=dl$qcb

	nn=0.1
	#db2[db2<nn]=0

	col_fun = colorRamp2(c(0, nn, 1, 100), 
		c("grey98", "orange", "red2", "red4"))

	#col_bar = colorRamp2(c(0, 1e10), c("green", "green"))

	col_fun2 = colorRamp2(c(min(db2), 0, max(db2)), 
		c("blue", "white", "red"))
	
	col_fun3 = colorRamp2(c(0, max(db2)), 
		c("blue", "red"))

	col_fun4 = colorRamp2(c(min(db2), -18.42, 0), 
		c("blue", "white","red"))

	col_fun5 = colorRamp2(c(min(db2), -0.00001,0, max(db2)), 
		c("blue", "blue", "white","red"))

	col_fun6 = colorRamp2(c(0, 0.01, 0.1, 1, 100), 
		c("blue", "skyblue", "yellow", "orange", "red"))
	
	col_list=list(
		Urethraltype=c("1"="steelblue1","2"="deeppink1","3"="purple2",
			"4"="khaki1","5"="orange","6"="black","7"="grey"),
        Last_60_days=c("N"="grey","V"="blue","R"="green","VR"="red"),
        Last_year=c("N"="grey","V"="blue","R"="green","VR"="red"),
        Ever=c("N"="grey","V"="blue","R"="green","VR"="red"))

	bar1_b=add_bars_b(dl,
		"bar_60d.csv")
	bar0_b=add_bar_clus_b(dl)

	bar2_b=add_bars_b(dl,
		"bar_1yr.csv")
	
	bar3_b=add_bars_b(dl,
		"bar_ever.csv")
	
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

	#print(dim(db2))
	if (col_f==1) {
		plot_heatmap(paste0(dl$dname,".bar_trans_",transf_bar),
			db2,column_taxa,in_cdr,
			col_fun,6000,5000,30,50,in_li)
	} else if (col_f==2) {
		plot_heatmap(paste0(dl$dname,".bar_trans_",transf_bar),
			db2,column_taxa,in_cdr,
			col_fun2,6000,5000,30,50,in_li)
	} else if (col_f==3) {
		plot_heatmap(paste0(dl$dname,".bar_trans_",transf_bar),
			db2,column_taxa,in_cdr,
			col_fun3,6000,5000,30,50,in_li)
	} else if (col_f==4) {
		plot_heatmap(paste0(dl$dname,".bar_trans_",transf_bar),
			db2,column_taxa,in_cdr,
			col_fun4,6000,5000,30,50,in_li)
	} else if (col_f==5) {
		plot_heatmap(paste0(dl$dname,".bar_trans_",transf_bar),
			db2,column_taxa,in_cdr,
			col_fun5,6000,5000,30,50,in_li)
	} else if (col_f==6) {
		plot_heatmap(paste0(dl$dname,".bar_trans_",transf_bar),
			db2,column_taxa,in_cdr,
			col_fun6,6000,5000,30,50,in_li)
	}
}

plot_final = function(file0,thres0,transf0,nk0,transf_bar0,
	taxa_clust,dir0,sample_clust,re_norm,qc0,cot0,
	in_col_f,ra0,in_li0) {
	dir.create(dir0)
	preprocess(file0,
	qc0,cot0,
	thres0,transf0,nk0,dir0,sample_clust,re_norm,ra0)

	if (taxa_clust=="phylo") {
		plot_ht_all(data.list,phylo,transf_bar0,in_col_f,in_li0)
	} else if (taxa_clust=="euclid") {
		plot_ht_all(data.list,euclid,transf_bar0,in_col_f,in_li0)
	}

	#rm(data.list, pos = ".GlobalEnv")
}
