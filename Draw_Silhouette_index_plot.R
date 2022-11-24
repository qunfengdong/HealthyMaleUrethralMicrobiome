# Draw Silhouette index plot
# Author: Yue Xing (yxing4@luc.edu)

	# The functions from Draw_heatmap.R were used for Silhouette index
	library(taxize)
	library(circlize)
	library(ComplexHeatmap)
	library(vegan)
	library(cluster)
	library(phyloseq)
	require(clusterSim)
	
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
	
		ent=list("cluster" = b, "dis" = data.dist, silhouette = obs.silhouette)
		return(ent)
	
	}
	
	preprocess = function(file,qc,cot,in_thres,transf,in_nk,
		in_dir,in_cl_fun,renorm,ra) {
		db=read.csv(paste0(file,".csv"),row.names=1)
		print(db[1:5,1:5])
	
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
			"b_clus"=ent_b$cluster, "b_dist"=ent_b$dis,
			"silhouette"=ent_b$silhouette)
	
		assign("data.list", nl, envir = .GlobalEnv)
		#return(nl)
	
	}

# Draw Silhouette index plot

	qc0=readRDS("healthy_QC.rds")
	cot0=readRDS("healthy_count.rds")
	ra0=read.csv("health_all.metaphlan3_species.ra.csv",row.names=1)
	
	tab=numeric()
	for (ks in c(2:20)) {
		rm(data.list)
		print(ks)
		preprocess("hea_heatmap_bacteria_count",qc0,cot0,0,"human",ks,
			"test_index","euclid","none",ra0)
	
		sil=data.list$silhouette
		tab=rbind(tab,c(ks,sil))
	}

	nclusters=NULL
	for (k in 1:20) { 
		if (k==1) {
			nclusters[k]=NA 
		} else {
			nclusters[k]=tab[which(tab[,1]==k),2]
		}
	}

	png("Silhouette_index_plot.png",width=2000,height=2000,res=300)
	print(plot(nclusters, type="h", xlab="k clusters", ylab="Silhouette index"))
	dev.off()