# ALR/CLR transformation and heatmap
# Author: Yue Xing (yxing4@luc.edu)

# Order the samples

	source("scripts//Draw_heatmap.R")

	qc0=readRDS("healthy_QC.rds")
	cot0=readRDS("healthy_count.rds")
	ra0=read.csv("health_all.metaphlan3_species.ra.csv",row.names=1)

	# For CLR transformed data, use <transfs="clr">
	transfs="human"
	thress=0
	transf_bars="log"
	files="hea_heatmap_bacteria_count"
	spcl="euclid"
	ks=2

	plot_final(files,thress,transfs,ks,transf_bars,
	taxa_clust="phylo",
	"k2_draw_alr",
	sample_clust=spcl,
	re_norm="none",
	qc0=qc0,
	cot0=cot0,
	in_col_f=6,
	ra0=ra0)

	d1=read.csv("hea_heatmap_bacteria_count.thres_0.data_trans_human.k_2.cluster.csv",row.names=1)
	db=data.list$db

	dis=dist(t(db),method = "euclidean")
	g1=d1[d1[,2]==1,1]
	g2=d1[d1[,2]==2,1]

	n=2
	matrixColSize <- n
	#matrixRowSize <- n
	colnames <- c("g1","g2")
	resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
	colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)

	resultsMatrix[2,1]=dist_between_centroids(dis, g1, g2)

	as.dist(resultsMatrix)->resultsMatrix
	attr(resultsMatrix, "method") <- "dist"

	h=hclust(resultsMatrix)
	s01=h$labels
	s02=h$order
	aa=s01[s02]
	# change g1-g2 to 1-2 in aa
	aa=gsub("g","",aa)
	aa=as.numeric(as.character(aa))

# Draw h-clusters

	li=character()
	for (gr in aa) {
		d2=d1[d1[,2]==gr,]
		db2=db[,colnames(db)%in%d2[,1]]
		dis2=dist(t(db2), method = "euclidean")
		h2=hclust(dis2)

		if (gr==1) {
			wd=2000
		} else {
			wd=4000
		}

		png(paste0("euclid_hclust.",gr,".k2.png"),height=2000,width=wd,res=300)
		plot(h2,hang=-1)
		dev.off()
		s1=h2$labels
		s2=h2$order
		li=c(li,s1[s2])
	}

# Draw final heatmap

	source("scripts//Draw_heatmap_ordered.R")

	# For CLR transformed data, use <transfs="clr">
	transfs="human"
	thress=0
	transf_bars="log"
	files="hea_heatmap_bacteria_count"
	spcl="euclid"
	ks=2

	plot_final(files,thress,transfs,ks,transf_bars,
	taxa_clust="phylo",
	"k2_draw_alr",
	sample_clust=spcl,
	re_norm="none",
	qc0=qc0,
	cot0=cot0,
	in_col_f=6,
	ra0=ra0,
	in_li0=li
	)
