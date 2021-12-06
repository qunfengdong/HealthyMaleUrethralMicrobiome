# ALR/CLR transformation and draw heatmap
# Author: Yue Xing (yxing4@luc.edu)

	source("scripts//heatmap.R")

	qc0=readRDS("healthy_QC.rds")
	cot0=readRDS("healthy_count.rds")
	ra0=read.csv("health_all.metaphlan3_species.ra.csv",row.names=1)

	transfs="human"
	#transfs="clr"
	transf_bars="log"
	files="hea_heatmap_bacteria_count"
	ks=2

	plot_final(files,transfs,ks,transf_bars,
		"out",qc0=qc0,cot0=cot0,ra0=ra0)

# Reorder samples and draw dendograms for heatmap: not shown, can be sent upon request.

