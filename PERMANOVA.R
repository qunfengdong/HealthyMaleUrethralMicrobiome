# PERMANOVA
# Author: Yue Xing (yxing4@luc.edu)

# PERMANOVA for meta data

	library(vegan)

	da=read.csv("hea_heatmap_bacteria_human.csv",row.names=1,na.strings=".")
	da=t(da)
	
	meta=read.csv("meta_for_permanova.csv",row.names=1,na.strings=".")
	meta=meta[rownames(meta)%in%rownames(da),]
	all(rownames(meta)==rownames(da))
	
	tab=character()
	for (i in 1:ncol(meta)) {
		rm(fit)
		da0=da
		x=meta[,i]
		x[x==""]=NA

		idx=which(is.na(x))
		if (length(idx)>0) {
			x=x[-idx]
			da0=da0[-idx,]
		}
		BC.dist=dist(da0, method = "euclidean")

		if (length(table(x))>1) {
			tryCatch(
				expr = {fit=adonis(BC.dist ~ x, data = meta, permutations = 1000)},
				error = function(e){print(paste("Error",colnames(meta)[i]))}
				)
				
			if (exists("fit")) {
				tab=rbind(tab,c(colnames(meta)[i],length(idx),
				fit$aov.tab["x","Df"],
				fit$aov.tab["x","R2"],
				fit$aov.tab["x","Pr(>F)"]))
			}
		}
	}
	
	colnames(tab)=c("Question","Num.missing","Df","R2","pvalue")
	tab=data.frame(tab)
	tab$adjp=p.adjust(tab$pvalue, method = "BH")
	write.csv(tab,"permanova_meta_euclid.csv",quote=FALSE,row.names=FALSE)

# PERMANOVA for sex behaviors
## Integrated sex behaviors

	meta=read.csv("meta.csv",row.names=1,na.strings=".")

	da=read.csv("hea_heatmap_bacteria_human.csv",row.names=1,na.strings=".")
	da=t(da)
	da=da[rownames(da)%in%rownames(meta),]
	all(rownames(meta)==rownames(da))

	dis="euclidean"
	BC.dist=dist(da, method = "euclidean")

	set.seed(1)

	####

	fit=adonis(BC.dist ~ X60d, data = meta, permutations = 1000)
	fit
	
	fit=adonis(BC.dist ~ X1yr, data = meta, permutations = 1000)
	fit
	
	fit=adonis(BC.dist ~ Xer, data = meta, permutations = 1000)
	fit

## Sex behavior split to categories

	meta=read.csv("meta_split.csv",row.names=1,na.strings=".")
	which(is.na(meta))

	da=read.csv("hea_heatmap_bacteria_human.csv",row.names=1,na.strings=".")
	da=t(da)
	all(rownames(meta)==rownames(da))

	dis="euclidean"
	BC.dist=dist(da, method = "euclidean")

	set.seed(1)

	####

	fit=adonis(BC.dist ~ X60d_V+X60d_O+X60d_R+X60d_OV+X60d_OR+X60d_OVR+X60d_N, data = meta, permutations = 1000)
	fit
	
	fit=adonis(BC.dist ~ X1yr_V+X1yr_O+X1yr_R+X1yr_OV+X1yr_OR+X1yr_OVR+X1yr_N, data = meta, permutations = 1000)
	fit

	fit=adonis(BC.dist ~ Xer_O+Xer_OV+Xer_OR+Xer_OVR, data = meta, permutations = 1000)
	fit

## Independent sex behaviors

	meta=read.csv("meta.csv",row.names=1,na.strings=".")

	da=read.csv("hea_heatmap_bacteria_human.csv",row.names=1,na.strings=".")
	da=t(da)
	da=da[rownames(da)%in%rownames(meta),]
	all(rownames(meta)==rownames(da))

	dis="euclidean"
	BC.dist=dist(da, method = "euclidean")

	set.seed(1)

	####

	fit=adonis(BC.dist ~ O_60d + V_60d + R_60d, data = meta, permutations = 1000)
	fit

	fit=adonis(BC.dist ~ O_1yr + V_1yr + R_1yr, data = meta, permutations = 1000)
	fit

	fit=adonis(BC.dist ~ O_Ever + V_Ever + R_Ever, data = meta, permutations = 1000)
	fit


