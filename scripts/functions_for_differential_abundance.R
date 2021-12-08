# Author: Yue Xing (yxing4@luc.edu)

# Wilcoxon (input ra, alr, etc.) with additional information
	wilcoxon_add_info = function(da,dra,gp) {
		print(all(rownames(da)==rownames(dra)))
		print(all(colnames(da)==colnames(dra)))

		da=as.matrix(da)
		dra=as.matrix(dra)
		tab=character()
		u1=gp[gp[,2]==1,]
		u2=gp[gp[,2]==2,]
		ut1=da[,colnames(da)%in%u1[,1]]
		ut2=da[,colnames(da)%in%u2[,1]]
		ut1r=dra[,colnames(dra)%in%u1[,1]]
		ut2r=dra[,colnames(dra)%in%u2[,1]]
		ut1=da[,u1[,1]]
		ut2=da[,u2[,1]]
		ut1r=dra[,u1[,1]]
		ut2r=dra[,u2[,1]]

		for (i in 1:nrow(da)) {
			tt=wilcox.test(ut1[i,],ut2[i,], paired = FALSE, exact = FALSE)
			tt2=t.test(ut1[i,],ut2[i,], paired = FALSE)

			s1=length(which(ut1r[i,]!=0))
			s2=length(which(ut2r[i,]!=0))
			s1=s1/length(ut1r[i,])
			s2=s2/length(ut2r[i,])

			# add mean and sd of input table, median
			# add ra and prevalence
			tab=rbind(tab,c(rownames(da)[i],tt$p.value,tt2$p.value,
				mean(ut1[i,]),sd(ut1[i,]),mean(ut2[i,]),sd(ut2[i,]),
				median(ut1[i,]),median(ut2[i,]),
				mean(ut1r[i,]),sd(ut1r[i,]),mean(ut2r[i,]),sd(ut2r[i,]),
				s1,s2))
		}
		colnames(tab)=c("Taxon","Wilcoxon","t_test","U1_mean","U1_sd",
			"U2_mean","U2_sd","U1_median","U2_median",
			"U1_ra_mean","U1_ra_sd","U2_ra_mean","U2_ra_sd",
			"U1_prevalence","U2_prevalence")

		tab=data.frame(tab)
		tab$wx.BH=p.adjust(tab$Wilcoxon,method="BH")
		tab$t.BH=p.adjust(tab$t_test,method="BH")

		return(tab)
	}

