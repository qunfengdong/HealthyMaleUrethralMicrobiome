# Author: Yue Xing (yxing4@luc.edu)

calc_OR_cov = function(nm) {

	d2=read.csv(nm,row.names=1,na.strings=".")
	d2$Race=factor(d2$Race)
	d2$Race=relevel(d2$Race,ref="Other")

	if (length(table(d2[,ncol(d2)-9])) <2) {
		print(paste0("Error: ",nm))
		tt=data.frame(t(rep(NA,9)))
		colnames(tt)=c("Type","Cov","OR","CI_2.5","CI_97.5","Estimate","Std_err","z","p_val")
		tabe=data.frame(t(rep(NA,6)))
		colnames(tabe)=c("Estimate","Std_err","z","p_val","Type","Cov")
		li = list(tt=tt,tabe=tabe)
		return(li)
	} else {
		flag_all=c("ars_l","ars_60")
		
		tab=data.frame()
		tabb=data.frame()
		tabe=data.frame()
		for (flag in flag_all) {
			print(flag)
			for (i in 1:(ncol(d2)-10)) {
				if (length(table(d2[,i])) <2) {
					print(colnames(d2)[i])
					next
				}

				x=factor(d2[,(ncol(d2)-9)])
				tp=colnames(d2)[(ncol(d2)-9)]
				leveln=names(table(x))
				if (length(leveln)<2) {
					print(tp)
					next
				}
			
				x=relevel(x,ref="0")
				y=factor(d2[,i])
				y=relevel(y,ref="1")

				if (flag=="ars_l") {
					lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"], family = "binomial")
				} else if (flag=="ars_60") {
					lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_60d"], family = "binomial")
				}
		
				rm(a)
				rm(b)
				tryCatch(
					expr={a=data.frame(exp(cbind(OR = coef(lmodel), confint(lmodel))))},
					error = function(e){print(paste("Error_a",colnames(d2)[i]))}
					)
				tryCatch(
					expr={b=data.frame(summary(lmodel)$coefficients)},
					error = function(e){print(paste("Error_b",colnames(d2)[i]))}
					)

				if (!exists("a")) {
					a=data.frame(t(rep(NA,5)))
					colnames(a)=c("tp","flag","OR","X2.5..","X97.5..")
				} else {
					a=a[rownames(a)%in%paste0("x",leveln),]
					if (nrow(a)==0) {
						print("glm error")
						exit(status=1)
					}
					print(rownames(a))
					a=cbind(tp,flag,a)
				}
				tab=rbind(tab,a)

				if (!exists("b")) {
					d=data.frame(t(rep(NA,6)))
					colnames(d)=c("tp","flag","Estimate","Std..Error","z.value","Pr...z..")
				} else {
					d=b[rownames(b)%in%paste0("x",leveln),]
					if (nrow(d)==0) {
						print("glm error")
						exit(status=1)
					}
					print(rownames(d))
					d=cbind(tp,flag,d)

					e=b[!(rownames(b)%in%paste0("x",leveln)),]
					e$tp=tp
					e$flag=flag
				}
				tabb=rbind(tabb,d)
				tabe=rbind(tabe,e)

			}
		}

		if(nrow(tab)==0) {
			tab=data.frame(t(rep(NA,5)))
		}

		if(nrow(tabb)==0) {
			tabb=data.frame(t(rep(NA,6)))
		}
		
		if(nrow(tabe)==0) {
			tabe=data.frame(t(rep(NA,6)))
		}

		colnames(tab)=c("Type","Cov","OR","CI_2.5","CI_97.5")
		colnames(tabb)=c("Type","Cov","Estimate","Std_err","z","p_val")
		colnames(tabe)=c("Estimate","Std_err","z","p_val","Type","Cov")

		rownames(tab)=paste(tab[,1],tab[,2],sep="_")
		rownames(tabb)=paste(tabb[,1],tabb[,2],sep="_")

		tt=merge(tab,tabb,by="row.names")

		tt=tt[,c("Type.x","Cov.x","OR","CI_2.5","CI_97.5","Estimate","Std_err","z","p_val")]
		colnames(tt)[1:2]=c("Type","Cov")

		li = list(tt=tt,tabe=tabe)
		return(li)
	}
}
