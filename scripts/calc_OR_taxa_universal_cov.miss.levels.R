# Author: Yue Xing (yxing4@luc.edu)

calc_OR_cov = function(nm,ct) {

	d1=read.csv(nm,row.names=1,na.strings=".")
	d2=numeric()
	
	if (ct==0) {
		# Y: 1 = bacteria, 0 = no bacteria
		for (i in 1:(ncol(d1)-10)) {
			line=d1[,i]
			for (j in 1:length(line)) {
				if (line[j]!=0) {line[j]=1}
			}
			d2=cbind(d2,line)
		}
	} else if (ct==0.1) {
		# Y: 1 = bacteria abundance > 0.1%, 0 = bacteria abundance <= 0.1%
		for (i in 1:(ncol(d1)-10)) {
			line=d1[,i]
			for (j in 1:length(line)) {
				if (line[j]>0.1) {
					line[j]=1
				} else {
					line[j]=0
				}
			}
			d2=cbind(d2,line)
		}
	}  else if (ct==0.5) {
		# Y: 1 = bacteria abundance > 0.5%, 0 = bacteria abundance <= 0.5%
		for (i in 1:(ncol(d1)-10)) {
			line=d1[,i]
			for (j in 1:length(line)) {
				if (line[j]>0.5) {
					line[j]=1
				} else {
					line[j]=0
				}
			}
			d2=cbind(d2,line)
		}
	} else if (ct==1) {
		# Y: 1 = bacteria abundance > 1%, 0 = bacteria abundance <= 1%
		for (i in 1:(ncol(d1)-10)) {
			line=d1[,i]
			for (j in 1:length(line)) {
				if (line[j]>1) {
					line[j]=1
				} else {
					line[j]=0
				}
			}
			d2=cbind(d2,line)
		}
	} else {
		print("nm error")
		quit(status=1)
	}

	d2=data.frame(d2)
	d2=cbind(d2,d1[,((ncol(d1)-9):ncol(d1))])
	colnames(d2)=colnames(d1)
	rownames(d2)=rownames(d1)

	d2$Race=factor(d2$Race)
	
	if (length(table(d2[,ncol(d2)-9])) <2) {
		print(paste0("Error: ",nm))
		tt=data.frame(t(rep(NA,13)))
		colnames(tt)=c("Taxa","Type","Cutoff","Cov","Test","Reference",
			"OR","CI_2.5","CI_97.5","Estimate","Std_err","z","p_val")
		tabe=data.frame(t(rep(NA,9)))
		colnames(tabe)=c("Estimate","Std_err","z","p_val",
			"Taxa","Type","Cutoff","Cov","Reference")
		li = list(tt=tt,tabe=tabe)
		return(li)
	} else {

		flag_all=c("3l_a","3l_b","3l_c","3l_d","3l_e",
			"360_a","360_b","360_c","360_d","360_e",
			"3l_abcde","360_abcde","3l_be","360_be",
			"3l_n","360_n")
		
		tab=tabb=tabe=data.frame()
		for (flag in flag_all) {
			print(flag)
			for (i in 1:(ncol(d2)-10)) {
				if (length(table(d2[,i])) <2) {
					print(colnames(d2)[i])
					next
				}
		
				taxa=gsub(".*__","",colnames(d2)[i])

				x=factor(d2[,(ncol(d2)-9)])
				tp=colnames(d2)[(ncol(d2)-9)]
				#print(tp)
				leveln=names(table(x))
				if (length(leveln)<2) {
					print(tp)
					next
				}
			
				levelref1=leveln[which(startsWith(leveln,"n") | leveln=="0")]
				if (length(levelref1)==0) {
					levelref1=leveln[which(!startsWith(leveln,"b"))][1]
				}
				print(levelref1)
				x=relevel(x,ref=levelref1)

				y=factor(d2[,i])
				y=relevel(y,ref="0")

				rm(lmodel)
				tryCatch(
					expr={
						if (flag=="3l_a") {
							lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"]+d2[,"a"], family = "binomial")
						} else if (flag=="3l_b") {
							lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"]+d2[,"b"], family = "binomial")
						} else if (flag=="3l_c") {
							lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"]+d2[,"c"], family = "binomial")
						} else if (flag=="3l_d") {
							lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"]+d2[,"d"], family = "binomial")
						} else if (flag=="3l_e") {
							lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"]+d2[,"e"], family = "binomial")
						} else if (flag=="360_a") {
							lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_60d"]+d2[,"a"], family = "binomial")
						} else if (flag=="360_b") {
							lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_60d"]+d2[,"b"], family = "binomial")
						} else if (flag=="360_c") {
							lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_60d"]+d2[,"c"], family = "binomial")
						} else if (flag=="360_d") {
							lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_60d"]+d2[,"d"], family = "binomial")
						} else if (flag=="360_e") {
							lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_60d"]+d2[,"e"], family = "binomial")
						} else if (flag=="360_abcde") {
							lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_60d"]+d2[,"a"]+d2[,"b"]+d2[,"c"]+d2[,"d"]+d2[,"e"], family = "binomial")
						} else if (flag=="3l_abcde") {
							lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"]+d2[,"a"]+d2[,"b"]+d2[,"c"]+d2[,"d"]+d2[,"e"], family = "binomial")
						} else if (flag=="3l_be") {
							lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"]+d2[,"b"]+d2[,"e"], family = "binomial")
						} else if (flag=="360_be") {
							lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_60d"]+d2[,"b"]+d2[,"e"], family = "binomial")
						} else if (flag=="360_n") {
							lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_60d"], family = "binomial")
						} else if (flag=="3l_n") {
							lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"], family = "binomial")
						}
					},
					error = function(e){print(paste("Error_lm",colnames(d2)[i]))}
				)
		
				if (!exists("lmodel")) {
					next
				} else {
					# do nothing, continue
					flag_here=1
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
					a=data.frame(t(rep(NA,9)))
					colnames(a)=c("taxa","tp","ct","flag","rownames(a)","levelref1","OR","X2.5..","X97.5..")
				} else {
					a=a[rownames(a)%in%paste0("x",leveln),]
					if (nrow(a)==0) {
						print("glm error")
						exit(status=1)
					}
					a=cbind(taxa,tp,ct,flag,rownames(a),levelref1,a)
				}
				tab=rbind(tab,a)

				if (!exists("b")) {
					d=data.frame(t(rep(NA,10)))
					colnames(d)=c("taxa","tp","ct","flag","rownames(d)",
						"levelref1","Estimate","Std..Error","z.value","Pr...z..")
				} else {
					d=b[rownames(b)%in%paste0("x",leveln),]
					if (nrow(d)==0) {
						print("glm error")
						exit(status=1)
					}
					d=cbind(taxa,tp,ct,flag,rownames(d),levelref1,d)

					e=b[!(rownames(b)%in%paste0("x",leveln)),]
					#e=data.frame(t(e))
					e$taxa=taxa
					e$tp=tp
					e$ct=ct
					e$flag=flag
					e$ref=levelref1
				}
				tabb=rbind(tabb,d)
				tabe=rbind(tabe,e)

				leveln=leveln[-which(leveln==levelref1)]

				while (length(leveln)>1) {
					levelref1=leveln[which(!startsWith(leveln,"b"))][1]
					print(levelref1)
					x=relevel(x,ref=levelref1)

					y=factor(d2[,i])
					y=relevel(y,ref="0")

					# lmodel
					if (flag=="3l_a") {
						lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"]+d2[,"a"], family = "binomial")
					} else if (flag=="3l_b") {
						lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"]+d2[,"b"], family = "binomial")
					} else if (flag=="3l_c") {
						lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"]+d2[,"c"], family = "binomial")
					} else if (flag=="3l_d") {
						lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"]+d2[,"d"], family = "binomial")
					} else if (flag=="3l_e") {
						lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"]+d2[,"e"], family = "binomial")
					} else if (flag=="360_a") {
						lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_60d"]+d2[,"a"], family = "binomial")
					} else if (flag=="360_b") {
						lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_60d"]+d2[,"b"], family = "binomial")
					} else if (flag=="360_c") {
						lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_60d"]+d2[,"c"], family = "binomial")
					} else if (flag=="360_d") {
						lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_60d"]+d2[,"d"], family = "binomial")
					} else if (flag=="360_e") {
						lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_60d"]+d2[,"e"], family = "binomial")
					} else if (flag=="360_abcde") {
						lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_60d"]+d2[,"a"]+d2[,"b"]+d2[,"c"]+d2[,"d"]+d2[,"e"], family = "binomial")
					} else if (flag=="3l_abcde") {
						lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"]+d2[,"a"]+d2[,"b"]+d2[,"c"]+d2[,"d"]+d2[,"e"], family = "binomial")
					} else if (flag=="3l_be") {
						lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"]+d2[,"b"]+d2[,"e"], family = "binomial")
					} else if (flag=="360_be") {
						lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_60d"]+d2[,"b"]+d2[,"e"], family = "binomial")
					} else if (flag=="360_n") {
						lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_60d"], family = "binomial")
					} else if (flag=="3l_n") {
						lmodel <- glm(y ~ x+d2[,"Age"]+d2[,"Race"]+d2[,"STD_life"], family = "binomial")
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
						a=data.frame(t(rep(NA,9)))
						colnames(a)=c("taxa","tp","ct","flag","rownames(a)","levelref1","OR","X2.5..","X97.5..")
					} else {
						a=a[rownames(a)%in%paste0("x",leveln),]
						if (nrow(a)==0) {
							print("glm error")
							exit(status=1)
						}
						a=cbind(taxa,tp,ct,flag,rownames(a),levelref1,a)
					}
					tab=rbind(tab,a)

					if (!exists("b")) {
						d=data.frame(t(rep(NA,10)))
						colnames(d)=c("taxa","tp","ct","flag","rownames(d)","levelref1",
							"Estimate","Std..Error","z.value","Pr...z..")
					} else {
						d=b[rownames(b)%in%paste0("x",leveln),]
						if (nrow(d)==0) {
							print("glm error")
							exit(status=1)
						}
						d=cbind(taxa,tp,ct,flag,rownames(d),levelref1,d)

						e=b[!(rownames(b)%in%paste0("x",leveln)),]
						e$taxa=taxa
						e$tp=tp
						e$ct=ct
						e$flag=flag
						e$ref=levelref1
					}
					tabb=rbind(tabb,d)
					tabe=rbind(tabe,e)

					leveln=leveln[-which(leveln==levelref1)]
				}
			}
		}

		if(nrow(tab)==0) {
			tab=data.frame(t(rep(NA,9)))
		}

		if(nrow(tabb)==0) {
			tabb=data.frame(t(rep(NA,10)))
		}
		
		if(nrow(tabe)==0) {
			tabe=data.frame(t(rep(NA,10)))
		}

		colnames(tab)=c("Taxa","Type","Cutoff","Cov","Test","Reference","OR","CI_2.5","CI_97.5")
		colnames(tabb)=c("Taxa","Type","Cutoff","Cov","Test","Reference",
			"Estimate","Std_err","z","p_val")
		colnames(tabe)=c("Estimate","Std_err","z","p_val",
			"Taxa","Type","Cutoff","Cov","Reference")

		if (length(which(is.na(tab[,1])))!=0) {
			tab=tab[-which(is.na(tab[,1])),]
		}
		if (length(which(is.na(tabb[,1])))!=0) {
			tabb=tabb[-which(is.na(tabb[,1])),]
		}
		if (length(which(is.na(tabe[,1])))!=0) {
			tabe=tabe[-which(is.na(tabe[,1])),]
		}

		rownames(tab)=paste(tab[,1],tab[,2],tab[,3],tab[,4],tab[,5],tab[,6],sep="_")
		rownames(tabb)=paste(tabb[,1],tabb[,2],tabb[,3],tabb[,4],tabb[,5],tabb[,6],sep="_")

		tt=merge(tab,tabb,by="row.names")

		tt=tt[,c("Taxa.x","Type.x","Cutoff.x","Cov.x","Test.x","Reference.x",
			"OR","CI_2.5","CI_97.5","Estimate","Std_err","z","p_val")]
		colnames(tt)[1:6]=c("Taxa","Type","Cutoff","Cov","Test","Reference")

		tt$Test=gsub("^x","",tt$Test)

		li = list(tt=tt,tabe=tabe)
		return(li)

	}
}


