# Genomospecies
### Author: Yue Xing (yxing4@luc.edu)

## For the 9 GS, make blast databases for them

	makeblastdb -in assemblies.fa -parse_seqids -title "assemblies" -dbtype nucl

## BLAST

	cat file.list | while read line
	do

		cat rmHm_fq/${line}_1.fq | sed -n '1~4s/^@/>/p;2~4p' - > ${line}_1.fasta
		cat rmHm_fq/${line}_2.fq | sed -n '1~4s/^@/>/p;2~4p' - > ${line}_2.fasta

		blastn -num_threads 10 \
		-db assemblies.fa \
		-query ${line}_1.fasta \
		-out ${line}_1.blastn.out \
		-max_hsps 1 -max_target_seqs 3 -evalue 1e-6 \
		-outfmt "6 qseqid sseqid sscinames scomnames sblastnames pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qcovus"

		blastn -num_threads 10 \
		-db assemblies.fa \
		-query ${line}_2.fasta \
		-out ${line}_2.blastn.out \
		-max_hsps 1 -max_target_seqs 3 -evalue 1e-6 \
		-outfmt "6 qseqid sseqid sscinames scomnames sblastnames pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qcovus"

		rm ${line}_1.fasta
		rm ${line}_2.fasta

	done

## Find perfect matches in blast results

Only perfect matches (%identity=100, e-value<1e-6, no mismatch, no gapopen, query coverage=100). And only used queries that align to exactly 1 subject.

	fls=list.files(".","_1.blastn.out")
	fl=fls[1]
	tab=read.delim(fl,header=FALSE)

	nms=data.frame(table(tab[,1]))
	nms=nms[nms[,2]==1,]
	nms=nms[,1]
	tab=tab[tab[,1]%in%nms,]

	tab2=tab[tab[,6]==100 & tab[,8]==0 & tab[,9]==0 & tab[,16]==100,]
	if (nrow(tab2)!=0) {
		tab2$id=gsub("_1\\.blastn\\.out","",fl)
		tab2=tab2[,c(19,1,2,6:18)]
	}
	
	for (i in 2:length(fls)) {
		d1=read.delim(fls[i],header=FALSE)

		nms=data.frame(table(d1[,1]))
		nms=nms[nms[,2]==1,]
		nms=nms[,1]
		d1=d1[d1[,1]%in%nms,]

		d1=d1[d1[,6]==100 & d1[,8]==0 & d1[,9]==0 & d1[,16]==100,]
		if (nrow(d1)!=0) {
			d1$id=gsub("_1\\.blastn\\.out","",fls[i])
			d1=d1[,c(19,1,2,6:18)]
			tab2=rbind(tab2,d1)
		}
		print(i)
	}

	colnames(tab2)=c("Sample","qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend",
	"sstart","send","evalue","bitscore","qcovs","qcovhsp","qcovus")
	
	tabL=tab2
	saveRDS(tabL,"tabL.rds")

	####

	fls=list.files(".","_2.blastn.out")
	fl=fls[1]
	tab=read.delim(fl,header=FALSE)

	nms=data.frame(table(tab[,1]))
	nms=nms[nms[,2]==1,]
	nms=nms[,1]
	tab=tab[tab[,1]%in%nms,]

	tab2=tab[tab[,6]==100 & tab[,8]==0 & tab[,9]==0 & tab[,16]==100,]
	if (nrow(tab2)!=0) {
		tab2$id=gsub("_2\\.blastn\\.out","",fl)
		tab2=tab2[,c(19,1,2,6:18)]
	}
	
	for (i in 2:length(fls)) {
		d1=read.delim(fls[i],header=FALSE)

		nms=data.frame(table(d1[,1]))
		nms=nms[nms[,2]==1,]
		nms=nms[,1]
		d1=d1[d1[,1]%in%nms,]

		d1=d1[d1[,6]==100 & d1[,8]==0 & d1[,9]==0 & d1[,16]==100,]
		if (nrow(d1)!=0) {
			d1$id=gsub("_2\\.blastn\\.out","",fls[i])
			d1=d1[,c(19,1,2,6:18)]
			tab2=rbind(tab2,d1)
		}
		print(i)
	}

	colnames(tab2)=c("Sample","qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend",
	"sstart","send","evalue","bitscore","qcovs","qcovhsp","qcovus")

	tabR=tab2
	saveRDS(tabR,"tabR.rds")

## Get the fq headers of the perfect, unique matches

	d3=readRDS("tabL.rds")
	nm=names(table(d3$Sample))
	for (n1 in nm) {
		tab=d3[d3$Sample==n1,]
		tab=data.frame(tab$qseqid)
		write.table(tab,paste0(n1,"_1.fqheaders.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
	}

	##

	d3=readRDS("tabR.rds")
	nm=names(table(d3$Sample))
	for (n1 in nm) {
		tab=d3[d3$Sample==n1,]
		tab=data.frame(tab$qseqid)
		write.table(tab,paste0(n1,"_2.fqheaders.txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
	}

####

	cat *_1.fqheaders.txt | wc -l
	cat *_2.fqheaders.txt | wc -l

## BLAST to nt database

	cat file.list | while read line
	do
	
		seqtk subseq ${line}_1.fasta IUMP_Gvags_class/${line}_1.fqheaders.txt > ${line}_1.sub.fa
		seqtk subseq ${line}_2.fasta IUMP_Gvags_class/${line}_2.fqheaders.txt > ${line}_2.sub.fa
		
		blastn -num_threads 10 \
		-db blast_db/nt \
		-query ${line}_1.sub.fa \
		-out ${line}_1.sub.blastn.out \
		-max_hsps 1 -max_target_seqs 3 -evalue 1e-6 \
		-outfmt "6 qseqid sseqid sscinames scomnames sblastnames pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qcovus"
		
		blastn -num_threads 10 \
		-db blast_db/nt \
		-query ${line}_2.sub.fa \
		-out ${line}_2.sub.blastn.out \
		-max_hsps 1 -max_target_seqs 3 -evalue 1e-6 \
		-outfmt "6 qseqid sseqid sscinames scomnames sblastnames pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp qcovus"

	done

## Get the perfect matches to the nt database

	fls=list.files(".","_1.sub.blastn.out")
	fl=fls[1]
	tab=read.delim(fl,header=FALSE)

	nms=data.frame(table(tab[,1]))
	nms=nms[nms[,2]==1,]
	nms=nms[,1]
	tab=tab[tab[,1]%in%nms,]

	tab2=tab[tab[,6]==100 & tab[,8]==0 & tab[,9]==0 & tab[,16]==100,]
	if (nrow(tab2)!=0) {
		tab2$id=gsub("_1\\.blastn\\.out","",fl)
		tab2=tab2[,c(19,1:18)]
	}
	
	for (i in 2:length(fls)) {
		if (file.info(fls[i])$size>0) {
			d1=read.delim(fls[i],header=FALSE)

			nms=data.frame(table(d1[,1]))
			nms=nms[nms[,2]==1,]
			nms=nms[,1]
			d1=d1[d1[,1]%in%nms,]

			d1=d1[d1[,6]==100 & d1[,8]==0 & d1[,9]==0 & d1[,16]==100,]
			if (nrow(d1)!=0) {
				d1$id=gsub("_1\\.blastn\\.out","",fls[i])
				d1=d1[,c(19,1:18)]
				tab2=rbind(tab2,d1)
			}
		}
		print(i)
	}

	colnames(tab2)=c("Sample","qseqid","sseqid","n1","n2","n3","pident","length",
	"mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qcovs","qcovhsp","qcovus")

	tabL=tab2
	saveRDS(tabL,"tabL_sub.rds")

####

	fls=list.files(".","_2.sub.blastn.out")
	fl=fls[1]
	tab=read.delim(fl,header=FALSE)

	nms=data.frame(table(tab[,1]))
	nms=nms[nms[,2]==1,]
	nms=nms[,1]
	tab=tab[tab[,1]%in%nms,]

	tab2=tab[tab[,6]==100 & tab[,8]==0 & tab[,9]==0 & tab[,16]==100,]
	if (nrow(tab2)!=0) {
		tab2$id=gsub("_2\\.blastn\\.out","",fl)
		tab2=tab2[,c(19,1:18)]
	}
	
	for (i in 2:length(fls)) {
		if (file.info(fls[i])$size>0) {
			d1=read.delim(fls[i],header=FALSE)

			nms=data.frame(table(d1[,1]))
			nms=nms[nms[,2]==1,]
			nms=nms[,1]
			d1=d1[d1[,1]%in%nms,]

			d1=d1[d1[,6]==100 & d1[,8]==0 & d1[,9]==0 & d1[,16]==100,]
			if (nrow(d1)!=0) {
				d1$id=gsub("_2\\.blastn\\.out","",fls[i])
				d1=d1[,c(19,1:18)]
				tab2=rbind(tab2,d1)
			}
		}
		print(i)
	}

	colnames(tab2)=c("Sample","qseqid","sseqid","n1","n2","n3","pident","length","mismatch",
	"gapopen","qstart","qend","sstart","send","evalue","bitscore","qcovs","qcovhsp","qcovus")

	tabR=tab2
	saveRDS(tabR,"tabR_sub.rds")

## Merge to a count table
	
	d1=readRDS("tabL_sub.rds")
	d2=readRDS("tabR_sub.rds")
	d3=rbind(d1,d2)
	
	d3$Sample=gsub("\\-.*","",d3$Sample)
	d3$Sample=gsub("\\..*","",d3$Sample)
	tab=data.frame(table(d3$Sample,d3$n1))
	tab=tab[tab$Freq!=0,]
	
	sp=read.table("Gvag_post.txt",header=FALSE)[,1]
	sid=names(table(tab$Var2))
	
	tab[,1]=as.character(tab[,1])
	tab[,2]=as.character(tab[,2])
	
	tab=tab[tab[,1]%in%sp,]

	df=data.frame(matrix(ncol=length(sp),nrow=length(sid)))
	colnames(df)=sp
	rownames(df)=sid
	
	for (i in 1:nrow(tab)) {
		df[tab[i,2],tab[i,1]]=tab[i,3]
	}
	df[is.na(df)]=0

	write.csv(df,"Gvags_sub_blast.csv",quote=FALSE)

## Extract the non G vags reads

	s=c("Gardnerella vaginalis","Gardnerella vaginalis 409-05","Gardnerella vaginalis HMP9231")
	d3=d3[!(d3$n1%in%s),]
	
	saveRDS(d3,"noiseSeqs.rds")
	
	dn=d3

## Remove the above reads from G vag blast results

	d1=readRDS("tabL.rds")
	d2=readRDS("tabR.rds")
	d3=rbind(d1,d2)
	d3$Sample=gsub("\\-.*","",d3$Sample)
	d3$Sample=gsub("\\..*","",d3$Sample)

	d3=d3[!(d3$qseqid%in%dn$qseqid),]

	tab=data.frame(table(d3$Sample,d3$sseqid))
	tab=tab[tab$Freq!=0,]
	
	sp=read.table("sublist/Gvag_post.txt",header=FALSE)[,1]
	sid=names(table(tab$Var2))
	
	tab[,1]=as.character(tab[,1])
	tab[,2]=as.character(tab[,2])
	
	tab=tab[tab[,1]%in%sp,]

	df=data.frame(matrix(ncol=length(sp),nrow=length(sid)))
	colnames(df)=sp
	rownames(df)=sid
	
	for (i in 1:nrow(tab)) {
		df[tab[i,2],tab[i,1]]=tab[i,3]
	}
	df[is.na(df)]=0

	write.csv(df,"Gvags_blast_filtered.csv",quote=FALSE)

## Match the counts to the 9 species

	d1=read.csv("Gvags_blast_filtered.csv",row.names=1)
	d2=read.delim("GS.list",header=FALSE,row.names=1)
	rownames(d1)=gsub("ref\\|","",rownames(d1))
	rownames(d1)=gsub("\\|","",rownames(d1))

	d3=merge(d2,d1,all.y=TRUE,by="row.names")

	d3=d3[,-1]

	library(dplyr)
	s=d3 %>%
	    group_by(V2) %>%
	    summarise_all(funs(sum)) %>%
	    as.data.frame()
	colnames(s)[1]="ID"
	write.csv(s,"Gvags_blast_filtered_counts.csv",quote=FALSE,row.names=FALSE)

## Find proportions

	s=read.csv("Gvags_blast_filtered_counts.csv",row.names=1)
	for (i in 1:ncol(s)) {
		s[,i]=s[,i]/sum(s[,i])*100
	}
	colSums(s)
	s[is.na(s)]=0
	s=round(s,4)
	write.csv(s,"Gvags_blast_filtered_proportions.csv",quote=FALSE)

## Plot

	library(circlize)
	library(ComplexHeatmap)

	d1=read.csv("Gvags_blast_filtered_proportions.csv",row.names=1)
	#colnames(d1)=gsub("\\..*","",colnames(d1))

	# Gvag_pos.csv: G. vags positive samples from MP3 output
	nms=read.csv("Gvag_pos.csv")
	all(nms[,1]==colnames(d1))

	ut=read.csv("UTs.csv")
	ut=ut[!is.na(ut[,2]),]
	ut=ut[ut[,1]%in%nms[,1],]
	rownames(ut)=ut[,1]
	ut=ut[colnames(d1),]
	print(all(rownames(ut)==colnames(d1)))

	col_fun = colorRamp2(c(0, 0.01, 0.1, 1, 100), 
		c("blue", "skyblue", "yellow", "orange", "red"))

	col_list=list(
		Urethrotype=c("UT1"="steelblue1","UT2"="deeppink1","3"="purple2","4"="khaki1","5"="orange","6"="black","7"="grey")
	)

	column_taxa = HeatmapAnnotation(
		Urethrotype = ut$UT,
		col = col_list)

	png(paste0("Gvags_blast_filtered.png"),width=5000,height=2500,res=300)

	ht = Heatmap(d1, name = "Abundance",
		top_annotation = column_taxa,
	    #column_names_rot = 45,
	    row_names_gp = gpar(fontsize = 10),
	    #row_names_side = "left",
	   	column_names_gp = gpar(fontsize = 10),
	   	row_title_gp = gpar(fontsize = 12),
	   	column_title_gp = gpar(fontsize = 12),
	    column_title = "Subjects", 
	    column_names_side = "top",
	    row_title = "Genomospecies",
	    #row_dend_width = unit(rdw, "mm"),
	    #column_dend_height = unit(cdw, "mm"),
	    row_dend_side = "right",
	    col=col_fun
	)
	ht=draw(ht)
	print(ht)

	dev.off()

## Make tables

	for (nm in c("UT1","UT2")) {
		nms=read.csv("Gvag_pos.csv")
	
		ut=read.csv("UTs.csv")
		rownames(ut)=ut[,1]
		ut=ut[ut[,1]%in%nms[,1],]
		ut=ut[(!(is.na(ut[,2])) & ut[,2]==nm),]
	
		d1=read.csv("Gvags_blast_counts.csv",row.names=1)
		colnames(d1)=gsub("\\..*","",colnames(d1))
		d1=d1[,ut[,1]]
	
		write.csv(d1,paste0("Gvags_blast_filtered_counts_",nm,".csv"),quote=FALSE)

		s=d1
		for (i in 1:ncol(s)) {
			s[,i]=s[,i]/sum(s[,i])*100
		}
		print(colSums(s))
		s[is.na(s)]=0
		s=round(s,4)
		write.csv(s,paste0("Gvags_blast_filtered_proportions_",nm,".csv"),quote=FALSE)

	}



