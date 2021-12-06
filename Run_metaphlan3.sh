# Run Metaphlan3 on the sequence data
# Author: Yue Xing (yxing4@luc.edu)

# Run Metaphlan3

	#!/bin/bash
	source activate mpa3
	module load perl/5.30.1
	metaphlan -v

	cat file.list | while read line
	do

		NAME=${line}
		NAME2=$(echo ${NAME} | cut -d"/" -f9-)
		FQ1=$(ls ${NAME}*R1*.fastq.gz)
		FQ2=$(ls ${NAME}*R2*.fastq.gz)

		echo $FQ1
		echo $FQ2
		echo $NAME2

		metaphlan ${FQ1},${FQ2} --input_type fastq --nproc 8 \
		-t rel_ab_w_read_stats --add_viruses \
		--bowtie2out ${dir}/metaphlan3/${NAME2}.bowtie2.bz2 \
		-o ${dir}/metaphlan3/${NAME2}_profile.txt

	done

# Replace ' and # in taxon names

	a=$(ls *_profile.txt)
	for file in $a
	do
		sed "s/'/_/g" ${file} > ${file}.2
		sed "s/#/_/g" ${file}.2 > ${file}.3
		rm ${file}.2
		mv ${file}.3 ${file}.2
	done

# Process results and merge results

	R

	source("~/scripts/process_metaphlan3_wct.R")
	process_metaphlan_all("metaphlan3")

	MDIR=~/tools/metaphlan2
	dir=metaphlan3

	${MDIR}/utils/merge_metaphlan_tables.py ${dir}/count/*_profile_count.txt > ${dir}/Metaphlan3_count_table.txt

	${MDIR}/utils/merge_metaphlan_tables.py ${dir}/ra/*_profile_ra.txt > ${dir}/Metaphlan3_ra_table.txt

	## Add QC data

	d1=read.table("Assigned_taxa.txt")
	d2=read.delim("multiqc_fastqc.txt")
	d2=d2[,c(1,5)]
	d2[,1]=gsub("_S._L..._R._001","",d2[,1])
	d2[,1]=gsub("_S.._L..._R._001","",d2[,1])
	d3=merge(d2,d1,all.x=TRUE,by.x="Sample",by.y="V1")
	colnames(d3)[2:3]=c("Total_reads","Assigned_reads")
	d3[,1]=gsub("\\-","\\.",d3[,1])

	d4=read.csv("human.counts.csv",row.names=1)
	d5=merge(d3,d4,all.x=TRUE,by="Sample")
	d5=d5[,c(1,2,4,3)]
	write.csv(d5,"QC.csv",quote=FALSE,row.names=FALSE)

	##

	fls=list.files(".","csv")

	for (fl in fls) {
		d1=read.csv(fl,row.names=1)
		d1=t(d1)
		dn=merge(d5,d1,by="row.names",all.x=TRUE)
		rownames(dn)=dn[,1]
		dn=dn[,-1]
		dn=t(dn)
		fl=gsub("\\.csv","",fl)
		write.csv(dn,paste0(fl,".QC.csv"),quote=FALSE)
	}







