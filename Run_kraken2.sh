# Run kraken2 on the sequence data
# Author: Yue Xing (yxing4@luc.edu)

# Remove human reads from fastqs

	#!/bin/bash
	cat file.list | while read line
	do

		NAME=${line}
		FQ1=${NAME}*R1*.fastq.gz
		FQ2=${NAME}*R2*.fastq.gz

		echo $FQ1
		echo $FQ2
		echo $NAME2

		kraken2 --paired --threads 10 --minimum-base-quality 30 \
		--db kraken_human_db --report ${NAME}.kraken2.report \
		--unclassified-out rmHm_fq/${NAME}#.fq ${FQ1} ${FQ2} > ${NAME}.kraken

	done

# Run kraken2 on human reads removed fastqs

	db_name=kraken_db
	out_dir=out
	in_fq_dir=rmHm_fq

	cat name.list | while read line
	do

		NAME=${line}
		echo $NAME

		kraken2 --paired --threads 10 --minimum-base-quality 30 \
		--db ${db_name} --output ${out_dir}/kraken/${NAME}.kraken2.txt \
		--use-mpa-style --report ${out_dir}/kraken/${NAME}.kraken2_mpa.txt \
		${in_fq_dir}/${NAME}_1.fq ${in_fq_dir}/${NAME}_2.fq

	done

# Process results

	a=$(ls *kraken2_mpa.txt)
	for file in $a
	do
		sed "s/'/_/g" ${file} > ${file}.2
		sed "s/#/_/g" ${file}.2 > ${file}.3
		sed "s/d__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Propionibacteriales|f__Propionibacteriaceae|g__Ponticoccus/d__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Propionibacteriales|f__Propionibacteriaceae|g__Ponticoccus2/g" ${file}.3 > ${file}.4
		rm ${file}.2
		rm ${file}.3
		mv ${file}.4 ${file}.2
	done

	##

	MDIR=MetaPhlAn2-master

	${MDIR}/utils/merge_metaphlan_tables.py ${dir}/*.kraken2_mpa.txt > ${dir}/kraken_raw_merged_abundance_table.txt

