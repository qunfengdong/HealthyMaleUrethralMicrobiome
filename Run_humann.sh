# Run humann
# Author: Yue Xing (yxing4@luc.edu)

cat file.list | while read line
do

	a=$(ls all/${line}_*_R1_001.fastq.gz)
	b=$(ls all/${line}_*_R2_001.fastq.gz)
	cat ${a} ${b} > ${line}.fastq.gz
	humann --input ${line}.fastq.gz --output ${line}

done

humann_join_tables -i out -o healthy_genefamilies.tsv --file_name genefamilies -s

# remove "unmapped"
sed -i '2d' healthy_genefamilies.tsv

# renorm by ra
humann_renorm_table -i healthy_genefamilies.tsv -u relab -p -o healthy_genefamilies.ra.tsv

# uniref90_level4ec
humann_regroup_table -i healthy_genefamilies.ra.tsv -o healthy_genefamilies.ra.ec.tsv -g uniref90_level4ec -e 6
humann_rename_table -i healthy_genefamilies.ra.ec.tsv -o healthy_genefamilies.ra.ec.ec.tsv -n ec
# uniref90_ko
humann_regroup_table -i healthy_genefamilies.ra.tsv -o healthy_genefamilies.ra.ko.tsv -g uniref90_ko -e 6
humann_rename_table -i healthy_genefamilies.ra.ko.tsv -o healthy_genefamilies.ra.ko.kegg_orthology.tsv -n kegg-orthology
# uniref90_go
humann_regroup_table -i healthy_genefamilies.ra.tsv -o healthy_genefamilies.ra.go.tsv -g uniref90_go -e 6
humann_rename_table -i healthy_genefamilies.ra.go.tsv -o healthy_genefamilies.ra.go.go.tsv -n go
# metacyc-rxn
humann_regroup_table -i healthy_genefamilies.ra.tsv -o healthy_genefamilies.ra.rxn.tsv -g uniref90_rxn -e 6
humann_rename_table -i healthy_genefamilies.ra.rxn.tsv -o healthy_genefamilies.ra.rxn.metacyc_rxn.tsv -n metacyc-rxn

# renorm by cpm	
humann_renorm_table -i healthy_genefamilies.tsv -u cpm -p -o healthy_genefamilies.cpm.tsv

# uniref90_level4ec
humann_regroup_table -i healthy_genefamilies.cpm.tsv -o healthy_genefamilies.cpm.ec.tsv -g uniref90_level4ec
humann_rename_table -i healthy_genefamilies.cpm.ec.tsv -o healthy_genefamilies.cpm.ec.ec.tsv -n ec
# uniref90_ko
humann_regroup_table -i healthy_genefamilies.cpm.tsv -o healthy_genefamilies.cpm.ko.tsv -g uniref90_ko
humann_rename_table -i healthy_genefamilies.cpm.ko.tsv -o healthy_genefamilies.cpm.ko.kegg_orthology.tsv -n kegg-orthology
# uniref90_go
humann_regroup_table -i healthy_genefamilies.cpm.tsv -o healthy_genefamilies.cpm.go.tsv -g uniref90_go
humann_rename_table -i healthy_genefamilies.cpm.go.tsv -o healthy_genefamilies.cpm.go.go.tsv -n go
# metacyc-rxn
humann_regroup_table -i healthy_genefamilies.cpm.tsv -o healthy_genefamilies.cpm.rxn.tsv -g uniref90_rxn
humann_rename_table -i healthy_genefamilies.cpm.rxn.tsv -o healthy_genefamilies.cpm.rxn.metacyc_rxn.tsv -n metacyc-rxn

