# Author: Yue Xing (yxing4@luc.edu)

process_metaphlan3 = function(fl) {
	a=character()
	for (name in fl) {
		print(name)
		d1=read.delim(paste0(name,"_profile.txt.2"),skip=4)
		#d1=as.matrix(d1)

		tot=read.table(paste0(name,"_profile.txt.2"),sep=":",skip=3,nrows=1)
		if (tot[1,1]!="_estimated_reads_mapped_to_known_clades") {
			print("error")
			quit(status=1)
		}
		tot=tot[1,2]
		a=rbind(a,c(name,tot))

		dr=d1[,c(1,3)]
		dc=d1[,c(1,5)]
		colnames(dr)[1:2]=c("#SampleID","Metaphlan_Analysis")
		colnames(dc)[1:2]=c("#SampleID","Metaphlan_Analysis")
		write.table(dr,paste0("ra/",name,"_profile_ra.txt"),sep="\t",quote=FALSE,row.names=FALSE)
		write.table(dc,paste0("count/",name,"_profile_count.txt"),sep="\t",quote=FALSE,row.names=FALSE)

	}
	write.table(a,"Assigned_taxa.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
}


##

process_metaphlan_all = function(in_dir) {

	setwd(in_dir)

	options(stringsAsFactors=FALSE)
			
	fl=list.files(".","*_profile.txt.2")
	fl=gsub("_profile\\.txt\\.2","",fl)
	print(fl)
	
	process_metaphlan3(fl)
}

