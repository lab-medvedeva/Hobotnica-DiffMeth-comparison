library(genomation)
library(usedist)
library(Rankcluster)
library(fossil)
library(foreach)
library(doParallel)
library(Hobotnica)
library(optparse)

option_list = list(make_option(c("-d", "--dataset_dir"), type="character", default=NULL,
        help="Working directory", metavar="character"))
    make_option(c("-g", "--group"), type="character", default=1,
        help="Group: 1 or 2"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
dataset_dir <- opt$dataset_dir
group <- opt$group

if(group == 1)
{
	dname <- 'default'
	samples <- paste0('sample', 0:31) 
	n1=16
	n2=16
	assembly = "hg38"
	files <- c(paste0(0:15, '_base.cov'), paste0(16:31, '_alt_', methyl_diff, '.sorted.cov'))
	regions_file_name = 'cpgIslandExt.hg38.chr22.sample_300.bed'
}
if(group == 2)
{
	dname <- 'GSE103886'
	samples <- c('SRR6040066', 'SRR6040067', 'SRR6040068', 'SRR6040069', 'SRR6040070', 'SRR6040071', 'SRR6040072', 'SRR6040073',
		'SRR6040074', 'SRR6040075', 'SRR6040076', 'SRR6040077', 'SRR6040078', 'SRR6040079', 'SRR6040080', 'SRR6040081', 'SRR6040082',
		'SRR6040083', 'SRR6040084', 'SRR6040085', 'SRR6040086', 'SRR6040087', 'SRR6040088')
	n1=11
	n2=12
	assembly = "mm39"
	files <- c(paste0(samples[1:n1], '_base.cov'),  paste0(samples[(n1+1):(n1+n2)], '_alt_', methyl_diff, '.sorted.cov'))
	regions_file_name = 'cpgIslandExt.mm39.chr19.sample_300.bed'

}

type = "RRBS"
methyl_diffs = c(10,15,20,30)
for (simul in 0:9)
{
	setwd(paste0(dataset_dir, '/simulation', simul, '/result_perm/'))
	for (methyl_diff in methyl_diffs)
	{
		dataset_name = paste(dname, methyl_diff, sep = '_')
		load(paste0(dataset_name,"_",methyl_diff,".rda"))
		ratio_table <- read.table(file = paste0("ratio_table_", methyl_diff), sep="\t")
		all <- GenomicRanges::granges(bsseq_stranded_filtered)
		cpg_islands<- readBed(paste0('../', regions_file_name))
		ratio_table_islands <- ratio_table[which( countOverlaps(all,  cpg_islands) %in% c(1)),]
		ground_truth <- rownames(ratio_table_islands)
		FPR_df <- data.frame(dataset=character(),
			method=character(),
			FPR=double(),
			stringsAsFactors=FALSE)

		method_name <- "radmeth"
		signature_csv <- read.csv(paste(dataset_name, method_name, "full", methyl_diff, sep = "_"), sep = '\t')
		if(dim(signature_csv)[1]>0)
		{
			signature<-paste0(signature_csv$V1,":",signature_csv$V2)
			false_positives <- signature
			FPR = length(false_positives)/length(rownames(ratio_table))
			
		}
		else
		{
			FPR = 0
		}

		new_row<-data.frame(method_name, FPR)
		names(new_row)<-c("method","FPR")
		FPR_df<-rbind(FPR_df, new_row)

		method_name <- "methylsig"
		signature_csv <- read.csv(paste(dataset_name, method_name, "full", methyl_diff, sep = "_"), sep = '\t')
		if(dim(signature_csv)[1]>0)
		{
			signature<-paste0(signature_csv$seqnames,":",signature_csv$start)
			false_positives <- signature
			FPR = length(false_positives)/length(rownames(ratio_table))
		}
		else
		{
			FPR = 0
		}

		new_row<-data.frame(method_name, FPR)
		names(new_row)<-c("method","FPR")
		FPR_df<-rbind(FPR_df, new_row)

		method_name <- "hmm"
		signature_csv <- read.csv(paste(dataset_name, method_name, "full", methyl_diff, sep = "_"), sep = '\t')
		if(dim(signature_csv)[1]>0)
		{
			signature<-paste0(signature_csv$chr,":",signature_csv$pos)
			false_positives <- signature
			FPR = length(false_positives)/length(rownames(ratio_table))
		}
		else
		{
			FPR = 0
		}

		new_row<-data.frame(method_name, FPR)
		names(new_row)<-c("method","FPR")
		FPR_df<-rbind(FPR_df, new_row)

		method_names <- c("methylkit", "methylkit_overdispersion")

		for (method_name in method_names)
		{
			signature_csv <- read.csv(paste(dataset_name, method_name, "full", methyl_diff, sep = "_"), sep = '\t')
			if(dim(signature_csv)[1]>0)
			{
				signature<-paste0(signature_csv$chr,":",signature_csv$start)
				false_positives <- signature
				FPR = length(false_positives)/length(rownames(ratio_table))
			}
			else
			{
				FPR = 0
			}

			new_row<-data.frame(method_name, FPR)
			names(new_row)<-c("method","FPR")
			FPR_df<-rbind(FPR_df, new_row)
		}

		method_names <- c("dss_no_smoothing", "dss_smoothing")
		for (method_name in method_names)
		{
			signature_csv <- read.csv(paste(dataset_name, method_name, "full", methyl_diff, sep = "_"), sep = '\t')
			if(dim(signature_csv)[1]>0)
			{
		        signature<-paste0(signature_csv$chr,":",signature_csv$pos)
				false_positives <- signature
				FPR = length(false_positives)/length(rownames(ratio_table))

			}
			else
			{
				FPR = 0
			}

			new_row<-data.frame(method_name, FPR)
			names(new_row)<-c("method","FPR")
			FPR_df<-rbind(FPR_df, new_row)
		}
		write.table(FPR_df, file = paste0("FPR_full_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)
	}
}
