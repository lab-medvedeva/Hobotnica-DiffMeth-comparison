library(optparse)
library(foreach)
library(doParallel)
library(Hobotnica)
library(genomation)
library(usedist)
library(Rankcluster)
library(fossil)

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
	trt<-rep(c('group1','group2'), times = c(n1,n2))
	samples_group1 <- samples[1:n1]
	samples_group2 <- samples[(n1+1):(n1+n2)]
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
	trt<-rep(c('group1','group2'), times = c(n1,n2))
	samples_group1 <- samples[1:n1]
	samples_group2 <- samples[(n1+1):(n1+n2)]
	assembly = "mm39"
	regions_file_name = 'cpgIslandExt.mm39.chr19.sample_300.bed'
}

type = "RRBS"
methyl_diffs = c(10,15,20,30)

get_H_score<-function(ratio_table, signature, trt) 
{
    sig_length<-length(signature)
    if (sig_length < 2)
    {
        return(list(H=0, pvalue=NA))
    }
    
    heat_matrix<-ratio_table[signature,]
    distMatrix<-dist(t(heat_matrix))
    if(sum(is.na(distMatrix))>0)
    {
        return(list(H=0, pvalue=NA))
    }

    H<-Hobotnica(distMatrix, trt)
    return(list(H=H, pvalue=NA))
}

for (simul in 0:9)
{
	setwd(paste0(dataset_dir, '/simulation', simul, '/result_alt/'))

	for (methyl_diff in methyl_diffs)
	{
		print(methyl_diff)
		dataset_name = paste(dname, methyl_diff, sep = '_')

		load(paste0(dataset_name,"_",methyl_diff,".rda"))
		ratio_table <- read.table(file = paste0("ratio_table_", methyl_diff), sep="\t")
		all <- GenomicRanges::granges(bsseq_stranded_filtered)
		cpg_islands<- readBed(paste0('../', regions_file_name))
		ratio_table_islands <- ratio_table[which( countOverlaps(all,  cpg_islands) %in% c(1)),]
		ground_truth <- rownames(ratio_table_islands)
		ground_truth_sig<-ratio_table[ground_truth, ]
		distMatrix<-dist(t(ground_truth_sig))
		Hobotnica(distMatrix, trt)

		not_ground_truth <- setdiff(rownames(ratio_table), ground_truth)
		not_ground_truth_sample <- sample(not_ground_truth, replace = FALSE)

		precs <-c()
		Hs <-c()
		pvalues <- c()
		for (i in 0:length(ground_truth))
		{
			prec <- (length(ground_truth) - i)/length(ground_truth)
			ground_truth[i] <- not_ground_truth_sample[i]
			ground_truth_sig <- ratio_table[ground_truth, ]
			distMatrix<-dist(t(ground_truth_sig))
			result <-  get_H_score(ratio_table, ground_truth, trt)
			H <- result$H
		    Hs <- c(Hs, H)
		    precs <- c(precs, prec)
		}

		df <- data.frame(precs, Hs)
		write.table(df, file = paste0("prec_H_", methyl_diff, "_", simul), sep="\t", col.names=TRUE, quote = FALSE)
	}
}
