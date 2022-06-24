#!/usr/bin/env Rscript
library(optparse)
library(bsseq)
library(BiocParallel)
library(GenomicFeatures)

option_list = list(
    make_option(c("-d", "--dataset"), type="character", default=NULL,
        help="Dataset GEO id", metavar="character"),
    make_option(c("-w", "--working_dir"), type="character", default=NULL,
        help="Working directory", metavar="character"),
    make_option(c("-c", "--cores"), type="integer", default=16,
        help="Number of cores"),
    make_option(c("-s", "--source_dir"), type="character", default=NULL,
        help="Source directory", metavar="character"))
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$dataset)){
    print_help(opt_parser)
    stop("Dataset id should be supplied", call.=FALSE)
}

if (is.null(opt$working_dir)){
    print_help(opt_parser)
    stop("Working directory should be supplied", call.=FALSE)
}

if (is.null(opt$source_dir)){
    print_help(opt_parser)
    stop("Source directory should be supplied", call.=FALSE)
}

# Set working directory and source directory
setwd(opt$working_dir)
source(file.path(opt$source_dir,"get_score.R"))

# Dataset initialization
dataset_name<-opt$dataset
switch(dataset_name, 
GSE149608={
    n1=10
    n2=10
    type = "WGBS"
    samples<-c('SRR11647648', 'SRR11647649', 'SRR11647650', 'SRR11647651', 'SRR11647652', 'SRR11647653',
        'SRR11647654', 'SRR11647655', 'SRR11647656', 'SRR11647657', 'SRR11647658', 'SRR11647659', 'SRR11647660',
        'SRR11647661', 'SRR11647662', 'SRR11647663', 'SRR11647664', 'SRR11647665', 'SRR11647666', 'SRR11647667')
},
GSE138598={
    n1=8
    n2=9
    type = "WGBS"
    samples <- c('T2D_1', 'T2D_2', 'T2D_3', 'T2D_4', 'T2D_5', 'T2D_6', 'T2D_7', 'T2D_8',
        'Control_1', 'Control_2', 'Control_3', 'Control_4', 'Control_5', 'Control_6', 'Control_7', 
        'Control_8', 'Control_9')
},
GSE119980={
    n1=6
    n2=6
    type = "WGBS"
    samples <-c('SRR7830270', 'SRR7830271', 'SRR7830272', 'SRR7830273', 'SRR7830274', 'SRR7830275', 'SRR7830276',
        'SRR7830277', 'SRR7830278', 'SRR7830279', 'SRR7830280', 'SRR7830281')
},
GSE117593={
    n1=25
    n2=25
    type = "WGBS"
    samples<-c('SRR7587054', 'SRR7587055', 'SRR7587056', 'SRR7587057', 'SRR7587058', 'SRR7587059', 'SRR7587060',
        'SRR7587061', 'SRR7587062', 'SRR7587063', 'SRR7587064', 'SRR7587065', 'SRR7587066', 'SRR7587067', 'SRR7587068',
        'SRR7587069', 'SRR7587070', 'SRR7587071', 'SRR7587072', 'SRR7587073', 'SRR7587074', 'SRR7587075', 'SRR7587076',
        'SRR7587077', 'SRR7587078', 'SRR7587079', 'SRR7587080', 'SRR7587081', 'SRR7587082', 'SRR7587083', 'SRR7587084',
        'SRR7587085', 'SRR7587086', 'SRR7587087', 'SRR7587088', 'SRR7587089', 'SRR7587090', 'SRR7587091', 'SRR7587092',
        'SRR7587093', 'SRR7587094', 'SRR7587095', 'SRR7587096', 'SRR7587097', 'SRR7587098', 'SRR7587099', 'SRR7587100',
        'SRR7587101', 'SRR7587102', 'SRR7587103')
},
{
    stop("Unknown dataset ID", call.=FALSE)
}
)

files<-paste0(samples, ".cov")
samples_group1<-samples[1:n1]
samples_group2<-samples[(n1+1):(n1+n2)]
trt<-rep(c('group1','group2'), times = c(n1,n2))

# Load filtered data
load(paste0(dataset_name,".rda"))

# Smooth initial dataset
BiocParallel::register(BiocParallel::MulticoreParam(workers = opt$cores, progressbar = TRUE))

bsseq_stranded = bsseq::read.bismark(
    files = files,
    colData = data.frame(row.names = samples),
    rmZeroCov = TRUE,
    strandCollapse = TRUE
)

bsseq_stranded<-sort(bsseq_stranded)
bsseq_stranded_smooth<-BSmooth(bsseq_stranded,  verbose = FALSE)

# Filter smoothed data by coverage
pos <- GenomicRanges::granges(bsseq_stranded_filtered)
bsseq_stranded_smooth_filtered<-subsetByOverlaps(bsseq_stranded_smooth, pos)
compl_cases<-bsseq_stranded_smooth_filtered[complete.cases(getMeth(bsseq_stranded_smooth_filtered)),]

# Save smoothed dataset
save(compl_cases, file = paste0(dataset_name,".smoothed.rda"))
# Save number of CpG sites for test
method_name<-"bsmooth"
cpg_num_table<-read.table(file = "cpg_num_table", sep="\t")
cpg_num<-dim(compl_cases)[1]
new_row<-data.frame(dataset_name, method_name, cpg_num)
names(new_row)<-c("dataset","method","cpg_num")
cpg_num_table<-rbind(cpg_num_table, new_row)
write.table(cpg_num_table, file = "cpg_num_table", sep="\t", col.names=TRUE, quote = FALSE)

# Get DMC signature
tstat<-BSmooth.tstat(compl_cases,
  group1=samples_group1,
  group2=samples_group2,
  estimate.var = "group2",
  local.correct = TRUE,
  verbose = TRUE)

tstats_df<-data.frame(getStats(tstat))
gr<-data.frame(tstat@gr)
rownames(tstats_df)<-paste0(gr$seqnames,":",gr$start)
q<-0.05
qcutoff = c(0+q, 1-q)
cutoffs<-quantile(tstats_df$tstat.corrected, qcutoff, na.rm = TRUE)
dmls<-tstats_df[tstats_df$tstat.corrected <= cutoffs[1] | tstats_df$tstat.corrected >= cutoffs[2],]
dmls_filtered<-dmls[!is.na(dmls$tstat.corrected),]
dmls_sorted<-dmls_filtered[order(-abs(dmls_filtered$tstat.corrected)),]
write.table(dmls_sorted, paste(dataset_name, method_name, "full", sep = "_"), sep = '\t')

# Filter by mean methylation difference (> 15%)
dmls_diff_cut<-dmls_sorted[abs(dmls_sorted$group2.means - dmls_sorted$group1.means) >= 0.15,]
write.table(dmls_diff_cut, paste(dataset_name, method_name, "cut", sep = "_"), sep = '\t')
dmls_diff_cut <-read.csv( paste(dataset_name, method_name, "cut", sep = "_"), sep ='\t')
# Get Hobotnica score
signature<-rownames(dmls_diff_cut)
sig_length<-length(signature)
result<-get_H_score(dataset_name, method_name, signature, trt, opt$cores)

# Write H score to result table
H_full<-read.table("H_full", header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, type, result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_full<-rbind(H_full, new_row)
write.table(H_full, file = "H_full", sep="\t", col.names=TRUE, quote = FALSE)

# Get Hobotnica score for top 100 signature
signature_top_100<-head(signature, 100)
sig_length<-length(signature_top_100)
result <-get_H_score(dataset_name, method_name, signature_top_100, trt, opt$cores) 

# Write H score to result table
H_100<-read.table("H_100", header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, type, result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_100<-rbind(H_100, new_row)
write.table(H_100, file = "H_100", sep="\t", col.names=TRUE, quote = FALSE)