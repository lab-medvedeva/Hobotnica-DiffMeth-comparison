#!/usr/bin/env Rscript
library(optparse)
library(methylSig)
library(bsseq)

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
    trt<-rep(c('group1','group2'), times = c(n1,n2))

},
GSE138598={
    n1=8
    n2=9
    type = "WGBS"
    trt<-rep(c('group1','group2'), times = c(n1,n2))

},
GSE119980={
    n1=6
    n2=6
    type = "WGBS"
    trt<-rep(c('group1','group2'), times = c(n1,n2))

},
GSE148060={
    n1=21
    n2=32
    type = "RRBS"
    trt<-rep(c('group1','group2'), times = c(n1,n2))

},
GSE103886={
    n1=11
    n2=12
    type = "RRBS"
    trt<-rep(c('group1','group2'), times = c(n1,n2))

},
 GSE150592=
 {
    samples <- c('SRR11790875','SRR11790876','SRR11790877','SRR11790878','SRR11790879','SRR11790880','SRR11790881','SRR11790882','SRR11790883','SRR11790884','SRR11790885','SRR11790886','SRR11790887','SRR11790888','SRR11790889','SRR11790890','SRR11790891','SRR11790892','SRR11790893','SRR11790894','SRR11790895','SRR11790896','SRR11790897','SRR11790898','SRR11790899','SRR11790900','SRR11790901','SRR11790902','SRR11790903','SRR11790904')
    trt<-c('group2','group1','group2','group1','group2','group1','group2','group1','group2','group1','group2','group1','group1','group2','group2','group1','group1','group2','group2','group1','group1','group2','group2','group1','group2','group1','group2','group1','group1','group2')
    type = "RRBS"

},
{
    stop("Unknown dataset ID", call.=FALSE)
}
)

 method_name <- 'methylsig'
# Load filtered data
load(paste0(dataset_name,".rda"))

# Get DMC signature
diff_gr<-diff_methylsig(
 bs = bsseq_stranded_filtered,
 group_column = 'Condition',
 comparison_groups = c('control' = 'group1', 'case' = 'group2'),
 disp_groups = c('case' = TRUE, 'control' = TRUE),
 local_window_size = 0,
 t_approx = TRUE,
 n_cores = opt$cores)

dmls_fdr<-diff_gr[!is.na(mcols(diff_gr)$fdr),]
dmls_fdr<-dmls_fdr[mcols(dmls_fdr)$fdr <= 0.05,]
dmls_fdr_sorted<-dmls_fdr[order(dmls_fdr$fdr),]
write.table(dmls_fdr_sorted, paste0(dataset_name,"_methylsig_full"), sep = '\t')

dmls_diff_sorted<-dmls_fdr_sorted[order(-abs(dmls_fdr_sorted$meth_diff)),]
write.table(dmls_diff_sorted, paste(dataset_name, method_name, "meth_diff_full", sep = "_"), sep = '\t')

# Filter by mean methylation difference (> 15%)
#dmls_fdr_sorted <- read.csv(paste0(dataset_name,"_methylsig_full"), sep = '\t')
dmls_fdr_diff_cut<-dmls_fdr_sorted[abs(dmls_fdr_sorted$meth_diff) >= 15,]
write.table(dmls_fdr_diff_cut, paste0(dataset_name,"_methylsig_cut"), sep = '\t')

# Get Hobotnica score
dmls_fdr_diff_cut <- read.csv(paste0(dataset_name,"_methylsig_cut"), sep = '\t')

dmls<-data.frame(dmls_fdr_diff_cut)
if(dim(dmls)[1] == 0)
{
    signature <- c()
} else {
    signature<-paste0(dmls$seqnames,":",dmls$start)
}

sig_length<-length(signature)
result <-get_H_score(dataset_name, "methylsig", signature, trt, opt$cores) 

# Write H score to result table
H_full<-read.table("H_full", header=TRUE, sep='\t')
new_row<-data.frame(dataset_name,"methylsig", result$H, type, result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_full<-rbind(H_full, new_row)
write.table(H_full, file = "H_full", sep="\t", col.names=TRUE, quote = FALSE)

# Get Hobotnica score for top 100 signature
signature_top_100<-head(signature, 100)
sig_length<-length(signature_top_100)
result <-get_H_score(dataset_name, "methylsig", signature_top_100, trt, opt$cores) 

# Write H score to result table
H_100<-read.table("H_100", header=TRUE, sep='\t')
new_row<-data.frame(dataset_name,"methylsig", result$H, type, result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_100<-rbind(H_100, new_row)
write.table(H_100, file = "H_100", sep="\t", col.names=TRUE, quote = FALSE)

# Get Hobotnica score for top 10 signature
signature_top_10<-head(signature, 10)
sig_length<-length(signature_top_10)
result <-get_H_score(dataset_name, "methylsig", signature_top_10, trt, opt$cores) 

# Write H score to result table
H_10<-read.table("H_10", header=TRUE, sep='\t')
new_row<-data.frame(dataset_name,"methylsig", result$H, type, result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_10<-rbind(H_10, new_row)
write.table(H_10, file = "H_10", sep="\t", col.names=TRUE, quote = FALSE)