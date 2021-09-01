#!/usr/bin/env Rscript
library(optparse)
library(bsseq)
library(stringr)
library(parallel)

option_list = list(
    make_option(c("-d", "--dataset"), type="character", default=NULL,
        help="Dataset GEO id", metavar="character"),
    make_option(c("-w", "--working_dir"), type="character", default=NULL,
        help="Working directory", metavar="character"),
    make_option(c("-c", "--cores"), type="integer", default=16,
        help="Number of cores"),
    make_option(c("-s", "--source_dir"), type="character", default=NULL,
        help="Source directory", metavar="character"),
    make_option(c("-k", "--hmm_dir"), type="character", default=NULL,
        help="HMM-DMM source directory", metavar="character"));

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

if (is.null(opt$hmm_dir)){
    print_help(opt_parser)
    stop("HMM-DM source directory should be supplied", call.=FALSE)
}

# Set working directory and source directory
setwd(opt$working_dir)
source(file.path(opt$source_dir,"get_score.R"))
source(file.path(opt$hmm_dir,"HMM.DM.R"))

# Dataset initialization
dataset_name<-opt$dataset
switch(dataset_name, 
GSE149608={
    n1=10
    n2=10
    type = "WGBS"
    chrs<-c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', 
      '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')
},
GSE138598={
    n1=8
    n2=9
    type = "WGBS"
    chrs<-c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', 
      '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')
},
GSE119980={
    n1=6
    n2=6
    type = "WGBS"
    chrs<-c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', 
      '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')
},
GSE117593={
    n1=25
    n2=25
    type = "WGBS"
     chrs<-c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', 
      '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')
},
GSE148060={
    n1=21
    n2=32
    type = "RRBS"
    chrs<-c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', 
      '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')
},
GSE103886={
    n1=11
    n2=12
    type = "RRBS"
    chrs<-c('1','2','3','4','5','6','7','8','9','10',
      '11','12','13','14','15','16','17','18','19','X','Y')
},
{
    stop("Unknown dataset ID", call.=FALSE)
}
)

trt<-rep(c('group1','group2'), times = c(n1,n2))

# Prepare input files for HMM-DM
meth<-read.table(file = "meth_table", sep="\t")
cov<-read.table(file = "cov_table", sep="\t")
split_pos<-str_split_fixed(rownames(cov), ":", 2)
meth<-data.frame(chr = split_pos[,1], pos = split_pos[,2], meth)
cov<-data.frame(chr = split_pos[,1], pos = split_pos[,2], cov)
cpg_num<-0

for(chr in chrs){
    cov_chr<-cov[cov$chr == paste0("chr",chr),]
    meth_chr<-meth[meth$chr == paste0("chr",chr),]
    rownames(cov_chr)<-cov_chr$pos
    cov_chr$pos<-NULL
    cov_chr$chr<-NULL
    rownames(meth_chr)<-meth_chr$pos
    meth_chr$pos<-NULL
    meth_chr$chr<-NULL
    meth_chr[is.na(meth_chr)]<-0
    write.table(cov_chr, file=paste0("HMM_DM_cov_chr",chr), sep="\t", col.names=FALSE, row.names = TRUE, quote = FALSE)
    write.table(meth_chr, file=paste0("HMM_DM_meth_chr",chr), sep="\t", col.names=FALSE, row.names = TRUE, quote = FALSE)
    cpg_num<-cpg_num+nrow(cov_chr)
}

# Save number of CpG sites for test
method_name<-"hmm"
cpg_num_table<-read.table(file = "cpg_num_table", sep="\t")
new_row<-data.frame(dataset_name, method_name, cpg_num)
names(new_row)<-c("dataset","method","cpg_num")
cpg_num_table<-rbind(cpg_num_table, new_row)
write.table(cpg_num_table, file = "cpg_num_table", sep="\t", col.names=TRUE, quote = FALSE)

# Run HMM-DM
run_hmm <- function(chr){
    total_reads<-read.table(paste0("HMM_DM_cov_chr",chr))
    meth_reads<-read.table(paste0("HMM_DM_meth_chr",chr))
    output_dir<-paste0("HMM_DM_chr",chr)
    system(paste0('mkdir ', output_dir))
    HMM.DM(total_reads, meth_reads, n1=n1, n2=n2, iterations=60, chromosome=chr, opt$hmm_dir, output_dir, meanDiff.cut = 0, min.percent = 0)
}

results <- mclapply(chrs, run_hmm, mc.cores = opt$cores)

cpgs<-data.frame(chr=character(), pos=numeric(), Hypo.pos=double(), EM.pos=double(), Hyper.pos=double(),
  max.p=double(), mCstatus=numeric(), meanDiff=double(), DM.status=numeric(), index=numeric(),
  meanCov.test=double(), meanCov.control=double(), stringsAsFactors=FALSE)

for(chr in chrs)
{
    filename<-file.path(paste0("HMM_DM_chr",chr),"DM.CG.txt")
    if (file.exists(filename))
    {
        cpgs_chr<-read.table(filename, header=TRUE)
        cpgs<-rbind(cpgs,cpgs_chr)
    }
}

# Get DMCs
dmls <-cpgs[cpgs$max.p >= 0.95,]
dmls_sorted<-dmls[order(-dmls$max.p),]
write.table(dmls_sorted, paste(dataset_name, method_name, "full", sep = "_"), sep = '\t')

# Filter by mean methylation difference (> 15%)
dmls_diff_cut<-dmls_sorted[abs(dmls_sorted$meanDiff) >= 0.15,]
write.table(dmls_diff_cut, paste(dataset_name, method_name, "cut", sep = "_"), sep = '\t')

# Get Hobotnica score
signature<-paste0(dmls_diff_cut$chr,":",dmls_diff_cut$pos)
sig_length<-length(signature)
result <-get_H_score(dataset_name, method_name, signature, trt, opt$cores)

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