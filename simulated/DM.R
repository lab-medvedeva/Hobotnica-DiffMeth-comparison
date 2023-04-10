library(methylSig)
library(bsseq)
library(foreach)
library(doParallel)
library(Hobotnica)
library(bsseq)
library(stringr)
library(parallel)
library(methylKit)
library(optparse)

option_list = list(
    make_option(c("-m", "--mode"), type="character", default=NULL,
        help="Mode: alt or perm", metavar="character"),
    make_option(c("-d", "--methyl_diff"), type="integer", default=10,
        help="Methylation difference: 10, 15, 20, 30"),
    make_option(c("-c", "--cores"), type="integer", default=16,
        help="Number of cores"),
    make_option(c("-h", "--hmm_code_folder"), type="character", default=NULL,
        help="HMM-DM source directory", metavar="character"))
    make_option(c("-g", "--group"), type="character", default=1,
        help="Group: 1 or 2"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
cores <- opt$cores
mode <- opt$mode
methyl_diff <- opt$methyl_diff
hmm_code_folder <- opt$hmm_code_folder
group <- opt$group
type = "RRBS"
dataset_name = paste(dname, methyl_diff, sep = '_')

############################# Init #############################

if(group == 1)
{
	dname <- 'default'
	samples <- paste0('sample', 0:31) 
	n1=16
	n2=16
	assembly = "hg38"
	files <- c(paste0(0:15, '_base.cov'), paste0(16:31, '_alt_', methyl_diff, '.sorted.cov')) 
    chrom = 22
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
    chrom = 19
}

if(mode =='alt')
{
	trt<-rep(c('group1','group2'), times = c(n1,n2))
	trt_bin<-trt
	trt_bin[trt_bin == 'group1'] <- 0
	trt_bin[trt_bin == 'group2'] <- 1
	trt_bin <- as.numeric(trt_bin)
	samples_group1 <- samples[1:n1]
	samples_group2 <- samples[(n1+1):(n1+n2)]

}
if(mode =='perm')
{
	trt <- sample(trt, n1+n2, replace = FALSE)
	trt_df <- data.frame(trt)
	write.table(trt, file = "trt_new.csv" , sep="\t", col.names=TRUE, quote = FALSE)	
	trt <- read.csv("trt_new.csv" , sep="\t")
	trt<- trt$x
	trt_bin<-trt
	trt_bin[trt_bin == 'group1'] <- 0
	trt_bin[trt_bin == 'group2'] <- 1
	trt_bin <- as.numeric(trt_bin)
	base <- rep(1, n1+n2)
	samples_group1<-samples[trt_bin == 0]
	samples_group2<-samples[trt_bin == 1]
}


####################### Hobotnica functions ################

get_random_H<-function(sig_length, trt, ratio_table,length_per_process)
{
    random_sig<-ratio_table[sample(nrow(ratio_table), sig_length), ]
    distMatrix<-dist(t(random_sig))
    if(sum(is.na(distMatrix))>0)
    {
        return(0)
    }
    return(Hobotnica(distMatrix, trt))
}

get_H_score<-function(dataset_name, method, signature, trt, cores, methyl_diff, sig_type) 
{
    sig_length<-length(signature)

    if (sig_length < 2)
    {
        return(list(H=0, pvalue=NA))
    }
    
    ratio_table<-read.table(file = paste0("ratio_table_", methyl_diff), sep="\t")
    heat_matrix<-ratio_table[signature,]
    distMatrix<-dist(t(heat_matrix))
    if(sum(is.na(distMatrix))>0)
    {
        return(list(H=0, pvalue=NA))
    }

    H<-Hobotnica(distMatrix, trt)

    # P-value calculation
    registerDoParallel(cores)
    H_random<-foreach(i = 1:5000, .combine = "c") %dopar% 
    {
        H_chunk<-get_random_H(sig_length, trt, ratio_table)
        H_chunk
    }
    H_random<-do.call(c, as.data.frame(H_random))
    H_random <- H_random[H_random != 0]

    # Save random values
    random_table_name<-paste(dataset_name, method, "H_random", sig_length, methyl_diff, sig_type, sep = "_")
    write.table(H_random, random_table_name, sep="\t", col.names=FALSE, row.names = FALSE, quote = FALSE)

    if(length(H_random) < 5000)
    {
        return(list(H=H, pvalue=NA))
    }

    # pvalue calculation with pseudo count
    pvalue<-(sum(H_random >= H)+1)/(length(H_random)+1)
    return(list(H=H, pvalue=pvalue))
}

get_H_score <- function(dataset_name, method, signature, trt, cores, methyl_diff, sig_type) 
{
    sig_length<-length(signature)

    if (sig_length < 2)
    {
        return(list(H=0, pvalue=NA))
    }
    
    ratio_table<-read.table(file = paste0("ratio_table_", methyl_diff), sep="\t")
    heat_matrix<-ratio_table[signature,]

    distMatrix<-dist(t(heat_matrix))
    if(sum(is.na(distMatrix))>0)
    {
        return(list(H=0, pvalue=NA))
    }

    H<-Hobotnica(distMatrix, trt)

    # P-value calculation
    registerDoParallel(cores)
    H_random<-foreach(i = 1:5000, .combine = "c") %dopar% 
    {
        H_chunk<-get_random_H(sig_length, trt, ratio_table)
        H_chunk
    }
    H_random<-do.call(c, as.data.frame(H_random))
    H_random <- H_random[H_random != 0]

    # Save random values
    random_table_name<-paste(dataset_name, method, "H_random", sig_length, methyl_diff, sig_type, sep = "_")
    write.table(H_random, random_table_name, sep="\t", col.names=FALSE, row.names = FALSE, quote = FALSE)

    if(length(H_random) < 5000)
    {
        return(list(H=H, pvalue=NA))
    }

    # pvalue calculation with pseudo count
    pvalue<-(sum(H_random >= H)+1)/(length(H_random)+1)
    return(list(H=H, pvalue=pvalue))
}

dataset_name = paste(dname, methyl_diff, sep = '_')
BiocParallel::register(BiocParallel::MulticoreParam(workers = cores, progressbar = TRUE))

# Coverage filtration
bsseq_stranded = bsseq::read.bismark(
    files = files,
    colData = data.frame(row.names = samples),
    rmZeroCov = TRUE,
    strandCollapse = TRUE
)

pData(bsseq_stranded)$Condition <- trt

# Low coverage loci (< 5) are marked
bsseq_stranded_filtered = filter_loci_by_coverage(bsseq_stranded, min_count = 5)

# Require at least two samples from group1 and two samples from group2
bsseq_stranded_filtered = filter_loci_by_group_coverage(
    bs = bsseq_stranded_filtered,
    group_column = 'Condition',
    c('group1' = 2, 'group2' = 2))

save(bsseq_stranded_filtered, file = paste0(dataset_name, "_", methyl_diff, ".rda"))

### Save methylation ratio and coverage tables ###
meth <- bsseq::getCoverage(bsseq_stranded_filtered, type = 'M')
names(meth) <- samples
cov <- bsseq::getCoverage(bsseq_stranded_filtered, type = 'Cov')
names(cov) <- samples
pos <- GenomicRanges::granges(bsseq_stranded_filtered)
pos_df<-data.frame(pos)
rownames(cov) <-paste0(pos_df$seqnames,":",pos_df$start)
rownames(meth) <-paste0(pos_df$seqnames,":",pos_df$start)
ratio_table = meth/cov

write.table(ratio_table, file = paste0("ratio_table_", methyl_diff), sep="\t", col.names=TRUE, row.names = TRUE, quote = FALSE)
write.table(meth, file = paste0("meth_table_", methyl_diff), sep="\t", col.names=TRUE, row.names = TRUE, quote = FALSE)
write.table(cov, file = paste0("cov_table_", methyl_diff), sep="\t", col.names=TRUE, row.names = TRUE, quote = FALSE)
write.table(pos_df, file = paste0("pos_table_", methyl_diff), sep="\t", col.names=TRUE, row.names = TRUE, quote = FALSE)

# Save number of CpG sites for test for several methods
cpg_num <- dim(bsseq_stranded_filtered)[1]
methods <- c("methylsig", "dss_no_smoothing", "dss_smoothing","methylkit", "methylkit_overdispersion", "radmeth")
cpg_num_table <- data.frame(dataset = dataset_name, method = methods, cpg_num = cpg_num)
write.table(cpg_num_table, file = "cpg_num_table", sep="\t", col.names=TRUE, quote = FALSE)

# Create empty Hobotnica scores tables
H_full <- data.frame(dataset=character(),
	method=character(),
	H=double(),
	type=character(),
	pvalue=double(),
	length=numeric(),
	stringsAsFactors=FALSE)
write.table(H_full, file = paste0("H_full_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

H_100 <- data.frame(dataset=character(),
	method=character(),
	H=double(),
	type=character(),
	pvalue=double(),
	length=numeric(),
	stringsAsFactors=FALSE)
write.table(H_100, file = paste0("H_100_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

H_10 <- data.frame(dataset=character(),
	method=character(),
	H=double(),
	type=character(),
	pvalue=double(),
	length=numeric(),
	stringsAsFactors=FALSE)
write.table(H_10, file = paste0("H_10_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

################################## methylsig #########
 
method_name <- 'methylsig'
# Load filtered data
load(paste0(dataset_name,"_",methyl_diff,".rda"))

# # Get DMC signature
diff_gr<-diff_methylsig(
 bs = bsseq_stranded_filtered,
 group_column = 'Condition',
 comparison_groups = c('control' = 'group1', 'case' = 'group2'),
 disp_groups = c('case' = TRUE, 'control' = TRUE),
 local_window_size = 0,
 t_approx = TRUE,
 n_cores = cores)


dmls_fdr<-diff_gr[!is.na(mcols(diff_gr)$fdr),]
dmls_fdr<-dmls_fdr[mcols(dmls_fdr)$fdr <= 0.05,]
dmls_fdr_sorted<-dmls_fdr[order(dmls_fdr$fdr),]
write.table(dmls_fdr_sorted, paste(dataset_name, method_name, "full", methyl_diff, sep = "_"), sep = '\t')

dmls_diff_sorted<-dmls_fdr_sorted[order(-abs(dmls_fdr_sorted$meth_diff)),]
write.table(dmls_diff_sorted, paste(dataset_name, method_name, "meth_diff_full", methyl_diff, sep = "_"), sep = '\t')

# Filter by mean methylation difference (> 15%)
#dmls_fdr_sorted <- read.csv(paste0(dataset_name,"_methylsig_full"), sep = '\t')
dmls_fdr_diff_cut<-dmls_fdr_sorted[abs(dmls_fdr_sorted$meth_diff) >= 15,]
write.table(dmls_fdr_diff_cut, paste(dataset_name, method_name, "cut", methyl_diff, sep = "_"), sep = '\t')

# Get Hobotnica score

dmls_fdr_diff_cut <- read.csv(paste(dataset_name, method_name, "full", methyl_diff, sep = "_"), sep = '\t')

dmls<-data.frame(dmls_fdr_diff_cut)
if(dim(dmls)[1] == 0)
{
    signature <- c()
} else {
    signature<-paste0(dmls$seqnames,":",dmls$start)
}

sig_length<-length(signature)
result <-get_H_score(dataset_name, method_name, signature, trt, cores, methyl_diff, 'full') 

# Write H score to result table
H_full<-read.table(paste0("H_full_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_full<-rbind(H_full, new_row)
write.table(H_full, file = paste0("H_full_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

# Get Hobotnica score for top 100 signature
signature_top_100<-head(signature, 100)
sig_length<-length(signature_top_100)
result <-get_H_score(dataset_name, method_name, signature_top_100, trt, cores, methyl_diff, 'full') 

# Write H score to result table
H_100<-read.table(paste0("H_100_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_100<-rbind(H_100, new_row)
write.table(H_100, file = paste0("H_100_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

# Get Hobotnica score for top 10 signature
signature_top_10<-head(signature, 10)
sig_length<-length(signature_top_10)
result <-get_H_score(dataset_name, method_name, signature_top_10, trt, cores, methyl_diff, 'full') 

# Write H score to result table
H_10<-read.table(paste0("H_10_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_10<-rbind(H_10, new_row)
write.table(H_10, file = paste0("H_10_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

################################### DSS ###############################  
 
method_name = "dss_no_smoothing"
smoothing=FALSE

#Get DMC signature
dmlTest = DMLtest(bsseq_stranded_filtered, group1= samples_group1, group2=samples_group2, smoothing=smoothing)
dmls <- callDML(dmlTest, p.threshold=1)
dmls_fdr  <- dmls[which(dmls$fdr <= 0.05),]
dmls_fdr_sorted<-dmls_fdr[order(dmls_fdr$fdr),]
write.table(dmls_fdr_sorted, paste(dataset_name, method_name, "full", methyl_diff, sep = "_"), sep = '\t')

dmls_diff_sorted<-dmls_fdr_sorted[order(-abs(dmls_fdr_sorted$diff)),]
write.table(dmls_diff_sorted, paste(dataset_name, method_name, "meth_diff_full", methyl_diff, sep = "_"), sep = '\t')

# Filter by mean methylation difference (> 15%)
dmls_fdr_diff_cut<-dmls_fdr_sorted[abs(dmls_fdr_sorted$diff) >= 0.15,]
write.table(dmls_fdr_diff_cut, paste(dataset_name, method_name, "cut", methyl_diff, sep = "_"), sep = '\t')

# Get Hobotnica score

dmls_fdr_sorted <- read.csv(paste(dataset_name, method_name, "full", methyl_diff, sep = "_"), sep = '\t')
if(dim(dmls_fdr_sorted)[1] == 0)
{
    signature <- c()
} else {
    signature<-paste0(dmls_fdr_sorted$chr,":",dmls_fdr_sorted$pos)
}

sig_length <- length(signature)
result <-get_H_score(dataset_name, method_name, signature, trt, cores, methyl_diff, "full") 

# Write H score to result table
H_full <- read.table(paste0("H_full_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_full <- rbind(H_full, new_row)
write.table(H_full, file = paste0("H_full_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

# Get Hobotnica score for top 100 signature
signature_top_100<-head(signature, 100)
sig_length <- length(signature_top_100)
result <-get_H_score(dataset_name, method_name, signature_top_100, trt, cores, methyl_diff, "full") 

# Write H score to result table
H_100 <- read.table(paste0("H_100_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_100 <- rbind(H_100, new_row)
write.table(H_100, file = paste0("H_100_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

# Get Hobotnica score for top 10 signature
signature_top_10<-head(signature, 10)
sig_length <- length(signature_top_10)
result <-get_H_score(dataset_name, method_name, signature_top_10, trt, cores, methyl_diff, "full") 

# Write H score to result table
H_10 <- read.table(paste0("H_10_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_10 <- rbind(H_10, new_row)
write.table(H_10, file = paste0("H_10_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

method_name = "dss_smoothing"

smoothing=TRUE

# Get DMC signature
dmlTest = DMLtest(bsseq_stranded_filtered, group1= samples_group1, group2=samples_group2, smoothing=smoothing)
dmls <- callDML(dmlTest, p.threshold=1)
dmls_fdr  <- dmls[which(dmls$fdr <= 0.05),]
dmls_fdr_sorted<-dmls_fdr[order(dmls_fdr$fdr),]
write.table(dmls_fdr_sorted, paste(dataset_name, method_name, "full", methyl_diff, sep = "_"), sep = '\t')

dmls_diff_sorted<-dmls_fdr_sorted[order(-abs(dmls_fdr_sorted$diff)),]
write.table(dmls_diff_sorted, paste(dataset_name, method_name, "meth_diff_full", methyl_diff, sep = "_"), sep = '\t')

# Filter by mean methylation difference (> 15%)
dmls_fdr_diff_cut<-dmls_fdr_sorted[abs(dmls_fdr_sorted$diff) >= 0.15,]
write.table(dmls_fdr_diff_cut, paste(dataset_name, method_name, "cut", methyl_diff, sep = "_"), sep = '\t')

# Get Hobotnica score
dmls_fdr_sorted <- read.csv(paste(dataset_name, method_name, "full", methyl_diff, sep = "_"), sep = '\t')
if(dim(dmls_fdr_sorted)[1] == 0)
{
    signature <- c()
} else {
    signature<-paste0(dmls_fdr_sorted$chr,":",dmls_fdr_sorted$pos)
}

sig_length <- length(signature)
result <-get_H_score(dataset_name, method_name, signature, trt, cores, methyl_diff, "full") 

# Write H score to result table
H_full <- read.table(paste0("H_full_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_full <- rbind(H_full, new_row)
write.table(H_full, file = paste0("H_full_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

# Get Hobotnica score for top 100 signature
signature_top_100<-head(signature, 100)
sig_length <- length(signature_top_100)
result <-get_H_score(dataset_name, method_name, signature_top_100, trt, cores, methyl_diff, "full") 

# Write H score to result table
H_100 <- read.table(paste0("H_100_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_100 <- rbind(H_100, new_row)
write.table(H_100, file = paste0("H_100_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

# Get Hobotnica score for top 10 signature
signature_top_10<-head(signature, 10)
sig_length <- length(signature_top_10)
result <-get_H_score(dataset_name, method_name, signature_top_10, trt, cores, methyl_diff, "full") 

# Write H score to result table
H_10 <- read.table(paste0("H_10_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_10 <- rbind(H_10, new_row)
write.table(H_10, file = paste0("H_10_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

# Load data
pos_df<-read.table(file = paste0("pos_table_", methyl_diff), sep="\t")
ratio_table<-read.table(file = paste0("ratio_table_", methyl_diff), sep="\t")
cov<-read.table(file = paste0("cov_table_", methyl_diff), sep="\t")

# Prepare input files for methylkit
methylkit_filenames<-c()

for(sample in samples)
{
    sample_table<-data.frame(V1 = paste0(pos_df$seqnames,".",pos_df$start),V2 = pos_df$seqnames, V3=pos_df$start)
    sample_table$V4<-'*'
    sample_table$V5<-cov[,sample]
    sample_table$V6<-round(as.numeric(ratio_table[,sample])*100,2)
    sample_table$V7<-round(100-as.numeric(ratio_table[,sample])*100,2)
    sample_table<-sample_table[sample_table$V5 > 0,]
    methylkit_filename<-paste("methylkit",sample, methyl_diff, sep="_")
    print(paste("Writing to",methylkit_filename))
    write.table(sample_table, file=methylkit_filename, sep="\t", col.names=FALSE, row.names = FALSE, quote = FALSE)
    methylkit_filenames <-c(methylkit_filenames,methylkit_filename)
}

method_name = "methylkit"
overdispersion <- FALSE
methylkit_filenames<-c()
for(sample in samples)
{
    methylkit_filename<-paste("methylkit",sample, methyl_diff, sep="_")
    methylkit_filenames <-c(methylkit_filenames,methylkit_filename)
}

# Get DMC signature
meth_obj=methRead(as.list(methylkit_filenames),
    sample.id=as.list(samples),
    assembly=assembly,
    treatment=trt_bin,
    context="CpG" )

meth<-unite(meth_obj, mc.cores=cores)
if (overdispersion) {
    meth_diff<-calculateDiffMeth(meth, overdispersion="MN", test="Chisq", mc.cores = cores)
} else {
    meth_diff<-calculateDiffMeth(meth,  mc.cores = cores)
}

dmls_qvalue <-getMethylDiff(meth_diff,difference=0,qvalue=0.05)
dmls_sorted<-dmls_qvalue[order(dmls_qvalue$qvalue),]
write.table(dmls_sorted, paste(dataset_name, method_name, "full", methyl_diff,  sep = "_"), sep = '\t')

dmls_diff_sorted<-dmls_sorted[order(-abs(dmls_sorted$meth.diff)),]
write.table(dmls_diff_sorted, paste(dataset_name, method_name, "meth_diff_full", methyl_diff, sep = "_"), sep = '\t')


# Filter by mean methylation difference (> 15%)
dmls_diff_cut<-dmls_sorted[abs(dmls_sorted$meth.diff) >= 15,]
write.table(dmls_diff_cut, paste(dataset_name, method_name, "cut", methyl_diff, sep = "_"), sep = '\t')

# Get Hobotnica score
dmls_diff_cut <- read.csv(paste(dataset_name, method_name, "full", methyl_diff,  sep = "_"), sep = '\t')
if(dim(dmls_diff_cut)[1] == 0)
{
    signature <- c()
} else {
    signature<-paste0(dmls_diff_cut$chr,":",dmls_diff_cut$start)
}

sig_length<-length(signature)
result <-get_H_score(dataset_name, method_name, signature, trt, cores, methyl_diff, "full")

# Write H score to result table
H_full<-read.table(paste0("H_full_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_full<-rbind(H_full, new_row)
write.table(H_full, file = paste0("H_full_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

# Get Hobotnica score for top 100 signature
signature_top_100<-head(signature, 100)
sig_length<-length(signature_top_100)
result <-get_H_score(dataset_name, method_name, signature_top_100, trt, cores, methyl_diff, "full") 

# Write H score to result table
H_100<-read.table(paste0("H_100_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_100<-rbind(H_100, new_row)
write.table(H_100, file = paste0("H_100_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

# Get Hobotnica score for top 10 signature
signature_top_10<-head(signature, 10)
sig_length<-length(signature_top_10)
result <-get_H_score(dataset_name, method_name, signature_top_10, trt, cores, methyl_diff, "full") 

# Write H score to result table
H_10<-read.table(paste0("H_10_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_10<-rbind(H_10, new_row)
write.table(H_10, file = paste0("H_10_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

dataset_name = paste(dname, methyl_diff, sep = '_')
method_name = "methylkit_overdispersion"
overdispersion <- TRUE
load(paste0(dataset_name,"_",methyl_diff,".rda"))

methylkit_filenames<-c()

for(sample in samples)
{
    methylkit_filename<-paste("methylkit",sample, methyl_diff, sep="_")
    methylkit_filenames <-c(methylkit_filenames,methylkit_filename)
}


# Get DMC signature
meth_obj=methRead(as.list(methylkit_filenames),
    sample.id=as.list(samples),
    assembly=assembly,
    treatment=trt_bin,
    context="CpG" )

meth<-unite(meth_obj, mc.cores=cores)
if (overdispersion) {
    meth_diff<-calculateDiffMeth(meth, overdispersion="MN", test="Chisq", mc.cores = cores)
} else {
    meth_diff<-calculateDiffMeth(meth,  mc.cores = cores)
}

dmls_qvalue <-getMethylDiff(meth_diff,difference=0,qvalue=0.05)
dmls_sorted<-dmls_qvalue[order(dmls_qvalue$qvalue),]
write.table(dmls_sorted, paste(dataset_name, method_name, "full", methyl_diff,  sep = "_"), sep = '\t')

dmls_diff_sorted<-dmls_sorted[order(-abs(dmls_sorted$meth.diff)),]
write.table(dmls_diff_sorted, paste(dataset_name, method_name, "meth_diff_full", methyl_diff, sep = "_"), sep = '\t')

# Filter by mean methylation difference (> 15%)
dmls_diff_cut<-dmls_sorted[abs(dmls_sorted$meth.diff) >= 15,]
write.table(dmls_diff_cut, paste(dataset_name, method_name, "cut", methyl_diff, sep = "_"), sep = '\t')

# Get Hobotnica score

dmls_diff_cut <- read.csv(paste(dataset_name, method_name, "full", methyl_diff,  sep = "_"), sep = '\t')
if(dim(dmls_diff_cut)[1] == 0)
{
    signature <- c()
} else {
    signature<-paste0(dmls_diff_cut$chr,":",dmls_diff_cut$start)
}

sig_length<-length(signature)
result <-get_H_score(dataset_name, method_name, signature, trt, cores, methyl_diff, "full")

# Write H score to result table
H_full<-read.table(paste0("H_full_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_full<-rbind(H_full, new_row)
write.table(H_full, file = paste0("H_full_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

# Get Hobotnica score for top 100 signature
signature_top_100<-head(signature, 100)
sig_length<-length(signature_top_100)
result <-get_H_score(dataset_name, method_name, signature_top_100, trt, cores, methyl_diff, "full") 

# Write H score to result table
H_100<-read.table(paste0("H_100_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_100<-rbind(H_100, new_row)
write.table(H_100, file = paste0("H_100_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

# Get Hobotnica score for top 10 signature
signature_top_10<-head(signature, 10)
sig_length<-length(signature_top_10)
result <-get_H_score(dataset_name, method_name, signature_top_10, trt, cores, methyl_diff, "full") 

# Write H score to result table
H_10<-read.table(paste0("H_10_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_10<-rbind(H_10, new_row)
write.table(H_10, file = paste0("H_10_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)


####################### hmm #######################

source(file.path(hmm_code_folder,"HMM.DM.R"))
method_name = "hmm"

# Prepare input files for HMM-DM
meth<-read.table(file = paste0("meth_table_",methyl_diff), sep="\t")
cov<-read.table(file =  paste0("cov_table_",methyl_diff), sep="\t")
split_pos<-str_split_fixed(rownames(cov), ":", 2)
meth<-data.frame(chr = split_pos[,1], pos = split_pos[,2], meth)
cov<-data.frame(chr = split_pos[,1], pos = split_pos[,2], cov)
cpg_num<-0

cov_chr<-cov
meth_chr<-meth
rownames(cov_chr)<-cov_chr$pos
cov_chr$pos<-NULL
cov_chr$chr<-NULL
rownames(meth_chr)<-meth_chr$pos
meth_chr$pos<-NULL
meth_chr$chr<-NULL
meth_chr[is.na(meth_chr)]<-0
write.table(cov_chr, file=paste0("HMM_DM_cov_",methyl_diff), sep="\t", col.names=FALSE, row.names = TRUE, quote = FALSE)
write.table(meth_chr, file=paste0("HMM_DM_meth_",methyl_diff), sep="\t", col.names=FALSE, row.names = TRUE, quote = FALSE)
cpg_num<-cpg_num+nrow(cov_chr)

cpg_num_table<-read.table(file = "cpg_num_table", sep="\t")
new_row<-data.frame(dataset_name, method_name, cpg_num)
names(new_row)<-c("dataset","method","cpg_num")
cpg_num_table<-rbind(cpg_num_table, new_row)
write.table(cpg_num_table, file = "cpg_num_table", sep="\t", col.names=TRUE, quote = FALSE)

# Run HMM-DM
total_reads<-read.table(paste0("HMM_DM_cov_",methyl_diff))
total_reads_fixed <- total_reads[,c(1, which(trt == "group1") + 1,which(trt == "group2") + 1)]
meth_reads<-read.table(paste0("HMM_DM_meth_",methyl_diff))
meth_reads_fixed <- meth_reads[,c(1, which(trt == "group1") + 1,which(trt == "group2") + 1)]
output_dir<-paste0("HMM_DM_",methyl_diff)
system(paste0('mkdir ', output_dir))
HMM.DM(total_reads, meth_reads, n1=n1, n2=n2, iterations=60, chromosome=chrom, hmm_code_folder, output_dir, meanDiff.cut = 0, min.percent = 0)

dataset_name = paste(dname, methyl_diff, sep = '_')
# Save number of CpG sites for test

filename<-file.path(paste0("HMM_DM_",methyl_diff), "DM.CG.txt")
cpgs<-read.table(filename, header=TRUE)

# Get DMCs
dmls <-cpgs[cpgs$max.p >= 0.95,]
dmls <- rbind(dmls[dmls$Hypo.pos == 1,], dmls[dmls$Hyper.pos == 1,])
dmls_sorted<-dmls[order(-dmls$max.p),]

write.table(dmls_sorted, paste(dataset_name, method_name, "full", methyl_diff, sep = "_"), sep = '\t')


dmls_diff_sorted<-dmls_sorted[order(-abs(dmls_sorted$meanDiff)) ,]
write.table(dmls_diff_sorted, paste(dataset_name, method_name, "meth_diff_full", methyl_diff, sep = "_"), sep = '\t')

# Filter by mean methylation difference (> 15%)
dmls_diff_cut<-dmls_sorted[abs(dmls_sorted$meanDiff) >= 0.15,]
write.table(dmls_diff_cut, paste(dataset_name, method_name, "cut", methyl_diff, sep = "_"), sep = '\t')

# Get Hobotnica score

dmls_diff_cut <- read.csv(paste(dataset_name, method_name, "full", methyl_diff, sep = "_"), sep = '\t')
if(dim(dmls_diff_cut)[1] == 0)
{
    signature <- c()
} else {
    signature<-paste0(dmls_diff_cut$chr,":",dmls_diff_cut$pos)
}

sig_length<-length(signature)
result <-get_H_score(dataset_name, method_name, signature, trt, cores, methyl_diff, "full")

# Write H score to result table
H_full<-read.table(paste0("H_full_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_full<-rbind(H_full, new_row)
write.table(H_full, file = paste0("H_full_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

# Get Hobotnica score for top 100 signature
signature_top_100<-head(signature, 100)
sig_length<-length(signature_top_100)
result <-get_H_score(dataset_name, method_name, signature_top_100, trt, cores, methyl_diff, "full") 

# Write H score to result table
H_100<-read.table(paste0("H_100_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_100<-rbind(H_100, new_row)
write.table(H_100, file = paste0("H_100_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

# Get Hobotnica score for top 10 signature
signature_top_10<-head(signature, 10)
sig_length<-length(signature_top_10)
result <-get_H_score(dataset_name, method_name, signature_top_10, trt, cores, methyl_diff, "full") 

# Write H score to result table
H_10<-read.table(paste0("H_10_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_10<-rbind(H_10, new_row)
write.table(H_10, file = paste0("H_10_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

############################### radmeth ###############################

method_name = "radmeth"

# Load data
meth<-read.table(file = paste0("meth_table_", methyl_diff), sep="\t")
cov<-read.table(file = paste0("cov_table_", methyl_diff), sep="\t")

# Prepare input files for radmeth
# Data table
radmeth_data<-data.frame(matrix(ncol = 2*length(samples), nrow = nrow(cov)))
i<-1
for (col in samples){
  radmeth_data[,i]<-cov[,col]
  radmeth_data[,i+1]<-meth[,col]
  i<-i+2
}
rownames(radmeth_data)<-paste0(rownames(cov),":+:CpG")
write.table(radmeth_data, file= paste0("radmeth_table_", methyl_diff), sep="\t", col.names=FALSE, row.names = TRUE, quote = FALSE)
radmeth_header<-paste(samples, collapse = '\t')
system(paste0("sed -i '1i", radmeth_header, "\' ", "radmeth_table_", methyl_diff))

# Design matrix
design_matrix<-data.frame(base,trt_bin)
rownames(design_matrix)<-samples
write.table(design_matrix, file="radmeth_design_matrix", sep="\t", col.names=TRUE, row.names = TRUE, quote = FALSE)
system(paste0('radmeth regression -factor trt_bin radmeth_design_matrix radmeth_table_', methyl_diff, ' > radmeth_result_', methyl_diff, '.bed'))
system(paste0('radmeth adjust -bins 1:200:1 radmeth_result_', methyl_diff, '.bed > radmeth_result_adjusted_', methyl_diff, '.bed'))
cpgs<-read.table(paste0("radmeth_result_adjusted_", methyl_diff, ".bed"), header=FALSE)
cpgs<-cpgs[cpgs$V7 != "-1-1",]
write.table(cpgs, file=paste0("radmeth_result_adjusted_filtered_", methyl_diff, ".bed"), sep="\t", col.names=FALSE, row.names = FALSE, quote = FALSE)

# Get DMC signature
cpgs<-read.table(paste0("radmeth_result_adjusted_filtered_", methyl_diff, ".bed"),header=FALSE)
dmls<-cpgs[which(cpgs$V7 <= 0.05),]
dmls_sorted<-dmls[order(dmls$V7),]
write.table(dmls_sorted, paste(dataset_name, method_name, "full", methyl_diff, sep = "_"), sep = '\t')


# Filter by mean methylation difference (> 15%)
ratio_table<-read.table(file = paste0("ratio_table_", methyl_diff), sep="\t")
dmls_sorted_meth<-ratio_table[paste0(dmls_sorted$V1,":",dmls_sorted$V2),]
meandif = abs(rowMeans(dmls_sorted_meth[,samples_group1], na.rm = TRUE) - rowMeans(dmls_sorted_meth[,samples_group2], na.rm = TRUE))

dmls_diff_sorted<-dmls_sorted[order(-abs(meandif)) ,]
write.table(dmls_diff_sorted, paste(dataset_name, method_name, "meth_diff_full", methyl_diff, sep = "_"), sep = '\t')

if (dim(dmls_sorted)[1] > 0)
{
	dmls_diff_cut<-dmls_sorted[abs(meandif) >= 0.15,]
	write.table(dmls_diff_cut, paste(dataset_name, method_name, "cut", methyl_diff, sep = "_"), sep = '\t')
}
else
{
	dmls_diff_cut<-dmls_sorted
	write.table(dmls_diff_cut, paste(dataset_name, method_name, "cut", methyl_diff, sep = "_"), sep = '\t')
}

# Get Hobotnica score
dmls_diff_cut <- read.csv(paste(dataset_name, method_name, "full", methyl_diff,  sep = "_"), sep = '\t')
if(dim(dmls_diff_cut)[1] == 0)
{
    signature <- c()
} else {
    signature<-paste0(dmls_diff_cut$V1,":",dmls_diff_cut$V2)
}

sig_length<-length(signature)
result <-get_H_score(dataset_name, method_name, signature, trt, cores, methyl_diff, "full")

# Write H score to result table
H_full<-read.table(paste0("H_full_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_full<-rbind(H_full, new_row)
write.table(H_full, file = paste0("H_full_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

# Get Hobotnica score for top 100 signature
signature_top_100<-head(signature, 100)
sig_length<-length(signature_top_100)
result <-get_H_score(dataset_name, method_name, signature_top_100, trt, cores, methyl_diff, "full") 

# Write H score to result table
H_100<-read.table(paste0("H_100_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_100<-rbind(H_100, new_row)
write.table(H_100, file = paste0("H_100_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)

# Get Hobotnica score for top 10 signature
signature_top_10<-head(signature, 10)
sig_length<-length(signature_top_10)
result <-get_H_score(dataset_name, method_name, signature_top_10, trt, cores, methyl_diff, "full") 

# Write H score to result table
H_10<-read.table(paste0("H_10_", methyl_diff), header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, "full", result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_10<-rbind(H_10, new_row)
write.table(H_10, file = paste0("H_10_", methyl_diff), sep="\t", col.names=TRUE, quote = FALSE)
 
