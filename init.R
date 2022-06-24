#!/usr/bin/env Rscript
library(optparse)
library(methylSig)
library(bsseq)

# Parse command line arguments
option_list = list(
	make_option(c("-d", "--dataset"), type="character", default=NULL,
		help="Dataset GEO id", metavar="character"),
	make_option(c("-w", "--working_dir"), type="character", default=NULL,
		help="Working directory", metavar="character")); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dataset)){
	print_help(opt_parser)
	stop("Dataset id should be supplied", call.=FALSE)
}

if (is.null(opt$working_dir)){
	print_help(opt_parser)
	stop("Working directory should be supplied", call.=FALSE)
}

# Set working directory
setwd(opt$working_dir)
# Dataset initialization
dataset_name <- opt$dataset
switch(dataset_name, 
GSE149608={
	n1=10
	n2=10
	samples <- c('SRR11647648', 'SRR11647649', 'SRR11647650', 'SRR11647651', 'SRR11647652', 'SRR11647653',
		'SRR11647654', 'SRR11647655', 'SRR11647656', 'SRR11647657', 'SRR11647658', 'SRR11647659', 'SRR11647660',
		'SRR11647661', 'SRR11647662', 'SRR11647663', 'SRR11647664', 'SRR11647665', 'SRR11647666', 'SRR11647667')
},
GSE138598={
	n1=8
	n2=9
	samples <- c('T2D_1', 'T2D_2', 'T2D_3', 'T2D_4', 'T2D_5', 'T2D_6', 'T2D_7', 'T2D_8',
		'Control_1', 'Control_2', 'Control_3', 'Control_4', 'Control_5', 'Control_6', 'Control_7', 
		'Control_8', 'Control_9')
},
GSE119980={
	n1=6
	n2=6
	samples <-c('SRR7830270', 'SRR7830271', 'SRR7830272', 'SRR7830273', 'SRR7830274', 'SRR7830275', 'SRR7830276',
		'SRR7830277', 'SRR7830278', 'SRR7830279', 'SRR7830280', 'SRR7830281')
},
GSE117593={
	n1=25
	n2=25
	samples <- c('SRR7587054', 'SRR7587055', 'SRR7587056', 'SRR7587057', 'SRR7587058', 'SRR7587059', 'SRR7587060',
		'SRR7587061', 'SRR7587062', 'SRR7587063', 'SRR7587064', 'SRR7587065', 'SRR7587066', 'SRR7587067', 'SRR7587068',
		'SRR7587069', 'SRR7587070', 'SRR7587071', 'SRR7587072', 'SRR7587073', 'SRR7587074', 'SRR7587075', 'SRR7587076',
		'SRR7587077', 'SRR7587078', 'SRR7587079', 'SRR7587080', 'SRR7587081', 'SRR7587082', 'SRR7587083', 'SRR7587084',
		'SRR7587085', 'SRR7587086', 'SRR7587087', 'SRR7587088', 'SRR7587089', 'SRR7587090', 'SRR7587091', 'SRR7587092',
		'SRR7587093', 'SRR7587094', 'SRR7587095', 'SRR7587096', 'SRR7587097', 'SRR7587098', 'SRR7587099', 'SRR7587100',
		'SRR7587101', 'SRR7587102', 'SRR7587103')
},
GSE148060={
	n1=21
	n2=32 
	samples <- c('SRR11477199', 'SRR11477200', 'SRR11477201', 'SRR11477202', 'SRR11477203', 'SRR11477204', 'SRR11477205',
		'SRR11477206', 'SRR11477207', 'SRR11477208', 'SRR11477209', 'SRR11477210', 'SRR11477211', 'SRR11477212', 'SRR11477213',
		'SRR11477214', 'SRR11477215', 'SRR11477216', 'SRR11477217', 'SRR11477218', 'SRR11477219', 'SRR11477220', 'SRR11477221',
		'SRR11477222', 'SRR11477223', 'SRR11477224', 'SRR11477225', 'SRR11477226', 'SRR11477227', 'SRR11477228', 'SRR11477229',
		'SRR11477230', 'SRR11477231', 'SRR11477232', 'SRR11477233', 'SRR11477234', 'SRR11477235', 'SRR11477236', 'SRR11477237',
		'SRR11477238', 'SRR11477239', 'SRR11477240', 'SRR11477241', 'SRR11477242', 'SRR11477243', 'SRR11477244', 'SRR11477245',
		'SRR11477246', 'SRR11477247', 'SRR11477248', 'SRR11477249', 'SRR11477250', 'SRR11477251')

},
GSE103886={
	n1=11
	n2=12
	samples <- c('SRR6040066', 'SRR6040067', 'SRR6040068', 'SRR6040069', 'SRR6040070', 'SRR6040071', 'SRR6040072', 'SRR6040073',
		'SRR6040074', 'SRR6040075', 'SRR6040076', 'SRR6040077', 'SRR6040078', 'SRR6040079', 'SRR6040080', 'SRR6040081', 'SRR6040082',
		'SRR6040083', 'SRR6040084', 'SRR6040085', 'SRR6040086', 'SRR6040087', 'SRR6040088')

},
{
	stop("Unknown dataset ID.n", call.=FALSE)
}
)

files <- paste0(samples, ".cov")
trt<-rep(c('group1','group2'), times = c(n1,n2))
samples_group1 <- samples[1:n1]
samples_group2 <- samples[(n1+1):(n1+n2)]

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

save(bsseq_stranded_filtered, file = paste0(dataset_name,".rda"))

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

write.table(ratio_table, file = "ratio_table", sep="\t", col.names=TRUE, row.names = TRUE, quote = FALSE)
write.table(meth, file = "meth_table", sep="\t", col.names=TRUE, row.names = TRUE, quote = FALSE)
write.table(cov, file = "cov_table", sep="\t", col.names=TRUE, row.names = TRUE, quote = FALSE)
write.table(pos_df, file = "pos_table", sep="\t", col.names=TRUE, row.names = TRUE, quote = FALSE)

# Save number of CpG sites for test for several methods
cpg_num <- dim(bsseq_stranded_filtered)[1]
methods <- c("methylsig", "dss_no_smoothing", "dss_smoothing","methylkit", "methylkit_overdispersion", "radmeth")
cpg_num_table <- data.frame(dataset = dataset_name, method = methods, cpg_num = cpg_num)
write.table(cpg_num_table, file = "cpg_num_table", sep="\t", col.names=TRUE, quote = FALSE)

# Save methylation differences
meandif = abs(rowMeans(ratio_table[,samples_group1], na.rm = TRUE) - rowMeans(ratio_table[,samples_group2], na.rm = TRUE))
write.table(meandif, file = paste0(dataset_name,"_meandif"), sep="\t", col.names=TRUE, quote = FALSE)

# Create empty Hobotnica scores tables
H_full <- data.frame(dataset=character(),
	method=character(),
	H=double(),
	type=character(),
	pvalue=double(),
	length=numeric(),
	stringsAsFactors=FALSE)
write.table(H_full, file = "H_full", sep="\t", col.names=TRUE, quote = FALSE)

H_100 <- data.frame(dataset=character(),
	method=character(),
	H=double(),
	type=character(),
	pvalue=double(),
	length=numeric(),
	stringsAsFactors=FALSE)
write.table(H_100, file = "H_100", sep="\t", col.names=TRUE, quote = FALSE)
