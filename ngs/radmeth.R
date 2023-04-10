#!/usr/bin/env Rscript
library(optparse)

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


trt<-rep(c('group1','group2'), times = c(n1,n2))
trt_bin<-rep(c(0, 1), times = c(n1,n2))
base<- rep(1, n1+n2)
samples_group1<-samples[1:n1]
samples_group2<-samples[(n1+1):(n1+n2)]

},
GSE138598={
    n1=8
    n2=9
    type = "WGBS"
    samples <- c('T2D_1', 'T2D_2', 'T2D_3', 'T2D_4', 'T2D_5', 'T2D_6', 'T2D_7', 'T2D_8',
        'Control_1', 'Control_2', 'Control_3', 'Control_4', 'Control_5', 'Control_6', 'Control_7', 
        'Control_8', 'Control_9')


trt<-rep(c('group1','group2'), times = c(n1,n2))
trt_bin<-rep(c(0, 1), times = c(n1,n2))
base<- rep(1, n1+n2)
samples_group1<-samples[1:n1]
samples_group2<-samples[(n1+1):(n1+n2)]

},
GSE119980={
    n1=6
    n2=6
    type = "WGBS"
    samples <-c('SRR7830270', 'SRR7830271', 'SRR7830272', 'SRR7830273', 'SRR7830274', 'SRR7830275', 'SRR7830276',
        'SRR7830277', 'SRR7830278', 'SRR7830279', 'SRR7830280', 'SRR7830281')


trt<-rep(c('group1','group2'), times = c(n1,n2))
trt_bin<-rep(c(0, 1), times = c(n1,n2))
base<- rep(1, n1+n2)
samples_group1<-samples[1:n1]
samples_group2<-samples[(n1+1):(n1+n2)]

},
GSE148060={
    n1=21
    n2=32
    type = "RRBS"
    samples<-c('SRR11477199', 'SRR11477200', 'SRR11477201', 'SRR11477202', 'SRR11477203', 'SRR11477204', 'SRR11477205',
        'SRR11477206', 'SRR11477207', 'SRR11477208', 'SRR11477209', 'SRR11477210', 'SRR11477211', 'SRR11477212', 'SRR11477213',
        'SRR11477214', 'SRR11477215', 'SRR11477216', 'SRR11477217', 'SRR11477218', 'SRR11477219', 'SRR11477220', 'SRR11477221',
        'SRR11477222', 'SRR11477223', 'SRR11477224', 'SRR11477225', 'SRR11477226', 'SRR11477227', 'SRR11477228', 'SRR11477229',
        'SRR11477230', 'SRR11477231', 'SRR11477232', 'SRR11477233', 'SRR11477234', 'SRR11477235', 'SRR11477236', 'SRR11477237',
        'SRR11477238', 'SRR11477239', 'SRR11477240', 'SRR11477241', 'SRR11477242', 'SRR11477243', 'SRR11477244', 'SRR11477245',
        'SRR11477246', 'SRR11477247', 'SRR11477248', 'SRR11477249', 'SRR11477250', 'SRR11477251')


trt<-rep(c('group1','group2'), times = c(n1,n2))
trt_bin<-rep(c(0, 1), times = c(n1,n2))
base<- rep(1, n1+n2)
samples_group1<-samples[1:n1]
samples_group2<-samples[(n1+1):(n1+n2)]

},
GSE103886={
    n1=11
    n2=12
    type = "RRBS"
    samples<-c('SRR6040066', 'SRR6040067', 'SRR6040068', 'SRR6040069', 'SRR6040070', 'SRR6040071', 'SRR6040072', 'SRR6040073',
        'SRR6040074', 'SRR6040075', 'SRR6040076', 'SRR6040077', 'SRR6040078', 'SRR6040079', 'SRR6040080', 'SRR6040081', 'SRR6040082',
        'SRR6040083', 'SRR6040084', 'SRR6040085', 'SRR6040086', 'SRR6040087', 'SRR6040088')

trt<-rep(c('group1','group2'), times = c(n1,n2))
trt_bin<-rep(c(0, 1), times = c(n1,n2))
base<- rep(1, n1+n2)
samples_group1<-samples[1:n1]
samples_group2<-samples[(n1+1):(n1+n2)]
},

 GSE150592=
 {
    n1=15
    n2=15
    type = "RRBS"
    samples <- c('SRR11790875','SRR11790876','SRR11790877','SRR11790878','SRR11790879','SRR11790880','SRR11790881','SRR11790882','SRR11790883','SRR11790884','SRR11790885','SRR11790886','SRR11790887','SRR11790888','SRR11790889','SRR11790890','SRR11790891','SRR11790892','SRR11790893','SRR11790894','SRR11790895','SRR11790896','SRR11790897','SRR11790898','SRR11790899','SRR11790900','SRR11790901','SRR11790902','SRR11790903','SRR11790904')
    trt<-c('group2','group1','group2','group1','group2','group1','group2','group1','group2','group1','group2','group1','group1','group2','group2','group1','group1','group2','group2','group1','group1','group2','group2','group1','group2','group1','group2','group1','group1','group2')
    trt_bin<-c(1,0,1,0,1,0,1,0,1,0,1,0,0,1,1,0,0,1,1,0,0,1,1,0,1,0,1,0,0,1)
    base<- rep(1, n1+n2)
    samples_group1 <- samples[which(trt == "group1")]
    samples_group2 <- samples[which(trt == "group2")]
},

{
    stop("Unknown dataset ID", call.=FALSE)
}
)

# Load data
meth<-read.table(file = "meth_table", sep="\t")
cov<-read.table(file = "cov_table", sep="\t")

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
write.table(radmeth_data, file="radmeth_table", sep="\t", col.names=FALSE, row.names = TRUE, quote = FALSE)
radmeth_header<-paste(samples, collapse = '\t')
system(paste0("sed -i '1i", radmeth_header, "\' radmeth_table"))

# Design matrix
design_matrix<-data.frame(base,trt_bin)
rownames(design_matrix)<-samples
write.table(design_matrix, file="radmeth_design_matrix", sep="\t", col.names=TRUE, row.names = TRUE, quote = FALSE)

# Run radmeth
system('radmeth regression -factor trt_bin radmeth_design_matrix radmeth_table > radmeth_result.bed')
system('radmeth adjust -bins 1:200:1 radmeth_result.bed > radmeth_result_adjusted.bed')
cpgs<-read.table("radmeth_result_adjusted.bed", header=FALSE)
cpgs<-cpgs[cpgs$V7 != "-1-1",]
write.table(cpgs, file="radmeth_result_adjusted_filtered", sep="\t", col.names=FALSE, row.names = FALSE, quote = FALSE)

# Get DMC signature
cpgs<-read.table("radmeth_result_adjusted_filtered",header=FALSE)
dmls<-cpgs[which(cpgs$V7 <= 0.05),]
dmls_sorted<-dmls[order(dmls$V7),]
method_name<-"radmeth"
write.table(dmls_sorted, paste(dataset_name, method_name, "full", sep = "_"), sep = '\t')

# Filter by mean methylation difference (> 15%)
ratio_table<-read.table(file = "ratio_table", sep="\t")
dmls_sorted_meth<-ratio_table[paste0(dmls_sorted$V1,":",dmls_sorted$V2),]
meandif = abs(rowMeans(dmls_sorted_meth[,samples_group1], na.rm = TRUE) - rowMeans(dmls_sorted_meth[,samples_group2], na.rm = TRUE))

dmls_diff_sorted<-dmls_sorted[order(-abs(meandif)) ,]
write.table(dmls_diff_sorted, paste(dataset_name, method_name, "meth_diff_full", sep = "_"), sep = '\t')


dmls_diff_cut<-dmls_sorted[abs(meandif) >= 0.15,]
write.table(dmls_diff_cut, paste(dataset_name, method_name, "cut", sep = "_"), sep = '\t')

# Get Hobotnica score
dmls_diff_cut <- read.csv(paste(dataset_name, method_name, "cut", sep = "_"), sep = '\t')
if(dim(dmls_diff_cut)[1] == 0)
{
    signature <- c()
} else {
    signature<-paste0(dmls_diff_cut$V1,":",dmls_diff_cut$V2)
}

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

# Get Hobotnica score for top 10 signature
signature_top_10<-head(signature, 10)
sig_length<-length(signature_top_10)
result <-get_H_score(dataset_name, method_name, signature_top_10, trt, opt$cores) 

# Write H score to result table
H_10<-read.table("H_10", header=TRUE, sep='\t')
new_row<-data.frame(dataset_name, method_name, result$H, type, result$pvalue, sig_length)
names(new_row)<-c("dataset","method","H","type","pvalue","length")
H_10<-rbind(H_10, new_row)
write.table(H_10, file = "H_10", sep="\t", col.names=TRUE, quote = FALSE)
