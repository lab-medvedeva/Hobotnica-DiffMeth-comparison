library(ChAMP)
library(GEOquery)
library(RnBeads)
library(minfi)
library(foreach)
library(doParallel)
library(Hobotnica)
library (ggplot2)

get_random_H<-function(random_sig, trt)
{
    distMatrix <- dist(t(random_sig))
    if(sum(is.na(distMatrix)) > 0)
    {
        return(0)
    }

    return(Hobotnica(distMatrix, trt))
}

get_H_score<-function(dataset_name, method, result_prefix, beta_values_matrix, signature, trt, cores) 
{
    sig_length <- length(signature)
    H <- get_H_score_int(signature, trt, beta_values_matrix)
    if(H == 0)
    {
        return(list(H=0, pvalue=NA))
    }
    
    l <- list()
    for (i in 1:5000)
    {
        random_sig <- beta_values_matrix[sample(nrow(beta_values_matrix), sig_length), ]
        l <- append(l, list(random_sig))
    }

    registerDoParallel(cores)
    H_random <- foreach(random_sig = l, .combine = "c") %dopar% 
    {
        H_chunk <- get_random_H(random_sig, trt)
        H_chunk
    }

    H_random <- do.call(c, as.data.frame(H_random))
    H_random <- H_random[H_random != 0]
    if(length(H_random) < 5000)
    {
        return(list(H=H, pvalue=NA))
    }
    
    # Save random values
    random_table_name<-paste(dataset_name, method, result_prefix, "H_random", sig_length, sep = "_")
    write.table(H_random, random_table_name, sep="\t", col.names=FALSE, row.names = FALSE, quote = FALSE)

    # pvalue calculation with pseudo count
    pvalue <- (sum(H_random >= H)+1)/(length(H_random)+1)
    return(list(H=H, pvalue=pvalue))
}

write_H_scores<-function(dataset, method, result_prefix, beta_values_matrix, signature, trt, cores, full_or_100) 
{
    sig_length <- length(signature)
    result <- get_H_score(dataset, method, result_prefix, beta_values_matrix, signature, trt, cores)

    # Write H score to result table
    H_full <- read.table(paste0("H_",full_or_100), header=TRUE, sep='\t')
    new_row <- data.frame(dataset, method, result_prefix, result$H, result$pvalue, sig_length)
    names(new_row) <- c("dataset", "method", "result_prefix","H","pvalue","length")
    H_full <- rbind(H_full, new_row)
    write.table(H_full, file = paste0("H_",full_or_100), sep="\t", col.names=TRUE, quote = FALSE)
}

get_limma_signature<-function(beta_values_matrix, n1, n2, trt, result_prefix, dataset, cores)
{
    limma_res <- computeDiffTab.extended.site(
        as.matrix(beta_values_matrix),
        1:n1,
        (n1+1):(n1+n2),
        diff.method = "limma",
        paired = FALSE,
        adjustment.table = NULL,
        eps = 0.01,
        imputed = FALSE)
    limma_res_fixed <- na.omit(limma_res)
    limma_res_sorted <- limma_res_fixed[order(limma_res_fixed$diffmeth.p.adj.fdr),]
    limma_res_sig <- limma_res_sorted[limma_res_sorted$diffmeth.p.adj.fdr < 0.05,]
    write.table(limma_res_sig, file = paste0(result_prefix, '_limma_full'), sep="\t", col.names=TRUE, quote = FALSE)
    which.diffmeth <- abs(limma_res_sig$mean.diff) > 0.15
    limma_res_final <- limma_res_sig[which.diffmeth, ]
    write.table(limma_res_final, file = paste0(result_prefix, '_limma_cut'), sep="\t", col.names=TRUE, quote = FALSE)
    signature <- rownames(limma_res_final)
    write_H_scores(dataset, "limma", result_prefix, beta_values_matrix, signature, trt, cores, 'full') 
    signature_100 <- head(signature, 100)
    write_H_scores(dataset, "limma", result_prefix, beta_values_matrix, signature_100, trt, cores, '100') 
    signature_10 <- head(signature, 10)
    write_H_scores(dataset, "limma", result_prefix, beta_values_matrix, signature_10, trt, cores, '10') 
}

get_ttest_signature<-function(beta_values_matrix, n1, n2, trt, result_prefix, dataset, cores)
{
    ttest_res <- computeDiffTab.extended.site(
        as.matrix(beta_values_matrix),
        1:n1,
        (n1+1):(n1+n2),
        diff.method = "ttest",
        paired = FALSE,
        adjustment.table = NULL,
        eps = 0.01,
        imputed = FALSE
    )

    ttest_res_fixed <- na.omit(ttest_res)
    ttest_res_sorted <- ttest_res_fixed[order(ttest_res_fixed$diffmeth.p.adj.fdr),]
    ttest_res_sig <- ttest_res_sorted[ttest_res_sorted$diffmeth.p.adj.fdr < 0.05,]
    write.table(ttest_res_sig, file = paste0(result_prefix, '_ttest_full'), sep="\t", col.names=TRUE, quote = FALSE)
    which.diffmeth <- abs(ttest_res_sig[, "mean.diff"]) > 0.15
    ttest_res_final <- ttest_res_sig[which.diffmeth, ]
    write.table(ttest_res_final, file = paste0(result_prefix, '_ttest_cut'), sep="\t", col.names=TRUE, quote = FALSE)
    signature <- rownames(ttest_res_final)
    write_H_scores(dataset, "ttest", result_prefix, beta_values_matrix, signature, trt, cores, 'full')
    signature_100 <- head(signature, 100)
    write_H_scores(dataset, "ttest", result_prefix, beta_values_matrix, signature_100, trt, cores, '100')
    signature_10 <- head(signature, 10)
    write_H_scores(dataset, "ttest", result_prefix, beta_values_matrix, signature_10, trt, cores, '10') 
}

get_dmpFinder_signature<-function(beta_values_matrix, n1, n2, trt, result_prefix, shrinkVar, dataset, cores)
{
    suffix = ''
    if(shrinkVar)
    {
        suffix = '_shrinkVar'
    }

    dmp <- dmpFinder(as.matrix(beta_values_matrix), pheno = trt, type = "categorical", shrinkVar = shrinkVar, qCutoff = 1)
    dmp_fixed <- na.omit(dmp)
    dmp_sorted <- dmp_fixed[order(dmp_fixed$qval),] 
    dmp_sig <- dmp_sorted[dmp_sorted$qval < 0.05,]
    write.table(dmp_sig, file = paste0(result_prefix, '_dmpFinder', suffix, '_full'), sep="\t", col.names=TRUE, quote = FALSE)
    which.diffmeth <- abs(as.vector(rowMeans(beta_values_matrix[,1:n1])) - as.vector(rowMeans(beta_values_matrix[,(n1+1):(n1+n2)]))) > 0.15
    cpgs <- rownames(beta_values_matrix[which.diffmeth,])
    signature <- rownames(dmp_sig)
    signature <- intersect(signature, cpgs)
    dmp_res_final <- dmp_sig[signature,]
    write.table(dmp_res_final, file = paste0(result_prefix, '_dmpFinder', suffix, '_cut'), sep="\t", col.names=TRUE, quote = FALSE)
    write_H_scores(dataset, paste0('dmpFinder', suffix), result_prefix, beta_values_matrix, signature, trt, cores, 'full') 
    signature_100 <- head(signature, 100)
    write_H_scores(dataset, paste0('dmpFinder', suffix), result_prefix, beta_values_matrix, signature_100, trt, cores, '100') 
    signature_10 <- head(signature, 10)
    write_H_scores(dataset, paste0('dmpFinder', suffix), result_prefix, beta_values_matrix, signature_10, trt, cores, '10') 
}

get_signatures<-function(columns1, columns2, beta_values, result_prefix, dataset, cores)
{
    beta_values_matrix <- beta_values[,c(columns1, columns2)]
    n1 <- length(columns1)
    n2 <- length(columns2)
    trt <- c(rep("group1", n1), rep("group2", n2))
    get_limma_signature(beta_values_matrix, n1, n2, trt, result_prefix, dataset, cores)
    get_ttest_signature(beta_values_matrix, n1, n2, trt, result_prefix, dataset, cores)
    get_dmpFinder_signature(beta_values_matrix, n1, n2, trt, result_prefix, FALSE, dataset, cores)
    get_dmpFinder_signature(beta_values_matrix, n1, n2, trt, result_prefix, TRUE, dataset, cores)
}

write_preprocessing_tables <- function(myImport, myLoad, result_prefix, control, case)
{
    ############### write cpg_num
    cpg_num_microarray_new<-read.table('cpg_num_microarray_new', header=TRUE, sep='\t')
    new_row<-data.frame(result_prefix, dim(myImport$beta)[1], 'before_filtering')
    names(new_row)<-c('datasets', 'cpg_num', 'filtering')
    cpg_num_microarray_new<-rbind(cpg_num_microarray_new, new_row)
    new_row<-data.frame(result_prefix, dim(myLoad$beta)[1], 'after_filtering')
    names(new_row)<-c('datasets', 'cpg_num', 'filtering')
    cpg_num_microarray_new<-rbind(cpg_num_microarray_new, new_row)
    write.table(cpg_num_microarray_new, file = 'cpg_num_microarray_new', sep="\t", col.names=TRUE, quote = FALSE)

    ############## write group size
    group_size_microarray_new<-read.table('group_size_microarray_new', header=TRUE, sep='\t')
    new_row<-data.frame(result_prefix, length(control), 'control')
    names(new_row)<-c('datasets', 'group_size', 'group')
    group_size_microarray_new<-rbind(group_size_microarray_new, new_row)
    new_row<-data.frame(result_prefix, length(case), 'case')
    names(new_row)<-c('datasets', 'group_size', 'group')
    group_size_microarray_new<-rbind(group_size_microarray_new, new_row)
    write.table(group_size_microarray_new, file = 'group_size_microarray_new', sep="\t", col.names=TRUE, quote = FALSE)
}


draw_plots<-function(dataset_name, method_name, type_sig, signature, sig_length, H_score, H_random, trt, ratio_table, dataset_title) 
{
    p <- ggplot(H_random, aes(x=H)) + 
    geom_histogram(color="darkblue", fill="lightblue", bins = 300, size = 0.05) +
    geom_point(aes(x=H_score[1], y = 0), colour = "red", size=3) +
    scale_y_continuous(limits = c(0, 5010)) +
    scale_x_continuous(limits = c(0, 1.1), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
    labs(x ='H-score', y = "Count") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title=element_text(size = 16,face="bold"),
        axis.text = element_text(size = 16)) +
    coord_cartesian(clip = "off") +
    ggtitle(paste0(dataset_title, ",\n", method_name,", length = ",sig_length))
    ggsave(paste0(dataset_name, '_', method_name, '_', type_sig, '_hist.png'),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")
    p <- NA
    heat_matrix<-ratio_table[signature,]
    signature_matrix_dist<-dist(t(heat_matrix))
    if(sum(is.na(signature_matrix_dist))>0)
    {
        print(paste0("NA: ", dataset_name, " ", method_name, " ", type_sig))
        return()
    }

    mds_fit <- cmdscale(signature_matrix_dist,eig=TRUE, k=2)
    mds_df <- as.data.frame(as.matrix(mds_fit$points))
    colnames(mds_df) <-c("V1", "V2")
    mds_df$sample_type <- trt
    p<-ggplot(mds_df, aes(x=V1, y=V2, color=sample_type)) + geom_point(size = 3) +
        labs(x ='Coordinate 1', y = "Coordinate 1", color = "Group") +
        ggtitle(paste0(dataset_title, ",\n", method_name, ", length = ",sig_length))+
        theme(plot.title = element_text(size = 20, face = "bold"),
            axis.title=element_text(size=16, face="bold"),
            axis.text = element_text(size = 16),
            legend.title=element_text(size=16, face="bold"),
            legend.text=element_text(size=16))

    ggsave(paste0(dataset_name, '_', method_name, '_', type_sig, '_mds.png'),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")
    p <- NA
}

draw_all_plots<-function(dataset, result_prefix, beta_values, control, case, dataset_title, trt, cores) 
{
    beta_values_matrix <- beta_values[,c(control, case)]
    n1 <- length(control)
    n2 <- length(case)
    trt <- c(rep("control", n1), rep("case", n2))
    means <- abs(as.vector(rowMeans(beta_values[,c(case)])) - as.vector(rowMeans(beta_values[,c(control)])))
    methods <- c("limma", "ttest", "dmpFinder", "dmpFinder_shrinkVar")
    sig_types <- c('full', '100', '10')
    for (method in methods)
    {
        signature_full <- rownames(read.csv(paste0(result_prefix, '_', method, '_cut'), sep="\t"))
        for (sig_type in sig_types)
        {
            if(sig_type == 'full')
            {
                signature = signature_full
                dataset_title_sig = paste0(dataset_title)
            }
            if(sig_type == '100')
            {
                signature = head(signature_full, 100)
                dataset_title_sig = paste0(dataset_title, ', top 100')
            }
            if(sig_type == '10')
            {
                signature = head(signature_full, 10)
                dataset_title_sig = paste0(dataset_title, ', top 10')
            }
            sig_length <- length(signature)
            if(sig_length < 2)
            {
                next
            }

            H_table <- read.csv(paste0("H_", sig_type), sep="\t")
            H_table_method <- H_table[H_table$method == method,]
            H_table_result_prefix <- H_table_method[H_table_method$result_prefix == result_prefix,]
            H <- H_table_result_prefix$H
            H_random <- read.csv(paste(dataset, method, result_prefix, "H_random", sig_length, sep = "_"), sep = '\t')
            colnames(H_random) <- c('H')
            draw_plots(result_prefix, method, sig_type, signature, sig_length, H, H_random, trt, beta_values_matrix, dataset_title_sig) 
        }
    }
}

H_full <- data.frame(dataset=character(),
    method=character(),
    result_prefix=character(),
    H=double(),
    pvalue=double(),
    length=numeric(),
    stringsAsFactors=FALSE)
write.table(H_full, file = "H_full", sep="\t", col.names=TRUE, quote = FALSE)

H_100 <- data.frame(dataset=character(),
    method=character(),
    result_prefix=character(),
    H=double(),
    pvalue=double(),
    length=numeric(),
    stringsAsFactors=FALSE)
write.table(H_100, file = "H_100", sep="\t", col.names=TRUE, quote = FALSE)

H_10 <- data.frame(dataset=character(),
    method=character(),
    result_prefix=character(),
    H=double(),
    pvalue=double(),
    length=numeric(),
    stringsAsFactors=FALSE)
write.table(H_10, file = "H_10", sep="\t", col.names=TRUE, quote = FALSE)

cpg_num_microarray_new <- data.frame(datasets=character(),
    cpg_num=character(),
    filtering=character(),
    stringsAsFactors=FALSE)
write.table(cpg_num_microarray_new, file = "cpg_num_microarray_new", sep="\t", col.names=TRUE, quote = FALSE)

group_size_microarray_new <- data.frame(datasets=character(),
    group_size=character(),
    group=character(),
    stringsAsFactors=FALSE)
write.table(group_size_microarray_new, file = "group_size_microarray_new", sep="\t", col.names=TRUE, quote = FALSE)

cores <- 16

##################### GSE210301 #####################

dataset <- 'GSE210301'
dataset_dir <- paste0(dataset, "_RAW")
gse_GSE210301 <- getGEO(dataset, GSEMatrix=TRUE)
#load(paste0('gse_',dataset))
data <- gse_GSE210301[[1]]
condition <- data$treatment
vehicle <- which(condition == "vehicle") 
cortisol <- which(condition == "cortisol")
relacorilant <- which(condition == "relacorilant")
cortisol_relacorilant <- which(condition == "cortisol and relacorilant")
myImport <- champ.import(paste0("champ/",dataset,"_RAW"),arraytype="EPIC")
myLoad <- champ.filter(beta=myImport$beta,
    pd=myImport$pd,
    detP=myImport$detP,
    beadcount=myImport$beadcount,
    autoimpute=FALSE,
    arraytype = "EPIC")

myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC",cores=cores)
beta_values <- myNorm
write.table(beta_values,paste0(dataset, '_beta_values'))

control <- vehicle
case <- cortisol
result_prefix <- paste0(dataset, "_cortisol")
beta_values <- read.csv(paste0(dataset, '_beta_values'), sep = ' ')
write_preprocessing_tables(myImport, myLoad, result_prefix, control, case)
get_signatures(control, case, beta_values, result_prefix, dataset, cores)
draw_all_plots(dataset, result_prefix, beta_values, control, case, 'GSE210301 (1)', trt, cores)

control <- vehicle
case <- relacorilant
result_prefix <- paste0(dataset, "_relacorilant")
write_preprocessing_tables(myImport, myLoad, result_prefix, control, case)
get_signatures(control, case, beta_values, result_prefix, dataset, cores)
draw_all_plots(dataset, result_prefix, beta_values, control, case, 'GSE210301 (2)', trt, cores) 

control <- vehicle
case <- cortisol_relacorilant
result_prefix <- paste0(dataset, "_cortisol_relacorilant")
write_preprocessing_tables(myImport, myLoad, result_prefix, control, case)
get_signatures(control, case, beta_values, result_prefix, dataset, cores)
draw_all_plots(dataset, result_prefix, beta_values, control, case, 'GSE210301 (3)', trt, cores) 

myImport <- NA
myLoad <- NA
myNorm <- NA

##################### GSE175458 #####################

dataset <- 'GSE175458'
dataset_dir <- paste0(dataset, "_RAW")
annotation_file <- paste0(dataset_dir,'//', dataset, "_files.csv")
gse <- getGEO(dataset, GSEMatrix=TRUE, getGPL=FALSE)
#load(paste0('gse_',dataset))
data <- gse_GSE175458[[1]]
condition <- data$case_control_status
myImport <- champ.import(paste0("champ/",dataset,"_RAW"),arraytype="EPIC")
myLoad <- champ.filter(beta=myImport$beta,
    pd=myImport$pd,
    detP=myImport$detP,
    beadcount=myImport$beadcount,
    autoimpute=FALSE,
    arraytype = "EPIC")

myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC",cores=cores)
beta_values <- myNorm
write.table(beta_values,paste0(dataset, '_beta_values'))
beta_values <- read.csv(paste0(dataset, '_beta_values'), sep = ' ')
control <- which(condition == "control") 
case <- which(condition == "case")
result_prefix <- dataset

write_preprocessing_tables(myImport, myLoad, result_prefix, control, case)
get_signatures(control, case, beta_values, result_prefix, dataset, cores)
draw_all_plots(dataset, result_prefix, beta_values, control, case, 'GSE175458', trt, cores) 

myImport <- NA
myLoad <- NA
myNorm <- NA

##################### GSE156994 #####################

dataset <- 'GSE156994'
dataset_dir <- paste0(dataset, "_RAW")
annotation_file <- paste0(dataset_dir,'//', dataset, "_files.csv")
gse <- getGEO(dataset, GSEMatrix=TRUE)
#load('gse_GSE156994')
data <- gse_GSE156994[[1]]
condition <- data$`sample_group`
myImport <- champ.import(paste0("champ/",dataset,"_RAW"),arraytype="450K")
myLoad <- champ.filter(beta=myImport$beta,
    pd=myImport$pd,
    detP=myImport$detP,
    beadcount=myImport$beadcount,
    autoimpute=FALSE,
    arraytype = "450K")

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=cores)
beta_values <- myNorm
write.table(beta_values,paste0(dataset, '_beta_values'))
beta_values <- read.table(paste0(dataset, '_beta_values'))
control <- which(condition ==  "CTRL") 
case <- which(condition == "sCJD")
result_prefix <- dataset

write_preprocessing_tables(myImport, myLoad, result_prefix, control, case)
get_signatures(control, case, beta_values, result_prefix, dataset, cores)
draw_all_plots(dataset, result_prefix, beta_values, control, case, 'GSE156994', trt, cores) 

myImport <- NA
myLoad <- NA
myNorm <- NA

##################### GSE157341 #####################

dataset <- 'GSE157341'
dataset_dir <- paste0(dataset, "_RAW")
annotation_file <- paste0(dataset_dir,'//', dataset, "_files.csv")
gse <- getGEO(dataset, GSEMatrix=TRUE)
#load('gse_GSE157341')
data <- gse_GSE157341[[1]]
condition <- data$disease
myImport <- champ.import(paste0("champ/",dataset,"_RAW"),arraytype="450K")
myLoad <- champ.filter(beta=myImport$beta,
    pd=myImport$pd,
    detP=myImport$detP,
    beadcount=myImport$beadcount,
    autoimpute=FALSE,
    arraytype = "450K")

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=cores)
beta_values <- myNorm
write.table(beta_values,paste0(dataset, '_beta_values'))
beta_values <- read.table(paste0(dataset, '_beta_values'))
control <- which(condition == "non-tumor liver tissue") 
case <- which(condition == "hepatocellular carcinoma") 
result_prefix <- dataset

write_preprocessing_tables(myImport, myLoad, result_prefix, control, case)
get_signatures(control, case, beta_values, result_prefix, dataset, cores)
draw_all_plots(dataset, result_prefix, beta_values, control, case, 'GSE157341', trt, cores) 

myImport <- NA
myLoad <- NA
myNorm <- NA

##################### GSE101764 #####################

dataset <- 'GSE101764'
dataset_dir <- paste0(dataset, "_RAW")
annotation_file <- paste0(dataset_dir,'//', dataset, "_files.csv")
gse <- getGEO(dataset, GSEMatrix=TRUE)
#load('gse_GSE101764')
data <- gse_GSE101764[[1]]
condition <- data$tissue
myImport <- champ.import(paste0("champ/",dataset,"_RAW"),arraytype="450K")
myLoad <- champ.filter(beta=myImport$beta,
    pd=myImport$pd,
    detP=myImport$detP,
    beadcount=myImport$beadcount,
    autoimpute=FALSE,
    arraytype = "450K")

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=cores)
beta_values <- myNorm
write.table(beta_values,paste0(dataset, '_beta_values'))
beta_values <- read.table(paste0(dataset, '_beta_values'))
control <- which(condition == "mucosa") 
case <- which(condition == "tumor") 
result_prefix <- dataset

write_preprocessing_tables(myImport, myLoad, result_prefix, control, case)
get_signatures(control, case, beta_values, result_prefix, dataset, cores)
draw_all_plots(dataset, result_prefix, beta_values, control, case, 'GSE101764', trt, cores) 

myImport <- NA
myLoad <- NA
myNorm <- NA

##################### GSE85845 #####################

dataset <- 'GSE85845'
dataset_dir <- paste0(dataset, "_RAW")
annotation_file <- paste0(dataset_dir,'//', dataset, "_files.csv")
gse <- getGEO(dataset, GSEMatrix=TRUE)
#load('gse_GSE85845')
data <- gse_GSE85845[[1]]
condition <- data$`disease state`
myImport <- champ.import(paste0("champ/",dataset,"_RAW"),arraytype="450K")
myLoad <- champ.filter(beta=myImport$beta,
    pd=myImport$pd,
    detP=myImport$detP,
    beadcount=myImport$beadcount,
    autoimpute=FALSE,
    arraytype = "450K")

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=cores)
beta_values <- myNorm
write.table(beta_values,paste0(dataset, '_beta_values'))
beta_values <- read.table(paste0(dataset, '_beta_values'))
control <- which(condition == "non-tumor") 
case <- which(condition == "lung adenocarcinoma") 
result_prefix <- dataset

write_preprocessing_tables(myImport, myLoad, result_prefix, control, case)
get_signatures(control, case, beta_values, result_prefix, dataset, cores)
draw_all_plots(dataset, result_prefix, beta_values, control, case, 'GSE85845', trt, cores) 

myImport <- NA
myLoad <- NA
myNorm <- NA

##################### GSE175399 #####################

dataset <- 'GSE175399'
dataset_dir <- paste0(dataset, "_RAW")
annotation_file <- paste0(dataset_dir,'//', dataset, "_files.csv")
gse_GSE175399 <- getGEO(dataset, GSEMatrix=TRUE)
#load('gse_GSE175399')
data <- gse_GSE175399[[1]]
condition <- data$disease
myImport <- champ.import(paste0("champ/",dataset,"_RAW"),arraytype="EPIC")
myLoad <- champ.filter(beta=myImport$beta,
    pd=myImport$pd,
    detP=myImport$detP,
    beadcount=myImport$beadcount,
    autoimpute=FALSE,
    arraytype = "EPIC")

myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC",cores=cores)
beta_values <- myNorm
write.table(beta_values,paste0(dataset, '_beta_values'))
beta_values <- read.table(paste0(dataset, '_beta_values'))
control <- which(condition == "normal") 
case <- which(condition == "thyroid-associated ophthalmopathy") 
result_prefix <- dataset

write_preprocessing_tables(myImport, myLoad, result_prefix, control, case)
get_signatures(control, case, beta_values, result_prefix, dataset, cores)
draw_all_plots(dataset, result_prefix, beta_values, control, case, 'GSE175399', trt, cores) 

myImport <- NA
myLoad <- NA
myNorm <- NA

##################### GSE210484 #####################

dataset <- 'GSE210484'
dataset_dir <- paste0(dataset, "_RAW")
annotation_file <- paste0(dataset_dir,'//', dataset, "_files.csv")
gse_GSE210484 <- getGEO(dataset, GSEMatrix=TRUE)
#load('gse_GSE210484')
data <- gse_GSE210484[[1]]
condition <- data$`disease state`
myImport <- champ.import(paste0("champ/",dataset,"_RAW"),arraytype="EPIC")
myLoad <- champ.filter(beta=myImport$beta,
    pd=myImport$pd,
    detP=myImport$detP,
    beadcount=myImport$beadcount,
    autoimpute=FALSE,
    arraytype = "EPIC")

myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC",cores=cores)
beta_values <- myNorm
write.table(beta_values,paste0(dataset, '_beta_values'))
beta_values <- read.table(paste0(dataset, '_beta_values'))
control <- which(condition == "CONTROL")
case <- which(condition == "Arboleda-Tham syndrome (ARTHS)")
result_prefix <- dataset

write_preprocessing_tables(myImport, myLoad, result_prefix, control, case)
get_signatures(control, case, beta_values, result_prefix, dataset, cores)
draw_all_plots(dataset, result_prefix, beta_values, control, case, 'GSE210484', trt, cores) 

myImport <- NA
myLoad <- NA
myNorm <- NA

##################### GSE196007 #####################

dataset <- 'GSE196007'
dataset_dir <- paste0(dataset, "_RAW")
annotation_file <- paste0(dataset_dir,'//', dataset, "_files.csv")
gse_GSE196007 <- getGEO(dataset, GSEMatrix=TRUE)
#load('gse_GSE196007')
data <- gse_GSE196007[[1]]
condition <- data$`disease status`
myImport <- champ.import(paste0("champ/",dataset,"_RAW"),arraytype="EPIC")
dim(myImport$beta)
myLoad <- champ.filter(beta=myImport$beta,
    pd=myImport$pd,
    detP=myImport$detP,
    beadcount=myImport$beadcount,
    autoimpute=FALSE,
    arraytype = "EPIC")

myNorm <- champ.norm(beta=myLoad$beta,arraytype="EPIC",cores=cores)
beta_values <- myNorm
write.table(beta_values,paste0(dataset, '_beta_values'))
beta_values <- read.table(paste0(dataset, '_beta_values'))
control <- which(condition == "Control") 
case <- which(condition == "Systemic sclerosis (SSc)") 
result_prefix <- dataset

write_preprocessing_tables(myImport, myLoad, result_prefix, control, case)
get_signatures(control, case, beta_values, result_prefix, dataset, cores)
draw_all_plots(dataset, result_prefix, beta_values, control, case, 'GSE196007', trt, cores) 

myImport <- NA
myLoad <- NA
myNorm <- NA

##################### GSE156669 #####################

dataset <- 'GSE156669'
dataset_dir <- paste0(dataset, "_RAW")
annotation_file <- paste0(dataset_dir,'//', dataset, "_files.csv")
gse_GSE156669 <- getGEO(dataset, GSEMatrix=TRUE)
#load('gse_GSE156669')
data <- gse_GSE156669[[1]]
condition <- data$tissue
myImport <- champ.import(paste0("champ/",dataset,"_RAW"),arraytype="450K")
myLoad <- champ.filter(beta=myImport$beta,
    pd=myImport$pd,
    detP=myImport$detP,
    beadcount=myImport$beadcount,
    autoimpute=FALSE,
    arraytype = "450K")

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=cores)
beta_values <- myNorm
write.table(beta_values,paste0(dataset, '_beta_values'))
beta_values <- read.table(paste0(dataset, '_beta_values'))
control <- which(condition == "normal buccal mucosa") 
case <- which(condition == "Oral Submucous fibrosis") 
result_prefix <- dataset

write_preprocessing_tables(myImport, myLoad, result_prefix, control, case)
get_signatures(control, case, beta_values, result_prefix, dataset, cores)
draw_all_plots(dataset, result_prefix, beta_values, control, case, 'GSE156669', trt, cores) 

myImport <- NA
myLoad <- NA
myNorm <- NA

##################### GSE178218 #####################

dataset <- 'GSE178218'
dataset_dir <- paste0(dataset, "_RAW")
annotation_file <- paste0(dataset_dir,'//', dataset, "_files.csv")
gse_GSE178218 <- getGEO(dataset, GSEMatrix=TRUE)
#load('gse_GSE178218')
data <- gse_GSE178218[[1]]
condition <- data$tissue
myImport <- champ.import(paste0("champ/",dataset,"_RAW"),arraytype="450K")
myLoad <- champ.filter(beta=myImport$beta,
    pd=myImport$pd,
    detP=myImport$detP,
    beadcount=myImport$beadcount,
    autoimpute=FALSE,
    arraytype = "450K")

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=cores)
beta_values <- myNorm
write.table(beta_values,paste0(dataset, '_beta_values'))
beta_values <- read.table(paste0(dataset, '_beta_values'))
control <- which(condition == "adjacent non-tumor") 
case <- which(condition == "tumor") 

result_prefix <- dataset
write_preprocessing_tables(myImport, myLoad, result_prefix, control, case)
get_signatures(control, case, beta_values, result_prefix, dataset, cores)
draw_all_plots(dataset, result_prefix, beta_values, control, case, 'GSE178218', trt, cores) 

myImport <- NA
myLoad <- NA
myNorm <- NA

##################### GSE178216 #####################

dataset <- 'GSE178216'
dataset_dir <- paste0(dataset, "_RAW")
annotation_file <- paste0(dataset_dir,'//', dataset, "_files.csv")
gse_GSE178216 <- getGEO(dataset, GSEMatrix=TRUE)
#load('gse_GSE178216')
data <- gse_GSE178216[[1]]
condition <- data$tissue
myImport <- champ.import(paste0("champ/",dataset,"_RAW"),arraytype="450K")
myLoad <- champ.filter(beta=myImport$beta,
    pd=myImport$pd,
    detP=myImport$detP,
    beadcount=myImport$beadcount,
    autoimpute=FALSE,
    arraytype = "450K")

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=cores)
beta_values <- myNorm
write.table(beta_values,paste0(dataset, '_beta_values'))
beta_values <- read.table(paste0(dataset, '_beta_values'))
control <- which(condition == "adjacent non-tumor") 
case <- which(condition == "tumor") 
result_prefix <- dataset

write_preprocessing_tables(myImport, myLoad, result_prefix, control, case)
get_signatures(control, case, beta_values, result_prefix, dataset, cores)
draw_all_plots(dataset, result_prefix, beta_values, control, case, 'GSE178216', trt, cores) 

myImport <- NA
myLoad <- NA
myNorm <- NA

##################### GSE178212 #####################

dataset <- 'GSE178212'
dataset_dir <- paste0(dataset, "_RAW")
annotation_file <- paste0(dataset_dir,'//', dataset, "_files.csv")
gse_GSE178212 <- getGEO(dataset, GSEMatrix=TRUE)
#load('gse_GSE178212')
data <- gse_GSE178212[[1]]
condition <- data$tissue
myImport <- champ.import(paste0("champ/",dataset,"_RAW"),arraytype="450K")
myLoad <- champ.filter(beta=myImport$beta,
    pd=myImport$pd,
    detP=myImport$detP,
    beadcount=myImport$beadcount,
    autoimpute=FALSE,
    arraytype = "450K")

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=cores)
beta_values <- myNorm
write.table(beta_values,paste0(dataset, '_beta_values'))
beta_values <- read.table(paste0(dataset, '_beta_values'))
control <- which(condition == "adjacent non-tumor") 
case <- which(condition == "tumor") 

result_prefix <- dataset
write_preprocessing_tables(myImport, myLoad, result_prefix, control, case)
get_signatures(control, case, beta_values, result_prefix, dataset, cores)
draw_all_plots(dataset, result_prefix, beta_values, control, case, 'GSE178212', trt, cores) 

myImport <- NA
myLoad <- NA
myNorm <- NA

##################### GSE157272 #####################

dataset <- 'GSE157272'
dataset_dir <- paste0(dataset, "_RAW")
annotation_file <- paste0(dataset_dir,'//', dataset, "_files.csv")
gse_GSE157272 <- getGEO(dataset, GSEMatrix=TRUE)
#load('gse_GSE157272')
data <- gse_GSE157272[[1]]
condition <- data$`disease state`
myImport <- champ.import(paste0("champ/",dataset,"_RAW"),arraytype="450K")
myLoad <- champ.filter(beta=myImport$beta,
    pd=myImport$pd,
    detP=myImport$detP,
    beadcount=myImport$beadcount,
    autoimpute=FALSE,
    arraytype = "450K")

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=cores)
beta_values <- myNorm
write.table(beta_values,paste0(dataset, '_beta_values'))
beta_values <- read.table(paste0(dataset, '_beta_values'))
control <- which(condition == "benign prostate tissue") 
case <- which(condition == "agressive prostate cancer tissue") 
result_prefix <- paste0(dataset)

write_preprocessing_tables(myImport, myLoad, result_prefix, control, case)
get_signatures(control, case, beta_values, result_prefix, dataset, cores)
draw_all_plots(dataset, result_prefix, beta_values, control, case, 'GSE157272', trt, cores) 
