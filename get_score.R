library(foreach)
library(doParallel)
library(Hobotnica)

get_random_H<-function(sig_length, trt, ratio_table,length_per_process)
{
    random_sig<-ratio_table[sample(nrow(ratio_table), sig_length), ]
    distMatrix<-dist(t(random_sig))
    return(Hobotnica(distMatrix, trt))
}

get_H_score<-function(dataset_name, method, signature, trt, cores) 
{
    sig_length<-length(signature)

    if (sig_length < 2)
    {
        return(list(H=0.5, pvalue=NA))
    }
    
    ratio_table<-read.table(file = "ratio_table", sep="\t")
    heat_matrix<-ratio_table[signature,]
    distMatrix<-dist(t(heat_matrix))
    H<-Hobotnica(distMatrix, trt)

    # P-value calculation
    registerDoParallel(cores)
    H_random<-foreach(i = 1:10000, .combine = "c") %dopar% 
    {
        H_chunk<-get_random_H(sig_length, trt, ratio_table)
        H_chunk
    }
    H_random<-do.call(c, as.data.frame(H_random))

    # Save random values
    random_table_name<-paste(dataset_name, method, "H_random", sig_length, sep = "_")
    write.table(H_random, random_table_name, sep="\t", col.names=FALSE, row.names = FALSE, quote = FALSE)

    # pvalue calculation with pseudo count
    pvalue<-(sum(H_random > H)+1)/(length(H_random)+1)
    return(list(H=H, pvalue=pvalue))
}