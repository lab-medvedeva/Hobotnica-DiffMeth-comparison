library(GEOquery)
library(Hobotnica)
library(foreach)
library(doParallel)
library (ggplot2)
library(minfi)
library(IlluminaHumanMethylation450kmanifest)

get_H<-function(signature_matrix, trt, sig_length)
{
  if (sig_length < 2)
    {
        return(0.5)
    }

    dist_matrix<-dist(t(signature_matrix))
    return(Hobotnica(dist_matrix, trt))
}

get_random_H <-function(H, sig, trt, data, sig_length, cores) 
{
    if (sig_length < 2)
    {
        return(list(H=c(), pvalue=NA))
    }
    
    # P-value calculation
    registerDoParallel(cores)
    H_random<-foreach(i = 1:100000, .combine = "c") %dopar% 
    {   
      random_signature_matrix<-data[sample(nrow(data), sig_length), ]
        H_chunk<-get_H(random_signature_matrix, trt, sig_length)
        H_chunk
    }

    # pvalue calculation with pseudo count
    pvalue<-(sum(H_random > H)+1)/(length(H_random)+1)
    return(list(H=H_random, pvalue=pvalue))
}

signature_df <- read.table("signatures//M-CLL_U-CLL.csv", header = F, sep = ' ')
signature <- signature_df$V1


gse <- getGEO('GSE144894', GSEMatrix=TRUE)
data <- gse[[1]]

ighv_mutation_status <- data$`ighv:ch1`
data_matrix <- exprs(data)

# truncate signature
signature_truncated <-intersect(signature, rownames(data_matrix) )
signature_matrix <- data_matrix[signature_truncated,]
signature_length <- length(signature_truncated)
H_score_untreated <- get_H(signature_matrix, ighv_mutation_status, signature_length)
result_untreated <- get_random_H(H_score_untreated, signature_truncated, ighv_mutation_status, data_matrix, signature_length, 120)
result_H_df <- as.data.frame(result_untreated$H)
colnames(result_H_df) <- c('H')

p <- ggplot(result_H_df, aes(x=H)) + 
  geom_histogram(color="darkblue", fill="lightblue", bins = 100) +
  geom_point(aes(x=H_score_untreated, y = 0), colour = "red", size=3) +
  scale_y_continuous(limits = c(0, 40000)) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x ='H-score', y = "Count") +
  ggtitle("H-score distribution for GSE144894, Kulis et al. (2012) signature")+
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(size = 10)) +
  coord_cartesian(clip = "off")
  
ggsave('plots//GSE144894_validation.png',  plot = p,  device = "png", width = 20,  height = 15,  units = "cm")

write.table(result_H_df,
  "plots//GSE144894_H.csv" ,
  quote=F,
  row.names=F,
  col.names = T,
  sep = '\t')

signature_matrix_dist <- dist(t(signature_matrix))
mds_fit <- cmdscale(signature_matrix_dist,eig=TRUE, k=2) # k is the number of dim

ighv_mutation_status <- replace(ighv_mutation_status, ighv_mutation_status=='M', 'mutated')
ighv_mutation_status <- replace(ighv_mutation_status, ighv_mutation_status=='U', 'unmutated')

mds_df <- as.data.frame(as.matrix(mds_fit$points))
colnames(mds_df) <-c("V1", "V2")
mds_df$sample_type <- ighv_mutation_status
p<-ggplot(mds_df, aes(x=V1, y=V2, color=sample_type)) + geom_point(size = 3) +
  labs(x ='Coordinate 1', y = "Coordinate 1", color = "Group") +
  ggtitle("MDS for GSE144894, Kulis et al. (2012) signature")+
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(size = 10),
        legend.title=element_text(size=10,face="bold"),
        legend.text=element_text(size=10))

ggsave(paste0("plots//GSE144894_validation_MDS.png"),  plot = p,  device = "png", width = 20,  height = 15,  units = "cm")
