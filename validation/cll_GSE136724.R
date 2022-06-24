library(GEOquery)
library(Hobotnica)
library(foreach)
library(doParallel)
library (ggplot2)

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
gse <- getGEO('GSE136724', GSEMatrix=TRUE)
data <- gse[[1]]
data_matrix <- exprs(data)

ighv_mutation_status <- data$`ighv_mutation_status:ch1`
untreated <- which(data$`treated:ch1` == "no") 
treated <- which(data$`treated:ch1` == "yes") 
ighv_mutation_status_untreated <- ighv_mutation_status[untreated]
data_matrix_untreated <- data_matrix[,untreated] #exprs(data)[,untreated]

# truncate signature
signature_truncated <-intersect(signature, rownames(data_matrix_untreated) )
signature_matrix <- data_matrix_untreated[signature_truncated,]
signature_length <- length(signature_truncated)
H_score_untreated <- get_H(signature_matrix, ighv_mutation_status_untreated, signature_length)
result_untreated <- get_random_H(H_score_untreated, signature_truncated, ighv_mutation_status_untreated, data_matrix_untreated, signature_length, 120)
result_H_df <- as.data.frame(result_untreated$H)

colnames(result_H_df) <- c('H')

p <- ggplot(result_H_df, aes(x=H)) + 
  geom_histogram(color="darkblue", fill="lightblue", bins = 100) +
  geom_point(aes(x=H_score_untreated, y = 0), colour = "red", size=3) +
  scale_y_continuous(limits = c(0, 40000)) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x ='H-score', y = "Count") +
  ggtitle("H-score distribution for GSE136724, Kulis et al. (2012) signature")+
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(size = 10)) +
  coord_cartesian(clip = "off")

ggsave('plots//GSE136724_untreated_validation.png',  plot = p,  device = "png", width = 20,  height = 15,  units = "cm")

write.table(result_H_df,
  "plots//GSE136724_H.csv" ,
  quote=F,
  row.names=F,
  col.names = T,
  sep = '\t')

signature_matrix_dist <- dist(t(signature_matrix))
mds_fit <- cmdscale(signature_matrix_dist,eig=TRUE, k=2) # k is the number of dim

mds_df <- as.data.frame(as.matrix(mds_fit$points))
colnames(mds_df) <-c("V1", "V2")
mds_df$sample_type <- ighv_mutation_status_untreated
p<-ggplot(mds_df, aes(x=V1, y=V2, color=sample_type)) + geom_point(size = 3) +
  labs(x ='Coordinate 1', y = "Coordinate 1", color = "Group") +
  ggtitle("MDS for GSE136724, Kulis et al. (2012) signature")+
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title=element_text(size=10,face="bold"),
        axis.text = element_text(size = 10),
        legend.title=element_text(size=10,face="bold"),
        legend.text=element_text(size=10))


ggsave(paste0("plots//GSE136724_untreated_validation_MDS.png"),  plot = p,  device = "png", width = 20,  height = 15,  units = "cm")

write.table(mds_df,
	"plots//GSE136724_untreated_validation_MDS.csv",
	sep = '\t',
	quote = F,
	row.names = F)

