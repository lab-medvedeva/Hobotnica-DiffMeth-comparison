library (ggplot2)
library(ggpattern)
library(ggpubr)
library(venn)
library(ggbeeswarm)
library(GEOquery)

c_datasets <- c('GSE210301_cortisol', 'GSE210301_relacorilant', 'GSE210301_cortisol_relacorilant', 'GSE175458', 'GSE175399', 'GSE210484', 'GSE196007', 'GSE156994', 'GSE157341', 'GSE101764', 'GSE85845', 'GSE156669', 'GSE178218', 'GSE178216', 'GSE178212', 'GSE157272')
dataset_names <- c('GSE210301_cortisol' = 'GSE210301 (1)',
                   'GSE210301_relacorilant' = 'GSE210301 (2)',
                   'GSE210301_cortisol_relacorilant' = 'GSE210301 (3)',
                   'GSE175458' = 'GSE175458',
                   'GSE175399' = 'GSE175399',
                   'GSE210484' = 'GSE210484',
                   'GSE196007' = 'GSE196007',
                   'GSE156994' = 'GSE156994',
                   'GSE157341' = 'GSE157341',
                   'GSE101764' = 'GSE101764',
                   'GSE85845' = 'GSE85845',
                   'GSE156669' = 'GSE156669',
                   'GSE178218' = 'GSE178218',
                   'GSE178216' = 'GSE178216',
                   'GSE178212' = 'GSE178212',
                   'GSE157272' = 'GSE157272')

c_datasets_names <- c('GSE210301 (1)', 'GSE210301 (2)', 'GSE210301 (3)', 'GSE175458', 'GSE175399', 'GSE210484', 'GSE196007', 'GSE156994', 'GSE157341', 'GSE101764', 'GSE85845', 'GSE156669', 'GSE178218', 'GSE178216', 'GSE178212', 'GSE157272')
c_methods_ids <- c("limma", "ttest", "dmpFinder", "dmpFinder_shrinkVar")
c_methods <-c("limma","T-test", "dmpFinder", "dmpFinder, vs")
method_names <- c(
  'limma' = "limma",
  'ttest' = "T-test",
  'dmpFinder' = "dmpFinder",
  'dmpFinder_shrinkVar' = "dmpFinder, vs")

q20 <-  c(
  RColorBrewer::brewer.pal(8,'Dark2'),
  RColorBrewer::brewer.pal(12,'Paired')
)

H_all_100 <- read.table("H_all_100",sep='\t', header = T)
H_all_100$length <- as.numeric(H_all_100$length)
H_all_100$result_prefix<-factor(H_all_100$result_prefix,levels = c_datasets)
H_all_100$method<-factor(H_all_100$method, levels = c_methods_ids)
H_all_100$tag<-factor(H_all_100$tag, levels = c('H_all', 'H_100', 'H_10'))

tag_names <- c(
  'H_all' = "Full",
  'H_100' = "Top 100 DMP",
  'H_10' = "Top 10 DMP")

c_tag_names <- c("Full",
  "Top 100 DMP",
  "Top 10 DMP")

################### Example Panel ################### 
###### Heatmap example
beta_values <- read.csv(paste0(dataset, 'beta_values'), sep = ' ')
limma <- read.csv('GSE178216_limma_cut', sep='\t')
signature <- rownames(limma)
heat_matrix<-beta_values[signature,]
distMatrix<-dist(t(heat_matrix))
distMatrix_df <- data.frame(as.matrix(distMatrix))
distMatrix_df_flat <- data.frame(x=character(),
                                 x_type=character(),
                                 y=character(),
                                 y_type=character(),
                                 d=double(),
                                 stringsAsFactors=FALSE)
for(r in 1:nrow(distMatrix_df))
{
  for(c in 1:ncol(distMatrix_df))
  {
    if(r %in% control)
    {
      x_type = 'control'
    }
    else
    {
      x_type = 'case'
    }
    
    if(c %in% control)
    {
      y_type = 'control'
    }
    else
    {
      y_type = 'case'
    }
    
    d <- distMatrix_df[r,c]
    new_row<-data.frame(rownames(distMatrix_df)[r], x_type, colnames(distMatrix_df)[c], y_type, d)
    names(new_row) <- c("x","x_type","y","y_type","distance")
    distMatrix_df_flat <- rbind(distMatrix_df_flat, new_row)
  }
}
write.table(distMatrix_df_flat, file = paste(dataset, 'limma', "distmatrix", sep = "_"), sep="\t", col.names=TRUE, quote = FALSE)

dataset <- 'GSE178216'
load('gse_GSE178216')
data <- gse_GSE178216[[1]]
condition <- data$tissue

condition <- c('adjacent non-tumor', 'tumor', 'tumor', 'adjacent non-tumor', 'tumor', 'adjacent non-tumor', 'tumor', 'adjacent non-tumor', 'tumor', 'tumor', 'tumor', 'tumor', 'tumor', 'adjacent non-tumor', 'tumor', 'tumor', 'adjacent non-tumor', 'tumor', 'adjacent non-tumor', 'tumor', 'tumor', 'tumor')
control <- which(condition == "adjacent non-tumor") 
case <- which(condition == "tumor") 

dist_matrix_df<-read.csv("GSE178216_limma_distmatrix", sep='\t')

colnames(dist_matrix_df) <- c('x','x_type','y','y_type','distance')

ord <- c(paste0('a_',1:length(control)), paste0('t_',1:length(case)))
dist_matrix_df$x<-factor(dist_matrix_df$x, levels = ord)
dist_matrix_df$y<-factor(dist_matrix_df$y, levels = ord)
dist_matrix_df$x_type<-factor(dist_matrix_df$x_type, levels = c('control','case'))
dist_matrix_df$y_type<-factor(dist_matrix_df$y_type, levels = c('case','control'))

p <- ggplot(dist_matrix_df, aes(x, y)) + 
  facet_grid(y_type~x_type, scales = "free",space="free") +
  geom_tile(aes(fill= distance)) +
  scale_fill_distiller(palette = "RdBu", direction = 1)+  #RdYlBu
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text = element_text(size= 16), #, face = "bold"
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.key.size = unit(0.6, "cm"),
        legend.key.width =  unit(0.8, "cm"),
        legend.position = "bottom",
        legend.key=element_rect(colour="black"),
        legend.margin = margin(t = -0.3, unit='cm') )


ggsave(paste0("GSE178216_example_heatmap.png"),  plot = p,  device = "png", width = 15,  height = 15,  units = "cm")

###### Hist example
H_GSE178216 <- H_all_100[H_all_100$result_prefix == 'GSE178216',]
H_GSE178216$length <- as.numeric(H_GSE178216$length)
H_GSE178216$result_prefix<-factor(H_GSE178216$result_prefix,levels = c_datasets)
H_GSE178216$method<-factor(H_GSE178216$method, levels = c_methods_ids)
H_GSE178216$tag<-factor(H_GSE178216$tag, levels = c('H_all', 'H_100', 'H_10'))

H_point<- H_GSE178216[H_GSE178216$method == 'limma' & H_GSE178216$tag == 'H_all',]$H
tag_names <- c(
  'H_all' = "Full DMP signatures",
  'H_100' = "Top 100 DMP signatures",
  'H_10' = "Top 10 DMP signatures")

H_random <- read.csv('GSE178216_limma_GSE178216_H_random_38878', header = F)
colnames(H_random) <-c('H')
p <- ggplot(H_random, aes(x=H)) + 
  geom_histogram(color="darkblue", fill="lightblue", bins = 300, size = 0.05) +
  geom_point(aes(x=H_point, y = 0), colour = "red", size=2) +
  scale_y_continuous(limits = c(0, 2000)) +
  scale_x_continuous(limits = c(0.5, 0.9), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  labs(x ='H-score', y = "Count") +
  theme_bw() +
  theme(#panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 15), #, face = "bold"
        axis.title=element_text(size = 15), #,face="bold"
        axis.text = element_text(size = 15)) +
  coord_cartesian(clip = "off")

ggsave("GSE178216_example_hist.png",  plot = p,  device = "png", width = 15,  height = 8,  units = "cm")

###### Barplot example
p <- ggplot(H_GSE178216, aes(x=method, y = H, fill=method, pattern = tag))+ 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        text = element_text(size=14), #, face = "bold"
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(t = -0.3, unit='cm'))+
  geom_bar_pattern(stat = "identity",
                   position='dodge',
                   aes(pattern  = tag),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.04,
                   pattern_key_scale_factor = 0.6)+
  guides(pattern = guide_legend(nrow = 3,override.aes = list(fill = "white")),
         fill = guide_legend(ncol=1, override.aes = list(pattern = "none"))) +
  labs(x ='Method', y = "H", pattern = "Signature type") + 
  scale_fill_brewer(name = "Method", labels = c_methods, palette="Set1", drop = F)+ 
  scale_pattern_manual(values = c('none', 'stripe', 'circle'),
                       labels = c_tag_names)+
  scale_x_discrete(drop=FALSE)

ggsave(paste0("GSE178216_example_barplot.png"),  plot = p,  device = "png", width = 15,  height = 8,  units = "cm")

################### Datasets ###################
p <- ggplot(H_all_100,aes(x=method, y = H, fill=method, pattern = tag))+
  facet_wrap(~ result_prefix,  ncol = 4, labeller = as_labeller(dataset_names))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        text = element_text(size=12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(t = -0.3, unit='cm'))+
  geom_bar_pattern(stat = "identity",
                   position='dodge',
                   aes(pattern  = tag),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.03,
                   pattern_key_scale_factor = 0.6)+
  guides(pattern = guide_legend(nrow = 3,override.aes = list(fill = "white")),
         fill = guide_legend(ncol=2, override.aes = list(pattern = "none"))) +
  labs(x ='Method', y = "H", pattern = "Signature type") + 
  scale_fill_brewer(name = "Method", labels = c_methods, palette="Set1", drop = F)+ 
  scale_pattern_manual(values = c('none', 'stripe', 'circle'),
                       labels = c_tag_names)+
  scale_x_discrete(drop=FALSE)

ggsave(paste0("datasets_microarray_all.png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")

################### Boxplots ###################
H_all_100_with_zero <- H_all_100
H_all_100_with_zero$is_H_zero <- TRUE
H_all_100_without_zero <- H_all_100
H_all_100_without_zero <- H_all_100_without_zero[H_all_100_without_zero$H!=0,]
H_all_100_without_zero$is_H_zero <- FALSE
  
data<- rbind(H_all_100_with_zero, H_all_100_without_zero )
data <- data[data$tag != 'H_10',]
data$is_H_zero <- factor(data$is_H_zero, levels = c(TRUE, FALSE))
data$tag <- factor(data$tag, levels = c('H_all', 'H_100'))
data$method <- factor(data$method, levels = c_methods_ids)
    
p<-ggplot(data, aes(x=method, y=H, fill=method, alpha = is_H_zero )) + 
  facet_wrap(~ tag,  ncol = 2, labeller = as_labeller(tag_names)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        legend.position="bottom",
        text = element_text(size=20)) + #, face = "bold"
  geom_boxplot() +  labs(x ='Method', y = "H") + 
  scale_y_continuous(limits = c(0, 1))+
  scale_x_discrete( drop=FALSE) + 
  scale_fill_brewer(name = "Method", labels = c_methods,palette="Set1",drop = F) +
  scale_alpha_manual(name = "Include H=0",values=c(1,0.3)) +
  guides(alpha=guide_legend(ncol=1, nrow = 2, override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(1,0.3)),
                                              colour=NA)),
         fill = guide_legend(ncol=2))
  
ggsave(paste0("methods_boxplot_microarray.png"),  plot = p,  device = "png", width = 30,  height = 10,  units = "cm")

################### Shapiro test ################### 
H_all <- H_all_100[H_all_100$tag == 'H_all',]
H_100 <- H_all_100[H_all_100$tag == 'H_100',]
H_10 <- H_all_100[H_all_100$tag == 'H_10',]

shapiro.test(H_all$H)
shapiro.test(H_100$H)
shapiro.test(H_10$H)

################### Friedman test ################### 

friedman_test<-function(matrix, df, methods, datasets)
{
  i <- 1
  j <- 1
  for (m in methods)
  {
    for (result_prefix in datasets)
    {
      matrix[i, j] <- H_all[H_all$method==m & H_all$result_prefix == result_prefix,]$H[1]
      i <- i+1
    }
    j <- j+1
  }

  friedman.test(matrix)
}

matrix <- matrix(, nrow = 16, ncol = 4)
friedman_test(matrix, H_all, c_methods_ids, c_datasets)
matrix <- matrix(, nrow = 16, ncol = 4)
friedman_test(matrix, H_100, c_methods_ids, c_datasets)
matrix <- matrix(, nrow = 16, ncol = 4)
friedman_test(matrix, H_10, c_methods_ids, c_datasets)
