library (ggplot2)
library(ggpattern)
library(ggpubr)
library(venn)

# Order: 'GSE149608','GSE138598','GSE119980','GSE103886', 'GSE150592', 'GSE148060'
c_datasets <- c('GSE149608','GSE138598','GSE119980','GSE103886', 'GSE150592', 'GSE148060')

c_methods <-c("methylSig", "DSS (no smoothing)", "DSS (smoothing)", "methylkit", "methylkit (overdispersion)", "RADMeth", "HMM-DM", "BSmooth")
c_methods_ids <- c("methylsig", "dss_no_smoothing", "dss_smoothing","methylkit", "methylkit_overdispersion", "radmeth", "hmm", "bsmooth")
method_names <- c(
  'dss_no_smoothing' = "DSS (no smoothing)",
  'dss_smoothing' = "DSS (smoothing)",
  'methylsig' = "methylSig",
  'methylkit' = "methylkit",
  'methylkit_overdispersion' = "methylkit (overdispersion)",
  'radmeth' = "RADMeth",
  'hmm' = "HMM-DM",
  'bsmooth' = "BSmooth")

c_methods_ids_wgbs <- c("methylsig", "dss_no_smoothing", "dss_smoothing","methylkit", "methylkit_overdispersion", "radmeth", "hmm", "bsmooth")
c_methods_ids_rrbs <- c("methylsig", "dss_no_smoothing", "dss_smoothing","methylkit", "methylkit_overdispersion", "radmeth", "hmm")

c_methods_wgbs <-c("methylSig", "DSS (no smoothing)", "DSS (smoothing)", "methylkit", "methylkit (overdispersion)", "RADMeth", "HMM-DM", "BSmooth")
c_methods_rrbs <-c("methylSig", "DSS (no smoothing)", "DSS (smoothing)", "methylkit", "methylkit (overdispersion)", "RADMeth", "HMM-DM")

c_datasets_wgbs <- c('GSE149608', 'GSE138598', 'GSE119980')
c_datasets_rrbs <- c('GSE103886', 'GSE150592', 'GSE148060')

dark2 <- RColorBrewer::brewer.pal(8,'Dark2')
fixed_dark2 <-  c( dark2[1:7])

################# read H-scores #######################

H_all_100 <- read.table("H_all_100",sep='\t', header = T)
H_all_100$tag<-factor(H_all_100$tag, levels = c('H_all', 'H_100', 'H_10'))
H_all_100$length <- as.numeric(H_all_100$length)
H_all_100$dataset<-factor(H_all_100$dataset,levels = c_datasets)
H_all_100$method<-factor(H_all_100$method, levels = c_methods_ids)

tag_names <- c(
  'H_all' = "Full",
  'H_100' = "Top 100 DMC",
  'H_10' = "Top 10 DMC")

c_tag_names <-  c("Full", "Top 100 DMC", "Top 10 DMC")
H_all_100$type<-factor(H_all_100$type,levels = c('WGBS','RRBS'))
H_all_100_wgbs <- H_all_100[H_all_100$type=='WGBS',]
H_all_100_rrbs <- H_all_100[H_all_100$type=='RRBS',]

H_all_100_rrbs$method<-factor(H_all_100_rrbs$method, levels = c_methods_ids_rrbs)
H_all_100_rrbs$dataset<-factor(H_all_100_rrbs$dataset,levels = c_datasets_rrbs)
H_all_100_wgbs$method<-factor(H_all_100_wgbs$method, levels = c_methods_ids_wgbs)
H_all_100_wgbs$dataset<-factor(H_all_100_wgbs$dataset,levels = c_datasets_wgbs)

####################### Example #######################
dataset <- 'GSE150592'

################ Heatmap

method_name <- 'dss_no_smoothing'
dss_no_smoothing <- read.csv(paste0(dataset, '_dss_no_smoothing_cut'),sep = '\t')
if(dim(dss_no_smoothing)[1] == 0)
{
  dss_no_smoothing_signature <- c()
} else {
  dss_no_smoothing_signature<-paste0(dss_no_smoothing$chr,":",dss_no_smoothing$pos)
}

sig_length<-length(dss_no_smoothing_signature)
ratio_table<-read.table(file = "ratio_table", sep="\t")
heat_matrix<-ratio_table[dss_no_smoothing_signature,]
distMatrix<-dist(t(heat_matrix))
distMatrix_df <- data.frame(as.matrix(distMatrix))
distMatrix_df_flat <- data.frame(x=character(),
                                 x_type=character(),
                                 y=character(),
                                 y_type=character(),
                                 d=double(),
                                 stringsAsFactors=FALSE)

trt<-c('group2','group1','group2','group1','group2','group1','group2','group1','group2','group1','group2','group1','group1','group2','group2','group1','group1','group2','group2','group1','group1','group2','group2','group1','group2','group1','group2','group1','group1','group2')
case <- which(trt == "group2")
control <- which(trt == "group1")
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
    names(new_row) <- c("x","x_type","y","y_type","d")
    distMatrix_df_flat <- rbind(distMatrix_df_flat, new_row)
  }
}

write.table(distMatrix_df_flat, file = paste(dataset, 'dss_no_smoothing', "distmatrix", sep = "_"), sep="\t", col.names=TRUE, quote = FALSE)
dist_matrix_df<-read.csv(paste(dataset, 'dss_no_smoothing', "distmatrix", sep = "_"), 
                         sep='\t')
colnames(dist_matrix_df) <- c('x','x_type','y','y_type','distance')
ord <- c('SRR11790875','SRR11790876','SRR11790877','SRR11790878','SRR11790879','SRR11790880','SRR11790881','SRR11790882','SRR11790883','SRR11790884','SRR11790885','SRR11790886','SRR11790887','SRR11790888','SRR11790889','SRR11790890','SRR11790891','SRR11790892','SRR11790893','SRR11790894','SRR11790895','SRR11790896','SRR11790897','SRR11790898','SRR11790899','SRR11790900','SRR11790901','SRR11790902','SRR11790903','SRR11790904')

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
        text = element_text(size= 15), #, face = "bold"
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.key.size = unit(0.6, "cm"),
        legend.key.width =  unit(0.8, "cm"),
        legend.position = "bottom",
        legend.key=element_rect(colour="black"),
        legend.margin = margin(t = -0.3, unit='cm') )

ggsave(paste0(dataset,"_example_heatmap.png"),  plot = p,  device = "png", width = 15,  height = 15,  units = "cm")

#############  Hist
H_GSE150592 <- H_all_100[H_all_100$dataset == 'GSE150592',]
H_GSE150592$length <- as.numeric(H_GSE150592$length)
H_GSE150592$dataset<-factor(H_GSE150592$dataset,levels = c_datasets)
H_GSE150592$method<-factor(H_GSE150592$method, levels = c_methods_ids)
H_GSE150592$tag<-factor(H_GSE150592$tag, levels = c('H_all', 'H_100', 'H_10'))

H_point<- H_GSE150592[H_GSE150592$method == 'dss_no_smoothing' & H_GSE150592$tag == 'H_all',]$H

H_random <- read.csv('GSE150592_dss_no_smoothing_H_random_445', header = F)
colnames(H_random) <-c('H')
p <- ggplot(H_random, aes(x=H)) + 
  geom_histogram(color="darkblue", fill="lightblue", bins = 300, size = 0.05) +
  geom_point(aes(x=H_point, y = 0), colour = "red", size=2) +
  scale_y_continuous(limits = c(0, 500)) +
  scale_x_continuous(limits = c(0.4, 1), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  labs(x ='H-score', y = "Count") +
  theme_bw() +
  theme(panel.border = element_blank(), # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(size = 15), #, face = "bold"
        axis.title=element_text(size = 15), #,face="bold"
        axis.text = element_text(size = 15)) +
  coord_cartesian(clip = "off")

ggsave(paste0(dataset,"_example_hist.png"),  plot = p,  device = "png", width = 15,  height = 15,  units = "cm")

##############  Barplot
H_all_100_selected <-H_all_100[H_all_100$dataset=='GSE150592',]
H_all_100_selected$method<-factor(H_all_100_selected$method, levels = c_methods_ids_rrbs)
H_all_100_selected$dataset<-factor(H_all_100_selected$dataset,levels = c_datasets_rrbs)
H_all_100_selected$tag<-factor(H_all_100_selected$tag, levels = c('H_all', 'H_100', 'H_10'))
H_all_100_selected$length <- as.numeric(H_all_100_selected$length)

p <- ggplot(H_all_100_selected,aes(x=method, y = H, fill=method, pattern = tag))+ 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        text = element_text(size=15), #, face = "bold"
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
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
         fill = guide_legend(ncol=2, override.aes = list(pattern = "none"))) +
  labs(x ='Method', y = "H", pattern = "Signature type") + 
  scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
  scale_pattern_manual(values = c('none', 'stripe', 'circle'),
                       labels = c_tag_names)+
  scale_x_discrete(drop=FALSE)

ggsave(paste0("GSE150592_example_barplot.png"),  plot = p,  device = "png", width = 30,  height = 10,  units = "cm")

############### datasets  ############### 

p1 <- ggplot(H_all_100_wgbs,aes(x=method, y = H, fill=method, pattern = tag))+
  facet_wrap(~ dataset,  ncol = 1)+
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
  scale_fill_brewer(name = "Method", labels = c_methods, palette="Dark2", drop = F)+ 
  scale_pattern_manual(values = c('none', 'stripe', 'circle'),
                       labels = c_tag_names)+
  scale_x_discrete(drop=FALSE)

p2 <- ggplot(H_all_100_rrbs,aes(x=method, y = H, fill=method, pattern = tag))+
  facet_wrap(~ dataset,  ncol = 1)+
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
  scale_fill_brewer(name = "Method", labels = c_methods, palette="Dark2", drop = F)+ 
  scale_pattern_manual(values = c('none', 'stripe', 'circle'),
                       labels = c_tag_names)+
  scale_x_discrete(drop=FALSE)


p3 <- ggarrange(p1, p2, ncol=2, nrow=1, legend="bottom", common.legend = TRUE,
                legend.grob = get_legend(p1))

ggsave(paste0("datasets_all.png"),  plot = p3,  device = "png", width = 40,  height = 20,  units = "cm")

###################### methods ###################### 
p <- ggplot(H_all_100,aes(x=dataset, y = H, fill=dataset, pattern = tag))+
  facet_wrap(~ method,  ncol = 2,labeller = as_labeller(method_names)) + 
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
  labs(x ='Dataset', y = "H", pattern = "Signature type") + 
  scale_fill_discrete(name = "Dataset", drop = F)+ 
  scale_pattern_manual(values = c('none', 'stripe', 'circle'),
                       labels = c_tag_names)+
  scale_x_discrete(drop=FALSE)


ggsave(paste0("methods_all.png"),  plot = p,  device = "png", width = 40,  height = 20,  units = "cm")

######################## H-score boxplots ######################## 

tag_names_reduced <- c(
  'H_all' = "Full",
  'H_100' = "Top 100 DMC")

H_all_100_wgbs_with_zero <- H_all_100_wgbs
H_all_100_wgbs_with_zero$is_H_zero <- TRUE
H_all_100_wgbs_without_zero <- H_all_100_wgbs
H_all_100_wgbs_without_zero <- H_all_100_wgbs_viol_without_zero[H_all_100_wgbs_without_zero$H!=0,]
H_all_100_wgbs_without_zero$is_H_zero <- FALSE

data<- rbind(H_all_100_wgbs_with_zero,H_all_100_wgbs_without_zero )
data <- data[data$tag != 'H_10',]
data$is_H_zero <- factor(data$is_H_zero, levels = c(TRUE, FALSE))
data$tag <- factor(data$tag, levels = c('H_all', 'H_100'))
data$method <-factor(data$method, levels = c_methods_ids_wgbs)

p1<-ggplot(data, aes(x=method, y=H, fill=method, alpha = is_H_zero )) + 
  facet_wrap(~ tag,  ncol = 2, labeller = as_labeller(tag_names_reduced)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        legend.position="bottom",
        text = element_text(size=20),# , face = "bold"
        title = element_text(size=20)) + # , face = "bold"
  geom_boxplot(fatten = 0.8) +  labs(x ='Method', y = "H") + 
  scale_y_continuous(limits = c(0, 1))+
  scale_x_discrete( drop=FALSE) + 
  scale_fill_brewer(name = "Method", labels = c_methods,palette="Dark2",drop = F) +
  scale_alpha_manual(name = "Include H=0",values=c(1,0.3)) +
  guides(alpha=guide_legend(ncol=1, nrow = 2, override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(1,0.3)),
                                                                colour=NA)),
         fill = guide_legend(ncol=2))+
  ggtitle("H-scores for WGBS datasets")

H_all_100_rrbs_with_zero <- H_all_100_rrbs
H_all_100_rrbs_with_zero$is_H_zero <- TRUE
H_all_100_rrbs_without_zero <- H_all_100_rrbs
H_all_100_rrbs_without_zero <- H_all_100_rrbs_without_zero[H_all_100_rrbs_without_zero$H!=0,]
H_all_100_rrbs_without_zero$is_H_zero <- FALSE

data<- rbind(H_all_100_rrbs_with_zero, H_all_100_rrbs_without_zero)
data <- data[data$tag != 'H_10',]
data$is_H_zero <- factor(data$is_H_zero, levels = c(TRUE, FALSE))
data$tag <- factor(data$tag, levels = c('H_all', 'H_100', 'H_10'))
data$method <-factor(data$method, levels = c_methods_ids_rrbs)

p2 <- ggplot(data, aes(x=method, y=H, fill=method, alpha = is_H_zero )) + 
  facet_wrap(~ tag,  ncol = 2, labeller = as_labeller(tag_names_reduced)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        legend.position="bottom",
        text = element_text(size=20), #, face = "bold"
        title = element_text(size=20)) + #, face = "bold"
  geom_boxplot(fatten = 0.8) +  labs(x ='Method', y = "H") + 
  scale_y_continuous(limits = c(0, 1))+
  scale_x_discrete( drop=FALSE) + 
  scale_fill_brewer(name = "Method", labels = c_methods,palette="Dark2",drop = F) +
  scale_alpha_manual(name = "Include H=0",values=c(1,0.3)) +
  guides(alpha=guide_legend(ncol=1, nrow = 2, override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(1,0.3)), colour=NA)),
         fill = guide_legend(ncol=2))+
  ggtitle("H-scores for RRBS datasets")

p <- ggarrange(p1, p2, ncol=1, nrow=2, common.legend = TRUE, legend="bottom")
ggsave(paste0("methods_boxplot_ngs.png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")


################# Shapiro test ################# 
H_all <- H_all_100[H_all_100$tag == 'H_all',]
H_100 <- H_all_100[H_all_100$tag == 'H_100',]
H_10 <- H_all_100[H_all_100$tag == 'H_10',]
H_all_wgbs <-  H_all[H_all$type=='WGBS',]
H_100_wgbs <-  H_100[H_100$type=='WGBS',]
H_10_wgbs <-  H_10[H_10$type=='WGBS',]

H_all_rrbs <-  H_all[H_all$type=='RRBS',]
H_100_rrbs <-  H_100[H_100$type=='RRBS',]
H_10_rrbs <-  H_10[H_10$type=='RRBS',]

shapiro.test(H_all_wgbs$H)
shapiro.test(H_100_wgbs$H)
shapiro.test(H_10_wgbs$H)

shapiro.test(H_all_rrbs$H)
shapiro.test(H_100_rrbs$H)
shapiro.test(H_10_rrbs$H)

H_all_nb <- H_all[H_all$method != 'bsmooth',]
H_100_nb <-H_100[H_100$method != 'bsmooth',]
H_10_nb <- H_10[H_10$method != 'bsmooth',]

shapiro.test(H_all_nb$H)
shapiro.test(H_100_nb$H)
shapiro.test(H_10_nb$H)

################# Friedman test ################# 

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

#WGBS, all
wgbs_matrix <- matrix(, nrow = 3, ncol = 8)
friedman_test(wgbs_matrix, H_all_wgbs, c_methods_ids_wgbs, c_datasets_wgbs)

#WGBS, 100
wgbs_matrix <- matrix(, nrow = 3, ncol = 8)\
friedman_test(wgbs_matrix, H_100_wgbs, c_methods_ids_wgbs, c_datasets_wgbs)

#WGBS, 10
wgbs_matrix <- matrix(, nrow = 3, ncol = 8)\
friedman_test(wgbs_matrix, H_10_wgbs, c_methods_ids_wgbs, c_datasets_wgbs)

#RRBS, all
rrbs_matrix <- matrix(, nrow = 3, ncol = 7)
friedman_test(rrbs_matrix, H_all_rrbs, c_methods_ids_rrbs, c_datasets_rrbs)

#RRBS, 100
rrbs_matrix <- matrix(, nrow = 3, ncol = 7)
friedman_test(rrbs_matrix, H_100_rrbs, c_methods_ids_rrbs, c_datasets_rrbs)

#RRBS, 10
rrbs_matrix <- matrix(, nrow = 3, ncol = 7)
friedman_test(rrbs_matrix, H_10_rrbs, c_methods_ids_rrbs, c_datasets_rrbs)

# all types, all
matrix <- matrix(, nrow = 6, ncol = 7)
friedman_test(matrix, H_all_nb, c_methods_ids_rrbs, c_datasets)

# all types, 100
matrix <- matrix(, nrow = 6, ncol = 7)
friedman_test(matrix, H_100_nb, c_methods_ids_rrbs, c_datasets)

# all types, 10
matrix <- matrix(, nrow = 6, ncol = 7)
friedman_test(matrix, H_10_nb, c_methods_ids_rrbs, c_datasets)
