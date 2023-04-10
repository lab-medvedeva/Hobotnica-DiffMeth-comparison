library(ggplot2)
library(ggpattern)
library(ggpubr)

c_methods_ids_rrbs <- c("methylsig", "dss_no_smoothing", "dss_smoothing","methylkit", "methylkit_overdispersion", "radmeth", "hmm")
c_methods_rrbs <-c("methylSig", "DSS (no smoothing)", "DSS (smoothing)", "methylkit", "methylkit (overdispersion)", "RADMeth", "HMM-DM")
c_methods <-c("methylSig", "DSS (no smoothing)", "DSS (smoothing)", "methylkit", "methylkit (overdispersion)", "RADMeth", "HMM-DM", "BSmooth")

################# alt ################# 
setwd('alt')
TPR_all_diffs_df <- data.frame(dataset=character(),
                                  method=character(),
                                  TPR=double(),
                                  simul = numeric(),
                                  methyl_diff = numeric(),
                                  stringsAsFactors=FALSE)
precision_all_diffs_df <- data.frame(dataset=character(),
                           method=character(),
                           precision=double(),
                           simul = numeric(),
                           methyl_diff = numeric(),
                           stringsAsFactors=FALSE)

H_all_diffs_df <-  data.frame(dataset=character(),
                    method=character(),
                    H=double(),
                    type=character(),
                    pvalue=double(),
                    length=numeric(),
                    simul = numeric(),
                    methyl_diff = numeric(),
                    stringsAsFactors=FALSE)

accuracy_all_diffs_df <- data.frame(dataset=character(),
                                     method=character(),
                                     accuracy=double(),
                                     simul = numeric(),
                                     methyl_diff = numeric(),
                                     stringsAsFactors=FALSE)
methyl_diffs = c(10,15,20,30)
for (methyl_diff in methyl_diffs)
{
  H_df <-  data.frame(dataset=character(),
                      method=character(),
                      H=double(),
                      type=character(),
                      pvalue=double(),
                      length=numeric(),
                      simul = numeric(),
                      stringsAsFactors=FALSE)
  
  TPR_df <- data.frame(dataset=character(),
                       method=character(),
                       TPR=double(),
                       simul = numeric(),
                       stringsAsFactors=FALSE)

  precision_df <- data.frame(dataset=character(),
                       method=character(),
                       precision=double(),
                       simul = numeric(),
                       stringsAsFactors=FALSE)  
  for (simul in 0:9)
  {
    input_folder = paste0('simulation', simul)
    H_full<-read.table(paste0(input_folder, "\\H_full_", methyl_diff), header=TRUE, sep='\t')
    H_full$simul <- simul
    print(dim(H_full))
    H_df <- rbind(H_df,H_full)

    H_full$methyl_diff <- methyl_diff
    H_all_diffs_df <- rbind(H_all_diffs_df,H_full)
    
    TPR_full<-read.table(paste0(input_folder, "\\TPR_full_", methyl_diff), header=TRUE, sep='\t')
    TPR_full$simul <- simul
    print(dim(TPR_full))
    TPR_df <- rbind(TPR_df,TPR_full)
    
    TPR_full$methyl_diff <- methyl_diff
    TPR_all_diffs_df <- rbind(TPR_all_diffs_df,TPR_full)
    
    precision_full<-read.table(paste0(input_folder, "\\precision_full_", methyl_diff), header=TRUE, sep='\t')
    precision_full$simul <- simul
    print(dim(precision_full))
    precision_df <- rbind(precision_df,precision_full)
    
    precision_full$methyl_diff <- methyl_diff
    precision_all_diffs_df <- rbind(precision_all_diffs_df,precision_full)
    
    accuracy_full<-read.table(paste0(input_folder, "\\accuracy_full_", methyl_diff), header=TRUE, sep='\t')
    accuracy_full$simul <- simul
    print(dim(accuracy_full))
    accuracy_df <- rbind(accuracy_df,accuracy_full)
    
    accuracy_full$methyl_diff <- methyl_diff
    accuracy_all_diffs_df <- rbind(accuracy_all_diffs_df,accuracy_full)
  }
  
  H_df$method<-factor(H_df$method, levels = c_methods_ids_rrbs)
  TPR_df$method<-factor(TPR_df$method, levels = c_methods_ids_rrbs)
  TPR_all_diffs_df$method<-factor(TPR_all_diffs_df$method, levels = c_methods_ids_rrbs)
  precision_df$method<-factor(precision_df$method, levels = c_methods_ids_rrbs)
  accuracy_df$method<-factor(accuracy_df$method, levels = c_methods_ids_rrbs)
  
  p <- ggplot(H_df,aes(x=method, y = H, fill=method))+ 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom",
          text = element_text(size=26, face = "bold"),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 22),
          legend.key.size = unit(0.5, "cm"),
          legend.margin = margin(t = -0.3, unit='cm'))+
    geom_boxplot()+
    labs(x ='Method', y = "H") + 
    scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
    scale_x_discrete(drop=FALSE) +
    ggtitle(paste0("Hobotnica, methylation difference ",methyl_diff,"%"))
  
  ggsave(paste0(paste0("H_full_", methyl_diff),".png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")
  
  H_df$H<-log(H_df$H)
  H_df <- H_df[!is.infinite(H_df$H),]
  
  p <- ggplot(H_df,aes(x=method, y = H, fill=method))+ 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom",
          text = element_text(size=26, face = "bold"),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 22),
          legend.key.size = unit(0.5, "cm"),
          legend.margin = margin(t = -0.3, unit='cm'))+
    geom_boxplot()+
    labs(x ='Method', y = "log(H)") + 
    scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
    scale_x_discrete(drop=FALSE) +
    ggtitle(paste0("Hobotnica, methylation difference ",methyl_diff,"%"))
  
  ggsave(paste0(paste0("H_full_log_", methyl_diff),".png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")
  
  H_df$length<-log(H_df$length)
  H_df <- H_df[!is.infinite(H_df$length),]
  
  p <- ggplot(H_df,aes(x=method, y = length, fill=method))+ 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom",
          text = element_text(size=26, face = "bold"),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 22),
          legend.key.size = unit(0.5, "cm"),
          legend.margin = margin(t = -0.3, unit='cm'))+
    geom_boxplot()+
    labs(x ='Method', y = "log(length)") + 
    scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
    scale_x_discrete(drop=FALSE) +
    ggtitle(paste0("Signature length, methylation difference ",methyl_diff,"%"))
  
  ggsave(paste0(paste0("H_full_length_", methyl_diff),".png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")
  
  p <- ggplot(TPR_df,aes(x=method, y = TPR, fill=method))+ 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom",
          text = element_text(size=26, face = "bold"),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 22),
          legend.key.size = unit(0.5, "cm"),
          legend.margin = margin(t = -0.3, unit='cm'))+
    geom_boxplot()+
    labs(x ='Method', y = "Recall") + 
    scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
    scale_x_discrete(drop=FALSE)+
    ggtitle(paste0("Recall, methylation difference ",methyl_diff,"%"))  +
    ylim(0, 1)
  
  ggsave(paste0(paste0("TPR_full_", methyl_diff),".boxplot.png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")
  
  
  p <- ggplot(precision_df,aes(x=method, y = precision, fill=method))+ 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom",
          text = element_text(size=26, face = "bold"),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 22),
          legend.key.size = unit(0.5, "cm"),
          legend.margin = margin(t = -0.3, unit='cm'))+
    geom_boxplot()+
    labs(x ='Method', y = "Precision") + 
    scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
    scale_x_discrete(drop=FALSE)+
    ggtitle(paste0("Precision, methylation difference ",methyl_diff,"%"))  +
    ylim(0, 1)
  
  ggsave(paste0(paste0("precision_full_", methyl_diff),".boxplot.png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")

  accuracy_df$accuracy<-log(accuracy_df$accuracy)
  accuracy_df <- accuracy_df[!is.infinite(accuracy_df$accuracy),]
  
  p <- ggplot(accuracy_df,aes(x=method, y = accuracy, fill=method))+ 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom",
          text = element_text(size=26, face = "bold"),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 22),
          legend.key.size = unit(0.5, "cm"),
          legend.margin = margin(t = -0.3, unit='cm'))+
    geom_boxplot()+
    labs(x ='Method', y = "log(Accuracy)") + 
    scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
    scale_x_discrete(drop=FALSE)+
    ggtitle(paste0("Accuracy, methylation difference ",methyl_diff,"%"))
  
  ggsave(paste0(paste0("accuracy_full_", methyl_diff),".boxplot.png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")

}

TPR_all_diffs_df$methyl_diff <- factor(TPR_all_diffs_df$methyl_diff, levels = c(10,15,20,30))
precision_all_diffs_df$methyl_diff <- factor(precision_all_diffs_df$methyl_diff, levels = c(10,15,20,30))
H_all_diffs_df$methyl_diff <- factor(H_all_diffs_df$methyl_diff, levels = c(10,15,20,30))
accuracy_all_diffs_df$methyl_diff <- factor(accuracy_all_diffs_df$methyl_diff, levels = c(10,15,20,30))

TPR_all_diffs_df$method<-factor(TPR_all_diffs_df$method, levels = c_methods_ids_rrbs)
precision_all_diffs_df$method<-factor(precision_all_diffs_df$method, levels = c_methods_ids_rrbs)
H_all_diffs_df$method<-factor(H_all_diffs_df$method, levels = c_methods_ids_rrbs)
accuracy_all_diffs_df$method<-factor(accuracy_all_diffs_df$method, levels = c_methods_ids_rrbs)

p_rec <- ggplot(TPR_all_diffs_df,aes(x=method, y = TPR, fill=method, alpha = methyl_diff ))+ 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none",
        text = element_text(size=28), #, face = "bold"
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(t = -0.3, unit='cm'),
        plot.background = element_rect(color = "black", size = 2))+
  geom_boxplot()+
  labs(x ='Method', y = "Recall") + 
  scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
  scale_x_discrete(drop=FALSE)+
  ggtitle(paste0("Recall"))  +
  ylim(0, 1) +
  scale_alpha_manual(name = "Methylation difference",values=c(0.1,0.4,0.7,1)) +
  guides(alpha=guide_legend(ncol=2, nrow = 4, override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(0.1,0.4,0.7,1)),
                                                                colour=NA)),
         fill = guide_legend(ncol=2))

ggsave("TPR_full_boxplot.png",  plot = p_rec,  device = "png", width = 30,  height = 20,  units = "cm")

p_prec <- ggplot(precision_all_diffs_df,aes(x=method, y = precision, fill=method, alpha = methyl_diff ))+ 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none",
        text = element_text(size=28), # , face = "bold"
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(t = -0.3, unit='cm'),
        plot.background = element_rect(color = "black", size = 2))+
  geom_boxplot()+
  labs(x ='Method', y = "Precision") + 
  scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
  scale_x_discrete(drop=FALSE)+
  ggtitle(paste0("Precision"))  +
  ylim(0, 1) +
  scale_alpha_manual(name = "Methylation difference",values=c(0.1,0.4,0.7,1)) +
  guides(alpha=guide_legend(ncol=2, nrow = 4, override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(0.1,0.4,0.7,1)),
                                                                colour=NA)),
         fill = guide_legend(ncol=2))

ggsave("precision_full_boxplot.png",  plot = p_prec,  device = "png", width = 30,  height = 20,  units = "cm")

p_h <- ggplot(H_all_diffs_df,aes(x=method, y = H, fill=method, alpha = methyl_diff ))+ 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none",
        text = element_text(size=28), #, face = "bold"
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(t = -0.3, unit='cm'),
        plot.background = element_rect(color = "black", size = 2))+
  geom_boxplot()+
  labs(x ='Method', y = "H") + 
  scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
  scale_x_discrete(drop=FALSE)+
  ggtitle(paste0("Hobotnica"))  +
  ylim(0, 1) +
  scale_alpha_manual(name = "Methylation difference",values=c(0.1,0.4,0.7,1)) +
  guides(alpha=guide_legend(ncol=2, nrow = 4, override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(0.1,0.4,0.7,1)),
                                                                colour=NA)),
         fill = guide_legend(ncol=2))

ggsave("H_full_boxplot.png",  plot = p_h,  device = "png", width = 30,  height = 20,  units = "cm")

H_all_diffs_df$length<-log(H_all_diffs_df$length)
H_all_diffs_df <- H_all_diffs_df[!is.infinite(H_all_diffs_df$length),]

p <- ggplot(H_all_diffs_df,aes(x=method, y = length, fill=method, alpha = methyl_diff ))+ 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none",
        text = element_text(size=28), #, face = "bold"
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(t = -0.3, unit='cm'),
        plot.background = element_rect(color = "black", size = 2))+
  geom_boxplot()+
  labs(x ='Method', y = "log(length)") + 
  scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
  scale_x_discrete(drop=FALSE)+
  scale_alpha_manual(name = "Methylation difference",values=c(0.1,0.4,0.7,1)) +
  guides(alpha=guide_legend(ncol=2, nrow = 4, override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(0.1,0.4,0.7,1)),
                                                                colour=NA)),
         fill = guide_legend(ncol=2))

ggsave("length_boxplot.png",  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")

accuracy_all_diffs_df$accuracy<-log(accuracy_all_diffs_df$accuracy)
accuracy_all_diffs_df <- accuracy_all_diffs_df[!is.infinite(accuracy_all_diffs_df$accuracy),]

p_ac <- ggplot(accuracy_all_diffs_df,aes(x=method, y = accuracy, fill=method, alpha = methyl_diff ))+ 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none",
        text = element_text(size=28), # , face = "bold"
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(t = -0.3, unit='cm'),
        plot.background = element_rect(color = "black", size = 2))+
  geom_boxplot()+
  labs(x ='Method', y = "log(Accuracy)") + 
  scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
  scale_x_discrete(drop=FALSE)+
  ggtitle(paste0("Accuracy"))  +
  scale_alpha_manual(name = "Methylation difference",values=c(0.1,0.4,0.7,1)) +
  guides(alpha=guide_legend(ncol=2, nrow = 4, override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(0.1,0.4,0.7,1)),
                                                                colour=NA)),
         fill = guide_legend(ncol=2))

ggsave("accuracy_full_boxplot.png",  plot = p_ac,  device = "png", width = 30,  height = 20,  units = "cm")


p1 <- ggarrange(p_rec, ggplot() + theme_void(), p_prec, nrow=1,
                widths = c(1, 0.05, 1),
                legend='none')

ggsave("row1.png",  plot = p1,  device = "png", width = 60,  height = 20,  units = "cm")

############################### perm  ############################### 
setwd('../perm')
FPR_all_diffs_df <- data.frame(dataset=character(),
                     method=character(),
                     FPR=double(),
                     simul = numeric(),
                     methyl_diff = numeric(),
                     stringsAsFactors=FALSE)

H_all_diffs_df <-  data.frame(dataset=character(),
                    method=character(),
                    H=double(),
                    type=character(),
                    pvalue=double(),
                    length=numeric(),
                    simul = numeric(),
                    methyl_diff = numeric(),
                    stringsAsFactors=FALSE)

methyl_diffs = c(10,15,20,30)
for (methyl_diff in methyl_diffs)
{
  H_df <-  data.frame(dataset=character(),
                      method=character(),
                      H=double(),
                      type=character(),
                      pvalue=double(),
                      length=numeric(),
                      simul = numeric(),
                      stringsAsFactors=FALSE)
  
  FPR_df <- data.frame(dataset=character(),
                       method=character(),
                       FPR=double(),
                       simul = numeric(),
                       stringsAsFactors=FALSE)
  
  for (simul in 0:9)
  {
    input_folder = paste0('simulation', simul)
    H_full<-read.table(paste0(input_folder, "\\H_full_", methyl_diff), header=TRUE, sep='\t')
    H_full$simul <- simul
    print(dim(H_full))
    H_df <- rbind(H_df,H_full)    
    H_full$methyl_diff <- methyl_diff
    H_all_diffs_df <- rbind(H_all_diffs_df,H_full)
    FPR_full<-read.table(paste0(input_folder, "\\FPR_full_", methyl_diff), header=TRUE, sep='\t')
    FPR_full$simul <- simul
    print(dim(FPR_full))
    FPR_df <- rbind(FPR_df,FPR_full)
    FPR_full$methyl_diff <- methyl_diff
    FPR_all_diffs_df <- rbind(FPR_all_diffs_df,FPR_full)
    
  }
  
  H_df$method<-factor(H_df$method, levels = c_methods_ids_rrbs)
  FPR_df$method<-factor(FPR_df$method, levels = c_methods_ids_rrbs)    
  
  p <- ggplot(H_df,aes(x=method, y = H, fill=method))+ 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom",
          text = element_text(size=28), # , face = "bold"
          legend.text = element_text(size = 28),
          legend.title = element_text(size = 28),
          legend.key.size = unit(0.5, "cm"),
          legend.margin = margin(t = -0.3, unit='cm'))+
    geom_boxplot()+
    labs(x ='Method', y = "H") + 
    scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
    scale_x_discrete(drop=FALSE) +
    ggtitle(paste0("Hobotnica, methylation difference ",methyl_diff,"%"))
  
  ggsave(paste0(paste0("H_false_full_", methyl_diff),".png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")
  
  H_df <- H_df[H_df$H!=0,]
  
  p <- ggplot(H_df,aes(x=method, y = H, fill=method))+ 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom",
          text = element_text(size=28), #, face = "bold"
          legend.text = element_text(size = 28),
          legend.title = element_text(size = 28),
          legend.key.size = unit(0.5, "cm"),
          legend.margin = margin(t = -0.3, unit='cm'))+
    geom_boxplot()+
    labs(x ='Method', y = "H") + 
    scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
    scale_x_discrete(drop=FALSE) +
    ggtitle(paste0("Hobotnica, methylation difference ", methyl_diff, "%"))
  
  ggsave(paste0(paste0("H_false_full_no_zero_", methyl_diff),".png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")
  
  H_df$H<-log(H_df$H)
  H_df <- H_df[!is.infinite(H_df$H),]
  
  p <- ggplot(H_df,aes(x=method, y = H, fill=method))+ 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom",
          text = element_text(size=28), #, face = "bold"
          legend.text = element_text(size = 28),
          legend.title = element_text(size = 28),
          legend.key.size = unit(0.5, "cm"),
          legend.margin = margin(t = -0.3, unit='cm'))+
    geom_boxplot()+
    labs(x ='Method', y = "log(H)") + 
    scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
    scale_x_discrete(drop=FALSE) +
    ggtitle(paste0("Hobotnica, methylation difference ",methyl_diff,"%"))
  
  ggsave(paste0(paste0("H_false_full_log_", methyl_diff),".png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")
  
  p <- ggplot(H_df,aes(x=method, y = length, fill=method))+ 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom",
          text = element_text(size=28), #, face = "bold"
          legend.text = element_text(size = 28),
          legend.title = element_text(size = 28),
          legend.key.size = unit(0.5, "cm"),
          legend.margin = margin(t = -0.3, unit='cm'))+
    geom_boxplot()+
    labs(x ='Method', y = "length") + 
    scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
    scale_x_discrete(drop=FALSE) +
    ggtitle(paste0("Signature length, methylation difference ",methyl_diff,"%"))
  
  ggsave(paste0(paste0("H_false_full_length_", methyl_diff),".png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")
  
  p <- ggplot(FPR_df,aes(x=method, y = FPR, fill=method))+ 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom",
          text = element_text(size=28), #, face = "bold"
          legend.text = element_text(size = 28),
          legend.title = element_text(size = 28),
          legend.key.size = unit(0.5, "cm"),
          legend.margin = margin(t = -0.3, unit='cm'))+
    geom_boxplot()+
    labs(x ='Method', y = "FPR") + 
    scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
    scale_x_discrete(drop=FALSE)+
    ggtitle(paste0("False positive rate, methylation difference ",methyl_diff,"%"))  +
    ylim(0, 1)
  
  ggsave(paste0(paste0("FPR_full_", methyl_diff),".boxplot.png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")
  
  FPR_df$FPR<-log(FPR_df$FPR)
  FPR_df <- FPR_df[!is.infinite(FPR_df$FPR),]
  
  p <- ggplot(FPR_df,aes(x=method, y = FPR, fill=method))+ 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom",
          text = element_text(size=28), # , face = "bold"
          legend.text = element_text(size = 28),
          legend.title = element_text(size = 28),
          legend.key.size = unit(0.5, "cm"),
          legend.margin = margin(t = -0.3, unit='cm'))+
    geom_boxplot()+
    labs(x ='Method', y = "log(FPR)") + 
    scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
    scale_x_discrete(drop=FALSE)+
    ggtitle(paste0("False positive rate, methylation difference ",methyl_diff,"%"))
  
  ggsave(paste0(paste0("FPR_full_log_", methyl_diff),".boxplot.png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")
}

FPR_all_diffs_df$method<-factor(FPR_all_diffs_df$method, levels = c_methods_ids_rrbs)
FPR_all_diffs_df$methyl_diff <- factor(FPR_all_diffs_df$methyl_diff, levels = c(10,15,20,30))

H_all_diffs_df$method<-factor(H_all_diffs_df$method, levels = c_methods_ids_rrbs)
H_all_diffs_df$methyl_diff <- factor(H_all_diffs_df$methyl_diff, levels = c(10,15,20,30)) 

FPR_all_diffs_df$FPR<-log(FPR_all_diffs_df$FPR)
FPR_all_diffs_df <- FPR_all_diffs_df[!is.infinite(FPR_all_diffs_df$FPR),]


p_fpr <- ggplot(FPR_all_diffs_df,aes(x=method, y = FPR, fill=method, alpha = methyl_diff ))+ 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none",
        text = element_text(size=28), # , face = "bold"
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(t = -0.3, unit='cm'),
        plot.background = element_rect(color = "black", size = 2))+
  geom_boxplot()+
  labs(x ='Method', y = "log(FPR)") + 
  scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
  scale_x_discrete(drop=FALSE)+
  ggtitle(paste0("False positive rate"))  +
  scale_alpha_manual(name = "Methylation difference",values=c(0.1,0.4,0.7,1)) +
  guides(alpha=guide_legend(ncol=2, nrow = 4, override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(0.1,0.4,0.7,1)),
                                                                colour=NA)),
         fill = guide_legend(ncol=2))

ggsave("FPR_full_log_boxplot.png",  plot = p_fpr,  device = "png", width = 30,  height = 20,  units = "cm")

p_h_perm_default <- ggplot(H_all_diffs_df,aes(x=method, y = H, fill=method, alpha = methyl_diff ))+ 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none",
        text = element_text(size=28), # , face = "bold"
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(t = -0.3, unit='cm'),
        plot.background = element_rect(color = "black", size = 2))+
  geom_boxplot()+
  labs(x ='Method', y = "H") + 
  scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
  scale_x_discrete(drop=FALSE)+ #     ggtitle(paste0("Hobotnica (permutation)"))  +
  ylim(0, 1) +
  scale_alpha_manual(name = "Methylation difference",values=c(0.1,0.4,0.7,1)) +
  guides(alpha=guide_legend(ncol=2, nrow = 4, override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(0.1,0.4,0.7,1)),
                                                                colour=NA)),
         fill = guide_legend(ncol=2))

ggsave("H_false_full_boxplot.png",  plot = p_h_perm_default,  device = "png", width = 30,  height = 20,  units = "cm")

p <- ggplot(H_all_diffs_df,aes(x=method, y = length, fill=method, alpha = methyl_diff ))+ 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none",
        text = element_text(size=28), # , face = "bold"
        legend.text = element_text(size = 28),
        legend.title = element_text(size = 28),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(t = -0.3, unit='cm'),
        plot.background = element_rect(color = "black", size = 2))+
  geom_boxplot()+
  labs(x ='Method', y = "length") + 
  scale_fill_brewer(name = "Method", labels = c_methods_rrbs, palette="Dark2", drop = F)+ 
  scale_x_discrete(drop=FALSE)+
  ggtitle(paste0("length (FPR)"))  +
  scale_alpha_manual(name = "Methylation difference",values=c(0.1,0.4,0.7,1)) +
  guides(alpha=guide_legend(ncol=2, nrow = 4, override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(0.1,0.4,0.7,1)),
                                                                colour=NA)),
         fill = guide_legend(ncol=2))
ggsave(paste0(paste0("H_false_full_nocut_length_boxplot"),".png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")

# draw precision again
setwd('../alt')
methyl_diff <- 10
prec_df <-  data.frame(precs=double(),
                       Hs=double(),
                       simul = character(),
                       stringsAsFactors=FALSE)
  
for (simul in 0:9)
{
  
  input_folder = paste0('simulation', simul)
  prec_full<-read.table(paste0(input_folder, "\\prec_H_", methyl_diff, "_", simul), header=TRUE, sep='\t')
  prec_full$simul <- paste0('S',simul+1)
  print(dim(prec_full))
  prec_df <- rbind(prec_df,prec_full)
}

prec_df$simul <- factor(prec_df$simul, levels = paste0('S', 1:10))

p_prec10<-ggplot(data=prec_df, aes(x=precs, y=Hs, group = simul, colour = simul)) +
  geom_line() +  labs(x ='precision', y = "H", colour = 'Simulation')+
  scale_fill_brewer(name = "Simulation", palette="Dark2", drop = F)+
  theme_bw() +
  theme(text = element_text(size=28), #, face = "bold"
    panel.border = element_blank(),
    legend.position = c(0.87, 0.25),
    legend.text = element_text(size = 28),
    legend.title = element_text(size = 28),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(color = "black", size = 2))+ # panel.border = element_rect(colour = "black", fill=NA)
  ylim(0.4, 1) +
  guides(colour=guide_legend(ncol=2)) +
  ggtitle("")

ggsave(paste0("prec_", methyl_diff, ".png"),  plot = p_prec10,  device = "png", width = 30,  height = 20,  units = "cm")

p2 <- ggarrange(p_ac,  ggplot() + theme_void(), p_h, nrow=1,  
                widths = c(1, 0.05, 1),
                legend = 'none')
ggsave("row2.png",  plot = p2, device = "png", width = 60,  height = 20,  units = "cm")

p3 <- ggarrange(p_fpr,  ggplot() + theme_void(), p_prec10, nrow=1,
                widths = c(1, 0.05, 1),
                legend = 'none')
ggsave("row3.png",  plot = p3, device = "png", width = 60,  height = 20,  units = "cm")

############################# Draw precision plots ############################# 
  
setwd('alt')
methyl_diffs = c(10,15,20,30)
for (methyl_diff in methyl_diffs)
{
  prec_df <-  data.frame(precs=double(),
                      Hs=double(),
                      simul = character(),
                      stringsAsFactors=FALSE)
  
  for (simul in 0:9)
  {
 
      input_folder = paste0('simulation', simul)
      prec_full<-read.table(paste0(input_folder, "\\prec_H_", methyl_diff, "_", simul), header=TRUE, sep='\t')
      prec_full$simul <- paste0('S',simul+1)
      print(dim(prec_full))
      prec_df <- rbind(prec_df,prec_full)
  }
  
  prec_df$simul <- factor(prec_df$simul, levels = paste0('S', 1:10))
  
  p<-ggplot(data=prec_df, aes(x=precs, y=Hs, group = simul, colour = simul)) +
    geom_line() +  labs(x ='precision', y = "H")+
    scale_fill_brewer(name = "Simulation", palette="Dark2", drop = F)+
    labs(x ='precision', y = "H", colour = 'Simulation')+
    theme_bw() +
    theme(
      text = element_text(size=24), #, face = "bold"
      legend.text = element_text(size = 24),
      legend.title = element_text(size = 24),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank())+
    ggtitle(paste0("Methylation difference = ", methyl_diff, "%"))  +
    ylim(0.4, 1) 
  ggsave(paste0("prec_H_", methyl_diff, ".png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")
}   

###################### Draw heatmaps #######################

dist_matrix_df<-read.csv("GSE103886_15_radmeth_full_15_distmatrix_rank", sep='\t')
colnames(dist_matrix_df) <- c('x','x_type','y','y_type','distance')
ord <- paste0("X", rep(1:23))
dist_matrix_df$x<-factor(dist_matrix_df$x, levels = ord)
dist_matrix_df$y<-factor(dist_matrix_df$y, levels = ord)
dist_matrix_df$x_type<-factor(dist_matrix_df$x_type, levels = c('control','case'))
dist_matrix_df$y_type<-factor(dist_matrix_df$y_type, levels = c('case','control'))
p_rank <- ggplot(dist_matrix_df, aes(x, y)) + 
  facet_grid(y_type~x_type, scales = "free",space="free") +
  geom_tile(aes(fill= distance)) +
  scale_fill_distiller(palette = "RdBu", direction = 1)+  #RdYlBu
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text = element_text(size= 24), #, face = "bold"
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.key.size = unit(1, "cm"),
        legend.key.width =  unit(3, "cm"),
        legend.position = "bottom",
        legend.key=element_rect(colour="black"),
        legend.margin = margin(t = -0.3, unit='cm') )

ggsave(paste0(dataset, "_rank_heatmap.png"),  plot = p_rank,  device = "png", width = 15,  height = 15,  units = "cm")

dist_matrix_df<-read.csv("GSE103886_15_radmeth_full_15_distmatrix", sep='\t')
colnames(dist_matrix_df) <- c('x','x_type','y','y_type','distance')
ord <- paste0("X", rep(1:23))
dist_matrix_df$x<-factor(dist_matrix_df$x, levels = ord)
dist_matrix_df$y<-factor(dist_matrix_df$y, levels = ord)
dist_matrix_df$x_type<-factor(dist_matrix_df$x_type, levels = c('control','case'))
dist_matrix_df$y_type<-factor(dist_matrix_df$y_type, levels = c('case','control'))

p_dist <- ggplot(dist_matrix_df, aes(x, y)) + 
  facet_grid(y_type~x_type, scales = "free",space="free") +
  geom_tile(aes(fill= distance)) +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        text = element_text(size= 24), #, face = "bold"
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24),
        legend.key.size = unit(1, "cm"),
        legend.key.width =  unit(1, "cm"),
        legend.position = "bottom",
        legend.key=element_rect(colour="black"),
        legend.margin = margin(t = -0.3, unit='cm') )

ggsave(paste0(dataset, "_heatmap.png"),  plot = p_dist,  device = "png", width = 15,  height = 15,  units = "cm")
  
  
  
  

  
  