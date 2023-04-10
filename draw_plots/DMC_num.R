library (ggplot2)
library(ggpubr)

########## NGS
# Order: 'GSE149608','GSE138598','GSE119980','GSE103886', 'GSE150592', 'GSE148060'
c_datasets <- c('GSE149608','GSE138598','GSE119980','GSE103886', 'GSE150592', 'GSE148060')

c_methods <-c("BSmooth","methylSig", "DSS (no smoothing)", "DSS (smoothing)", "methylkit", "methylkit (overdispersion)", "HMM-DM", "RADMeth")
c_methods_ids <- c("bsmooth", "methylsig", "dss_no_smoothing", "dss_smoothing","methylkit", "methylkit_overdispersion","hmm", "radmeth")
c_methods <-c("BSmooth","methylSig", "DSS (no smoothing)", "DSS (smoothing)", "methylkit", "methylkit (overdispersion)", "HMM-DM", "RADMeth")
method_names <- c(
  'bsmooth' = "BSmooth",
  'dss_no_smoothing' = "DSS (no smoothing)",
  'dss_smoothing' = "DSS (smoothing)",
  'methylsig' = "methylSig",
  'methylkit' = "methylkit",
  'methylkit_overdispersion' = "methylkit (overdispersion)",
  'hmm' = "HMM-DM",
  'radmeth' = "RADMeth")

# Methods evaluation

c_methods_ids_wgbs <- c("bsmooth", "methylsig", "dss_no_smoothing", "dss_smoothing","methylkit", "methylkit_overdispersion","hmm", "radmeth")
c_methods_ids_rrbs <- c("methylsig", "dss_no_smoothing", "dss_smoothing","methylkit", "methylkit_overdispersion","hmm", "radmeth")

c_methods_wgbs <-c("BSmooth","methylSig", "DSS (no smoothing)", "DSS (smoothing)", "methylkit", "methylkit (overdispersion)", "HMM-DM", "RADMeth")
c_methods_rrbs <-c("methylSig", "DSS (no smoothing)", "DSS (smoothing)", "methylkit", "methylkit (overdispersion)", "HMM-DM", "RADMeth")

c_datasets_wgbs <- c('GSE149608', 'GSE138598', 'GSE119980')
c_datasets_rrbs <- c('GSE103886', 'GSE150592', 'GSE148060')

## Number of DMCs
H_all <- read.table("H_all_ngs.csv",sep='\t', header = T)
H_all$length <- as.numeric(H_all$length)
H_all$dataset<-factor(H_all$dataset,levels = c_datasets)
H_all$method<-factor(H_all$method, levels = c_methods_ids)

H_all_wgbs <- H_all[H_all$type=='WGBS',]
H_all_rrbs <- H_all[H_all$type=='RRBS',]

H_all_wgbs$method<-factor(H_all_wgbs$method, levels = c_methods_ids_wgbs)
H_all_rrbs$method<-factor(H_all_rrbs$method, levels = c_methods_ids_rrbs)

H_all_wgbs$dataset<-factor(H_all_wgbs$dataset,levels = c_datasets_wgbs)
H_all_rrbs$dataset<-factor(H_all_rrbs$dataset,levels = c_datasets_rrbs)

H_all_wgbs$length <- as.numeric(H_all_wgbs$length)
H_all_rrbs$length <- as.numeric(H_all_rrbs$length)

dark2 <- RColorBrewer::brewer.pal(8,'Dark2')
fixed_dark2 <-  c( dark2[2:length(dark2)])

p21 <- ggplot(H_all_wgbs, aes(x=method, y = length, fill=method))  +
  facet_wrap(~ dataset, scales = "free",  ncol = 1) +
  geom_bar(stat="identity", position=position_dodge())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), legend.position="bottom", 
        text = element_text(size=20, face = "bold")) +
  labs(x ='Method', y = "log(Number of DMCs)") + 
  scale_fill_brewer(name = "Method", labels = c_methods_wgbs,palette="Dark2",drop = F) + 
  scale_x_discrete( drop=FALSE) + scale_y_log10()+
  ggtitle("Signature length, NGS")

p22 <- ggplot(H_all_rrbs, aes(x=method, y = length, fill=method))  +
  facet_wrap(~ dataset, scales = "free",  ncol = 1) +
  geom_bar(stat="identity", position=position_dodge())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), legend.position="bottom", 
        text = element_text(size=20, face = "bold")) +
  labs(x ='Method', y = "log(Number of DMCs)") + 
  scale_fill_manual(name = "Method", labels = c_methods_rrbs,values = fixed_dark2 ,drop = F) + 
  scale_x_discrete( drop=FALSE) + scale_y_log10() +
  ggtitle("")

p2 <- ggarrange(p21, p22, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

########## microarray
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
c_methods <-c("limma","Welch t-test", "minfi (dmpFinder)", "minfi (dmpFinder, variance shrinkage)")
method_names <- c(
  'limma' = "limma",
  'ttest' = "Welch T-test",
  'dmpFinder' = "minfi (dmpFinder)",
  'dmpFinder_shrinkVar' = "minfi (dmpFinder, variance shrinkage)")

H_all <- read.table("H_full_microarray", sep='\t', header = T)
H_all$length <- as.numeric(H_all$length)
H_all$result_prefix<-factor(H_all$result_prefix,levels = c_datasets)
H_all$method<-factor(H_all$method, levels = c_methods_ids)

p1 <- ggplot(H_all, aes(x=method, y = H_all$length, fill=method))  +
  facet_wrap(~ result_prefix,  ncol = 4, labeller = as_labeller(dataset_names)) + # scales = "free",
  geom_bar(stat="identity", position=position_dodge())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), legend.position="bottom", 
        text = element_text(size=20, face = "bold")) +
  labs(x ='Method', y = "log(Number of DMCs)") + 
  scale_fill_brewer(name = "Method", labels = c_methods, palette="Set1",drop = F) + 
  scale_x_discrete( drop=FALSE) + scale_y_log10()+
  ggtitle("Signature length, microarray")


p3 <- ggarrange(p1, p2, ncol=1, nrow=2, common.legend = FALSE)
ggsave(paste0("DMC_num_log.png"),  plot = p3,  device = "png", width = 35,  height = 40,  units = "cm")
