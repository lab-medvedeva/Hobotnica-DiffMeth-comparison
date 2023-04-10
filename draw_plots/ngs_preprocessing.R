library (ggplot2)

# Order: 'GSE149608','GSE138598','GSE119980','GSE103886', 'GSE150592', 'GSE148060'
c_datasets <- c('GSE149608','GSE138598','GSE119980','GSE103886', 'GSE150592', 'GSE148060')
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

## Sample size stack
sample_size <- read.csv(file = 'sample_num.csv', sep='\t')
sample_size$dataset<-factor(sample_size$dataset,levels = c_datasets)
sample_size$group <-factor(sample_size$group, levels = c('control', 'case'))
p <- ggplot(sample_size, aes(x=dataset, y=size, fill=group)) +
  geom_bar(stat = "identity", position = 'stack')+
  labs(x ='Dataset', y = "Group size", fill = 'Group' ) +
  scale_x_discrete(drop=FALSE, labels = c_datasets) +
  scale_fill_brewer( palette="Set1", drop = F)+ 
  theme(text = element_text(size=20, face = "bold"))

ggsave(paste0("sample_num.png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")

## Number of reads
total_reads<-read.csv(file = 'total_reads', sep='\t') 
total_reads$dataset<-factor(total_reads$dataset,levels = c_datasets)
p<-ggplot(total_reads, aes(x=dataset, y=total.sequences, fill=dataset)) +  geom_boxplot(show.legend = FALSE) +
  labs(x ='Dataset', y = "Total Sequences") + theme(text = element_text(size=20, face = "bold"))
ggsave(paste0("total_count.png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")

## Aligned reads
aligned_reads<-read.csv(file = 'aligned_reads', sep='\t') 
aligned_reads$dataset<-factor(aligned_reads$dataset,levels = c_datasets)
p<-ggplot(aligned_reads, aes(x=dataset, y=aligned, fill=dataset)) +  geom_boxplot(show.legend = FALSE) +
  labs(x ='Dataset', y = "Aligned Reads (%)") + theme(text = element_text(size=20, face = "bold"))
ggsave(paste0("aligned.png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")

## Deduplication
dedupl<-read.csv(file = 'dedupl', sep='\t') 
dedupl$dataset<-factor(dedupl$dataset,levels = c_datasets)
p<-ggplot(dedupl, aes(x=dataset, y=dedup_leftover, fill=dataset)) +  geom_boxplot(show.legend = FALSE) + 
  scale_x_discrete( drop=FALSE) + scale_fill_discrete(drop = F) +
  labs(x ='Dataset', y = "Deduplication Leftover (%)") +
  theme(text = element_text(size=20, face = "bold"))
ggsave(paste0("dedupl.png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")

## CpG number
cpg_num<-read.csv(file = 'cpg_num', sep='\t') 
cpg_num$dataset<-factor(cpg_num$dataset,levels = c_datasets)
p<-ggplot(cpg_num, aes(x=dataset, y=number, fill=dataset)) +  geom_boxplot(show.legend = FALSE) +
  labs(x ='Dataset', y = "Number of CpGs") +
  theme(text = element_text(size=20, face = "bold"))
ggsave(paste0("cpg_num.png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")

df <- data.frame(cpg_num$dataset, cpg_num$number)
means <- aggregate(.~cpg_num.dataset, data=df, mean)
colnames(means) <- c('dataset', 'num')
sds <- aggregate(.~cpg_num.dataset, data=df, sd)
colnames(sds) <- c('dataset', 'sd')
adding <- merge(x = means, y = sds, by = "dataset", all = TRUE)
adding$method <- 'before filtering'
adding <- adding[,c(2,4,1,3)]

## CpG number filtered
cpgnum_filtered<-read.csv(file = 'cpgnum_filtered', sep=' ')
cpgnum_filtered$dataset<-factor(cpgnum_filtered$dataset,levels = c_datasets)
cpgnum_filtered$method<-factor(cpgnum_filtered$method, levels = c('before filtering', 'all','bsmooth', 'hmm'))
cpgnum_filtered$sd <- NA
cpgnum_filtered <- rbind(cpgnum_filtered, adding)

c_methods_num_filtered <- c("before filtering", "methylSig, DSS, methylkit, RADMeth", "BSmooth", "HMM-DM")
p<-ggplot(cpgnum_filtered,aes(x=dataset, y = num, fill=method)) + 
  geom_bar(stat="identity", position=position_dodge()) +  
  geom_errorbar(aes(ymin=num-sd, ymax=num+sd), width=.2,
                position=position_dodge(.9)) +
  labs(x ='Dataset', y = "Number of CpGs") + 
  scale_fill_brewer(name = "Type", labels = c_methods_num_filtered,palette="Set2",drop = F) +
  scale_x_discrete( drop=FALSE) +
  theme(text = element_text(size=20, face = "bold"), legend.position="bottom")

ggsave(paste0("cpgnum_filtered.png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")
