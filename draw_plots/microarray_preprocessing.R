library (ggplot2)

datasets <- c('GSE210301', 'GSE175458', 'GSE180355', 'GSE174525', 'GSE156994', 'GSE145361', 'GSE60185', 'GSE136319', 'GSE157341', 'GSE49149', 'GSE101764', 'GSE85845')
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

## Sample size stack
sample_size <- read.csv(file = 'group_size_microarray_new', sep='\t')
sample_size$datasets <-factor(sample_size$dataset, levels = c_datasets)
sample_size$group <-factor(sample_size$group, levels = c('control', 'case'))

p <- ggplot(sample_size, aes(x=datasets , y=group_size  , fill=group)) +
  geom_bar(stat = "identity", position = 'stack')+
  labs(x ='Dataset', y = "log(Group size)", fill = 'Group' ) +
  scale_x_discrete(drop=FALSE, labels = c_datasets_names) +
  scale_fill_brewer( palette="Set1", drop = F)+ 
  theme(text = element_text(size=20, face = "bold"), axis.text.x = element_text(angle = 90)) + scale_y_log10()

ggsave(paste0("sample_num_microarray.png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")
  
## Cpg num before and after filtering
cpg_num_microarray <- read.csv(file = 'cpg_num_microarray_new', sep='\t')
cpg_num_microarray$datasets <-factor(cpg_num_microarray$dataset, levels = c_datasets)
cpg_num_microarray$filtering <-factor(cpg_num_microarray$filtering, levels = c('after_filtering',
                                                                              'before_filtering'))
filtering_labels = c('After filtering', 'Before filtering')

p <- ggplot(cpg_num_microarray,aes(x=datasets, y = cpg_num , fill=filtering)) +
  geom_bar(stat = "identity", position = 'stack') +
  scale_x_discrete(drop=FALSE, labels = c_datasets_names) +
  scale_fill_brewer( palette="Set2", drop = F, labels = filtering_labels) +
  labs(x ='Dataset', y = "Number of probes", fill = 'State') +
  theme(text = element_text(size=20, face = "bold"), axis.text.x = element_text(angle = 90))
ggsave(paste0("cpg_num_microarray.png"),  plot = p,  device = "png", width = 30,  height = 20,  units = "cm")
