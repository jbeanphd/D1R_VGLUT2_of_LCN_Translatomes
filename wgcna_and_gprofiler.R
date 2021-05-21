# Jonathan Bean, PhD
# 2021-05-20

#libraries
libs <- c("WGCNA", "limma", "gplots",
          "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr", "gprofiler2", "extrafont", "forcats")

lapply(libs, require, character.only = TRUE)

#import data
# genes read group tracking (rgt)
rgt <- read.delim('genes.read_group_tracking')

# gene exp diff (stats)
stats <- read.delim('gene_exp.diff')

# parse stats into meaningful groups
d1r_stats <- filter(stats, sample_1=='Drd1_N', sample_2=='Drd1_IP')
vg2_stats <- filter(stats, sample_1=='Vglut2_N', sample_2=='Vglut2_IP')
vg2ip_v_d1rip_stats <- filter(stats, sample_1=='Vglut2_IP', sample_2=='Drd1_IP')
vg2n_v_d1rn_stats <- filter(stats, sample_1=='Vglut2_N', sample_2=='Drd1_N')


# for WGCNA analysis 
# using normalized frags 
# and giving each sample a unique name
rgt_matrix_nfrags <- rgt[, c(1:3,6)] %>% unite(condition2, c(condition, replicate)) %>% spread(key = condition2, value = external_scaled_frags)
rgt_matrix_nfrags2 <- rgt_matrix_nfrags
rgt_matrix_nfrags <- rgt_matrix_nfrags[,-1]
row.names(rgt_matrix_nfrags) <- rgt_matrix_nfrags2[,1]
t_rgt_matrix_nfrags <- t(rgt_matrix_nfrags)


# Code borrowed from tutorials found at:
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
# Some variable names were changed

# remove genes with no counts on any sample and other troublesome genes 
gsg_nfrags = goodSamplesGenes(t_rgt_matrix_nfrags, verbose = 3);
gsg_nfrags$allOK


if (!gsg_nfrags$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg_nfrags$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(t_rgt_matrix_nfrags)[!gsg_nfrags$goodGenes], collapse = ", ")));
  if (sum(!gsg_nfrags$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(t_rgt_matrix_nfrags)[!gsg_nfrags$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  t_rgt_matrix_nfrags = t_rgt_matrix_nfrags[gsg_nfrags$goodSamples, gsg_nfrags$goodGenes]
}



sampleTree_nfrags = hclust(dist(t_rgt_matrix_nfrags), method = "average") 

plot(sampleTree_nfrags, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


powers <- c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft_nfrags = pickSoftThreshold(t_rgt_matrix_nfrags, powerVector = powers, verbose = 5)


# plot to pick a power
plot(sft_nfrags$fitIndices[,1], -sign(sft_nfrags$fitIndices[,3])*sft_nfrags$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_nfrags$fitIndices[,1], -sign(sft_nfrags$fitIndices[,3])*sft_nfrags$fitIndices[,2],
     labels=powers,cex=1,col="red");


plot(sft_nfrags$fitIndices[,1], sft_nfrags$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft_nfrags$fitIndices[,1], sft_nfrags$fitIndices[,5], labels=powers, cex=1,col="red")




# block-wise modules for all conditions at one time
# use power from picksoftthreshold
bwm_all <- blockwiseModules(t_rgt_matrix_nfrags, power = 10, TOMType = 'signed', networkType = 'signed',
                            minModuleSize = 100, reassignThreshold = 0, mergeCutHeight = 0.2,
                            detectCutHeight = 0.99, minKMEtoStay = 0.1, corType = 'bicor',
                            numericLabels = TRUE, pamStage = TRUE, maxBlockSize = 2500,
                            saveTOMs = TRUE, saveTOMFileBase = 'bwm_all', verbose = 10)




modLabs_bwm_all <- bwm_all$colors
modCols_bwm_all <- labels2colors(modLabs_bwm_all)
modMems_bwm_all <- data.frame(Gene = colnames(t_rgt_matrix_nfrags), Module =(modCols_bwm_all))



#examine modules
unique(modMems_bwm_all$Module) 
table(bwm_all$colors)

modnumbers <- modMems_bwm_all %>% group_by(Module) %>% count() %>% as.data.frame()
write.csv(modnumbers, file = '../20210415_wgcna_heatmap_etc/modnumbers.csv')


    # examine module relationships
plotDendroAndColors(bwm_all$dendrograms[[1]], modCols_bwm_all[bwm_all$blockGenes[[1]]],
                    "all conditions at one time modules",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


bwm_all_MEs <- c('turquoise','brown','darkturquoise','darkgreen','greenyellow','cyan','royalblue',     
                 'purple','blue','grey','lightcyan','red','salmon','tan','magenta',       
                 'orange','lightyellow','yellow','black','pink',         
                 'darkred','green','midnightblue','lightgreen','grey60','darkgrey') %>% as.list()


bwm_allMEs0 = moduleEigengenes(t_rgt_matrix_nfrags, modCols_bwm_all)$eigengenes
bwm_allorder_MEs = orderMEs(bwm_allMEs0)


color_scale = colorpanel(50, "Blue", "Black", "Yellow")
## label the gene expression level



plotEigengeneNetworks(bwm_allorder_MEs, "using all conditions", marDendro = c(0,4,2,0),
                      heatmapColors = color_scale,
                      signed = TRUE,
                      plotHeatmaps = TRUE)

# end code borrowed from
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html



## custom function to look at ME overlap with enriched genes or DE genes

percent_enr <- function(ME_list, modMem, de_list){
  p_list <- data.frame(ME =factor(levels = ME_list), percent = numeric())
  pos <- 1
  
  for(val in ME_list){
    
    n = filter(modMem, Module == val) %>% inner_join(de_list, by = "Gene") %>% count()
    d = filter(modMem, Module == val) %>% count()
    p = n/d
    p_list[pos,1] = factor(val)
    if(p > 0){
      p_list[pos,2] = p*100
      
    } else {
      p_list[pos,2] = 0
    }
    
    pos = pos + 1
    
    
  } 
  return(p_list)
  
} 


# make list of genes that are enriched or DE
vg2_enr1.5x_Gene <- vg2_stats %>% filter(log2.fold_change. > log2(1.5), q_value < 0.05) %>% 
  rename(Gene=test_id) %>% select(Gene)
d1r_enr1.5x_Gene <- d1r_stats %>% filter(log2.fold_change. > log2(1.5), q_value < 0.05) %>% 
  rename(Gene=test_id) %>% select(Gene)


vg2ip_de1.5x_Gene <- vg2ip_v_d1rip_stats %>% filter(log2.fold_change. < -log2(1.5), q_value < 0.05) %>% 
  rename(Gene=test_id) %>% select(Gene)
d1rip_de1.5x_Gene <- vg2ip_v_d1rip_stats %>% filter(log2.fold_change. > log2(1.5), q_value < 0.05) %>% 
  rename(Gene=test_id) %>% select(Gene)


# call function percent_enr on list
####

vg2_enr1.5x_ME <- percent_enr(bwm_all_MEs, modMems_bwm_all, vg2_enr1.5x_Gene)


d1r_enr1.5x_ME <- percent_enr(bwm_all_MEs, modMems_bwm_all, d1r_enr1.5x_Gene)


vg2ip_de1.5x_ME <- percent_enr(bwm_all_MEs, modMems_bwm_all, vg2ip_de1.5x_Gene)


d1rip_de1.5x_ME <- percent_enr(bwm_all_MEs, modMems_bwm_all, d1rip_de1.5x_Gene)


module_percentage <- modnumbers %>% rename(ME=Module) %>% 
  inner_join(vg2_enr1.5x_ME, by = 'ME') %>% 
  rename(vg2_enr1.5X=percent) %>% 
  inner_join(vg2ip_de1.5x_ME, by = 'ME') %>% 
  rename(vg2_de1.5x=percent) %>% 
  inner_join(d1r_enr1.5x_ME, by = 'ME') %>%
  rename(d1r_enr1.5x=percent) %>% 
  inner_join(d1rip_de1.5x_ME, by = 'ME') %>% 
  rename(d1r_de1.5x=percent)

#write_excel_csv(module_percentage, file = 'module_percentages.xlsx')


###gprofiler###

# between beginning analyses and completing them, gprofiler updated from Ensembl 102 to 103
# according to gprofiler website as of 2021-05-20:
#Since updating to Ensembl 103, part of the Gene Ontology annotations for 
#some species (including Mouse/Mus musculus/mmusculus) were lost. 
#While waiting for the annotations to hopefully reemerge in Ensembl 104, 
#we suggest using our archive of Ensembl 102 if you notice unexpected changes in your results.

#R users can access archived versions using instructions at the gprofiler2 vignette.

# setting base url to Ensemble archive 102 fixes lost annotations from 103

set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e102_eg49_p15")



# turqouise module for d1r+

turquoiseME <-modMems_bwm_all %>% filter(Module=='turquoise') 

# filter for genes in turq module that are both 50% enriched and 50% DE (1.5x) in d1r neurons
turq_top <- turquoiseME %>% inner_join(d1r_enr1.5x_Gene) %>% inner_join(d1rip_de1.5x_Gene)


turq_top_gost <- gost(query = (turq_top)[1] %>% as.list(), organism = 'mmusculus', ordered_query = FALSE, 
                         evcodes = TRUE, correction_method = 'gSCS', domain_scope = 'annotated', exclude_iea = TRUE,
                         sources = c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'REAC', 'TF', 'WP'))

turq_top_gostplot <- gostplot(turq_top_gost, interactive = FALSE, capped = TRUE) +
  scale_y_continuous(limits = c(NA,17), breaks = c(0, 4, 8, 12, 16),
                     labels = c(1, 0.0001, 0.00000001, 0.000000000001, 0.0000000000000001)) +
  labs(y = 'corrected p-value', text = element_text(family = "Arial", size = 16)) +
  theme(text = element_text(family = "Arial", size = 16),
        axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'))

turq_top_gostplot


turq_top_table <- turq_top_gost$result
turq_top_table <- turq_top_table[,c(10,9,11,5,6,3,16)]

#write_excel_csv(turq_top_table, file = '../20210514_RESULTS/turq_d1r_gprofiler.xlsx')

publish_gostplot(turq_top_gostplot, highlight_terms = turq_top_gost$result[c(130,137,146,150,151,103,108,109,116,122,4,7,8,11,14,160,164,336),])



d1r_stats_turq <- d1r_stats %>% rename(Gene=test_id) %>% inner_join(turquoiseME[1])
d1r_stats_minus_turq <- d1r_stats %>% rename(Gene=test_id) %>% anti_join(turquoiseME[1]) %>% filter(status=='OK')
d1r_stats_turqDE <- d1r_stats_turq %>% inner_join(d1r_enr1.5x_Gene)
d1r_stats_turq_nonDE <- d1r_stats_turq %>% anti_join(d1r_stats_turqDE) %>% filter(status=="OK")


ggplot(d1r_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = d1r_stats_minus_turq, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='white', shape=21, size=3) +
  geom_point(inherit.aes = FALSE, data = d1r_stats_turq_nonDE,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='turquoise', shape=21, alpha=0.5, size=3) +
  geom_point(inherit.aes = FALSE, data = d1r_stats_turqDE,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'turquoise', shape=21, alpha=0.8, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = log2(1.5), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 4.5), breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "Input  >   log2(fold)   <   IP", y = "corrected p-value", title = 'D1R+ Neurons', family = 'Arial') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

#ggsave(path = '../20210411_turqME/', filename = 'd1r_neurons_turq_volc.tif', 
#       width = 6.7, height = 3, units = 'in', device = 'tiff', dpi = 150)




ip_stats_turq <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% inner_join(turquoiseME[1])
ip_stats_minus_turq <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% anti_join(turquoiseME[1]) %>% filter(status=='OK')
ip_stats_turqDE <- ip_stats_turq %>% filter(log2.fold_change.>log2(1.5), q_value<0.05) 
ip_stats_turq_nonDE <- ip_stats_turq %>% anti_join(ip_stats_turqDE) %>% filter(status=='OK')


ggplot(vg2ip_v_d1rip_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = ip_stats_minus_turq, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='white', shape=21, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_turq_nonDE,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='turquoise', shape=21, alpha=0.5, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_turqDE,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'turquoise', shape=21, alpha=0.8, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = log2(1.5), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 10), breaks = c(-4, -2, 0, 2, 4, 6, 8, 10)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "VGLUT2+  >   log2(fold)   <   D1R+", y = "corrected p-value", title = 'Immunoprecipitants', family = 'Arial') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

#ggsave(path = '../20210411_turqME/', filename = 'ip_v_ip_turq_volc.tif', 
#       width = 6.7, height = 3, units = 'in', device = 'tiff', dpi = 150)






## midnight blue - module for d1r+
midnightblueME <-modMems_bwm_all %>% filter(Module=='midnightblue')
mnblued1rtop <- midnightblueME[1] %>% inner_join(d1r_enr1.5x_Gene) %>% inner_join(d1rip_de1.5x_Gene)

mnblued1rtop_gost <- gost(query = mnblued1rtop %>% as.list(), organism = 'mmusculus', ordered_query = FALSE, 
                          evcodes = TRUE, correction_method = 'gSCS', domain_scope = 'annotated', exclude_iea = TRUE,
                          sources = c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'REAC', 'TF', 'WP'))






mnblued1rtop_gostplot <- gostplot(mnblued1rtop_gost, interactive = FALSE, capped = FALSE) +
  scale_y_continuous(limits = c(NA,13), breaks = c(0, 2, 4, 6, 8, 10, 12),
                     labels = c(1, 0.01, 0.0001, 0.000001, 0.00000001, 0.0000000001, 0.000000000001)) +
  labs(y = 'corrected p-value', text = element_text(family = "Arial", size = 16)) +
  theme(text = element_text(family = "Arial", size = 16),
        axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'))

mnblued1rtop_gostplot

publish_gostplot(mnblued1rtop_gostplot, highlight_terms = mnblued1rtop_gost$result[c(118,123,131,133,136,87,91,93,94,100,3,9,10,15,19,140,141,142),])


mnblued1rtop_table <- mnblued1rtop_gost$result
mnblued1rtop_table <- mnblued1rtop_table[,c(10,9,11,5,6,3,16)]

#write_excel_csv(mnblued1rtop_table, file = '../20210514_RESULTS/mnblue_d1r_gprofiler.xlsx')



d1r_stats_mnblue <- d1r_stats %>% rename(Gene=test_id) %>% inner_join(midnightblueME[1])
d1r_stats_minus_mnblue <- d1r_stats %>% rename(Gene=test_id) %>% anti_join(midnightblueME[1]) %>% filter(status=='OK')
d1r_stats_mnblueDE <- d1r_stats_mnblue %>% inner_join(d1r_enr1.5x_Gene)
d1r_stats_mnblue_nonDE <- d1r_stats_mnblue %>% anti_join(d1r_stats_mnblueDE) %>% filter(status=="OK")


ggplot(d1r_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = d1r_stats_minus_mnblue, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='white', shape=21, size=3) +
  geom_point(inherit.aes = FALSE, data = d1r_stats_mnblue_nonDE,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='midnightblue', shape=21, alpha=0.5, size=3) +
  geom_point(inherit.aes = FALSE, data = d1r_stats_mnblueDE,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'midnightblue', shape=21, alpha=0.8, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = log2(1.5), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 4.5), breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "Input  >   log2(fold)   <   IP", y = "corrected p-value", title = 'D1R+ Neurons', family = 'Arial') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

#ggsave(path = '../20210406_wgcna_gprofiler/', filename = 'd1r_neurons_mnblue_volc.tif', 
#       width = 6.7, height = 3, units = 'in', device = 'tiff', dpi = 150)




ip_stats_mnblue <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% inner_join(midnightblueME[1])
ip_stats_minus_mnblue <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% anti_join(midnightblueME[1]) %>% filter(status=='OK')
ip_stats_mnblueDE <- ip_stats_mnblue %>% filter(log2.fold_change.>log2(1.5), q_value<0.05) 
ip_stats_mnblue_nonDE <- ip_stats_mnblue %>% anti_join(ip_stats_mnblueDE) %>% filter(status=='OK')


ggplot(vg2ip_v_d1rip_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = ip_stats_minus_mnblue, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='white', shape=21, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_mnblue_nonDE,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='midnightblue', shape=21, alpha=0.5, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_mnblueDE,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'midnightblue', shape=21, alpha=0.8, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = log2(1.5), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 10), breaks = c(-4, -2, 0, 2, 4, 6, 8, 10)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "VGLUT2+  >   log2(fold)   <   D1R+", y = "corrected p-value", title = 'Immunoprecipitants', family = 'Arial') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

#ggsave(path = '../20210406_wgcna_gprofiler/', filename = 'ip_v_ip_mnblue_volc.tif', 
#       width = 6.7, height = 3, units = 'in', device = 'tiff', dpi = 150)



# tan module for d1r 

tanME <-modMems_bwm_all %>% filter(Module=='tan') 

tan_d1rtop <- tanME %>% inner_join(d1r_enr1.5x_Gene) %>% inner_join(d1rip_de1.5x_Gene)


tan_top_gost <- gost(query = (tan_d1rtop)[1] %>% as.list(), organism = 'mmusculus', ordered_query = FALSE, 
                       evcodes = TRUE, correction_method = 'gSCS', domain_scope = 'annotated', exclude_iea = TRUE,
                       sources = c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'REAC', 'TF', 'WP'))

tan_top_gostplot <- gostplot(tan_top_gost, interactive = FALSE, capped = TRUE) +
  scale_y_continuous(limits = c(NA,11), breaks = c(0,2,4,6, 8,10),
                     labels = c(1, 0.01, 0.0001, 0.000001, 0.00000001, 0.0000000001)) +
  labs(y = 'corrected p-value', text = element_text(family = "Arial", size = 16)) +
  theme(text = element_text(family = "Arial", size = 16),
        axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'))

tan_top_gostplot


tan_top_table <- tan_top_gost$result
tan_top_table <- tan_top_table[,c(10,9,11,5,6,3,16)]

#write_excel_csv(tan_top_table, file='../20210514_RESULTS/tan_d1r_gprofiler.xlsx')

publish_gostplot(tan_top_gostplot, highlight_terms = tan_top_gost$result[c(86,87,89,60,78,2,6,17,95,99,100,101,104,107,109,115,117,126),])



d1r_stats_tan <- d1r_stats %>% rename(Gene=test_id) %>% inner_join(tanME[1])
d1r_stats_minus_tan <- d1r_stats %>% rename(Gene=test_id) %>% anti_join(tanME[1]) %>% filter(status=='OK')
d1r_stats_tanDE <- d1r_stats_tan %>% inner_join(d1r_enr1.5x_Gene)
d1r_stats_tan_nonDE <- d1r_stats_tan %>% anti_join(d1r_stats_tanDE) %>% filter(status=="OK")


ggplot(d1r_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = d1r_stats_minus_tan, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='white', shape=21, size=3) +
  geom_point(inherit.aes = FALSE, data = d1r_stats_tan_nonDE,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='tan', shape=21, alpha=0.5, size=3) +
  geom_point(inherit.aes = FALSE, data = d1r_stats_tanDE,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'tan', shape=21, alpha=0.8, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = log2(1.5), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 4.5), breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "Input  >   log2(fold)   <   IP", y = "corrected p-value", title = 'D1R+ Neurons', family = 'Arial') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

#ggsave(path = '../20210411_tan_d1r/', filename = 'd1r_neurons_tan_volc.tif', 
#       width = 6.7, height = 3, units = 'in', device = 'tiff', dpi = 150)




ip_stats_tan <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% inner_join(tanME[1])
ip_stats_minus_tan <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% anti_join(tanME[1]) %>% filter(status=='OK')
ip_stats_tanDE <- ip_stats_tan %>% filter(log2.fold_change.>log2(1.5), q_value<0.05) 
ip_stats_tan_nonDE <- ip_stats_tan %>% anti_join(ip_stats_tanDE) %>% filter(status=='OK')


ggplot(vg2ip_v_d1rip_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = ip_stats_minus_tan, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='white', shape=21, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_tan_nonDE,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='tan', shape=21, alpha=0.5, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_tanDE,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'tan', shape=21, alpha=0.8, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = log2(1.5), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 10), breaks = c(-4, -2, 0, 2, 4, 6, 8, 10)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "VGLUT2+  >   log2(fold)   <   D1R+", y = "corrected p-value", title = 'Immunoprecipitants', family = 'Arial') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

#ggsave(path = '../20210411_tan_d1r/', filename = 'ip_v_ip_tan_volc.tif', 
#       width = 6.7, height = 3, units = 'in', device = 'tiff', dpi = 150)


#### Tan module highlighting enriched in both
# excluding genes highlighted in first tan analysis

tan_bothtop <- tanME %>% inner_join(vg2_enr1.5x_Gene) %>% inner_join(d1r_enr1.5x_Gene) %>% 
  anti_join(d1rip_de1.5x_Gene)


tan_top_both_gost <- gost(query = (tan_bothtop)[1] %>% as.list(), organism = 'mmusculus', ordered_query = FALSE, 
                         evcodes = TRUE, correction_method = 'gSCS', domain_scope = 'annotated', exclude_iea = TRUE,
                         sources = c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'REAC', 'TF', 'WP'))

tan_top_both_gostplot <- gostplot(tan_top_both_gost, interactive = FALSE, capped = TRUE) +
  scale_y_continuous(limits = c(NA,11), breaks = c(0,2,4,6, 8,10),
                     labels = c(1, 0.01, 0.0001, 0.000001, 0.00000001, 0.0000000001)) +
  labs(y = 'corrected p-value', text = element_text(family = "Arial", size = 16)) +
  theme(text = element_text(family = "Arial", size = 16),
        axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'))

tan_top_both_gostplot


tan_top_both_table <- tan_top_both_gost$result
tan_top_both_table <- tan_top_both_table[,c(10,9,11,5,6,3,16)]

#write_excel_csv(tan_top_both_table, file = 'tan_both_gprofiler.xlsx')

publish_gostplot(tan_top_both_gostplot, highlight_terms = tan_top_both_gost$result[c(33,34,22,30,31,3,8,11),])



vg2_stats_tan <- vg2_stats %>% rename(Gene=test_id) %>% inner_join(tanME[1])
vg2_stats_minus_tan <- vg2_stats %>% rename(Gene=test_id) %>% anti_join(tanME[1]) %>% filter(status=='OK')
vg2_stats_tanDE <- vg2_stats_tan %>% inner_join(vg2_enr1.5x_Gene)
vg2_stats_tan_nonDE <- vg2_stats_tan %>% anti_join(vg2_stats_tanDE) %>% filter(status=="OK")


ggplot(vg2_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_minus_tan, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='white', shape=21, size=3) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_tan_nonDE,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='tan', shape=21, alpha=0.5, size=3) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_tanDE,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'tan', shape=21, alpha=0.8, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = log2(1.5), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 4.5), breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "Input  >   log2(fold)   <   IP", y = "corrected p-value", title = 'VGLUT2+ Neurons', family = 'Arial') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

#ggsave(path = '../20210412_tan_vg2/', filename = 'vg2_neurons_tan_volc.tif', 
#       width = 6.7, height = 3, units = 'in', device = 'tiff', dpi = 150)




#ip_stats_tan <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% inner_join(tanME[1])
#ip_stats_minus_tan <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% anti_join(tanME[1]) %>% filter(status=='OK')
ip_stats_tanDE_vg2 <- ip_stats_tan %>% inner_join(vg2_enr1.5x_Gene) %>% anti_join(ip_stats_tanDE[1]) 
ip_stats_tan_nonDE_vg2 <- ip_stats_tan %>% anti_join(ip_stats_tanDE_vg2) %>% filter(status=='OK')


ggplot(vg2ip_v_d1rip_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = ip_stats_minus_tan, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='white', shape=21, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_tan_nonDE_vg2,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='tan', shape=21, alpha=0.5, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_tanDE_vg2,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'tan', shape=21, alpha=0.8, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = log2(1.5), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 10), breaks = c(-4, -2, 0, 2, 4, 6, 8, 10)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "VGLUT2+  >   log2(fold)   <   D1R+", y = "corrected p-value", title = 'Immunoprecipitants', family = 'Arial') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

#ggsave(path = '../20210412_tan_vg2/', filename = 'ip_v_ip_tan_volc.tif', 
#       width = 6.7, height = 3, units = 'in', device = 'tiff', dpi = 150)



# black module for vglut2+

blackME <-modMems_bwm_all %>% filter(Module=='black')
blackvg2top <- blackME[1] %>% inner_join(vg2_enr1.5x_Gene) %>% inner_join(vg2ip_de1.5x_Gene)

blackvg2top_gost <- gost(query = blackvg2top %>% as.list(), organism = 'mmusculus', ordered_query = FALSE, 
                          evcodes = TRUE, correction_method = 'gSCS', domain_scope = 'annotated', exclude_iea = TRUE,
                          sources = c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'REAC', 'TF', 'WP'))



blackvg2top_gostplot <- gostplot(blackvg2top_gost, interactive = FALSE, capped = TRUE) +
  scale_y_continuous(limits = c(NA,17), breaks = c(0, 4, 8, 12, 16),
                     labels = c(1, 0.0001, 0.00000001, 0.000000000001, 0.0000000000000001)) +
  labs(y = 'corrected p-value', text = element_text(family = "Arial", size = 16)) +
  theme(text = element_text(family = "Arial", size = 16),
        axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'))

blackvg2top_gostplot

publish_gostplot(blackvg2top_gostplot, highlight_terms = blackvg2top_gost$result[c(139,144,86,93,95,106,119,8,22,25,28,47,161,163,164,165,166,189),])



blackvg2top_table <- blackvg2top_gost$result
blackvg2top_table <- blackvg2top_table[,c(10,9,11,5,6,3,16)]

#write_excel_csv(blackvg2top_table, file = '../20210514_RESULTS/black_vg2_gprofiler.xlsx')


vg2_stats_black <- vg2_stats %>% rename(Gene=test_id) %>% inner_join(blackME[1])
vg2_stats_minus_black <- vg2_stats %>% rename(Gene=test_id) %>% anti_join(blackME[1]) %>% filter(status=='OK')
vg2_stats_blackDE <- vg2_stats_black %>% inner_join(vg2_enr1.5x_Gene)
vg2_stats_black_nonDE <- vg2_stats_black %>% anti_join(vg2_stats_blackDE) %>% filter(status=="OK")


ggplot(vg2_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_minus_black, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='white', shape=21, size=3) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_black_nonDE,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='black', shape=21, alpha=0.5, size=3) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_blackDE,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'black', shape=21, alpha=0.8, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = log2(1.5), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 4.5), breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "Input  >   log2(fold)   <   IP", y = "corrected p-value", title = 'VGLUT2+ Neurons', family = 'Arial') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

#ggsave(path = '../20210412_blackME_vg2/', filename = 'vg2_neurons_black_volc.tif', 
#       width = 6.7, height = 3, units = 'in', device = 'tiff', dpi = 150)




ip_stats_black <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% inner_join(blackME[1])
ip_stats_minus_black <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% anti_join(blackME[1]) %>% filter(status=='OK')
ip_stats_blackDE_vg2 <- ip_stats_black %>% filter(log2.fold_change.< - log2(1.5), q_value<0.05)  
ip_stats_black_nonDE_vg2 <- ip_stats_black %>% anti_join(ip_stats_blackDE_vg2) %>% filter(status=='OK')


ggplot(vg2ip_v_d1rip_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = ip_stats_minus_black, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='white', shape=21, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_black_nonDE_vg2,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='black', shape=21, alpha=0.5, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_blackDE_vg2,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'black', shape=21, alpha=0.8, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = -log2(1.5), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 10), breaks = c(-4, -2, 0, 2, 4, 6, 8, 10)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "VGLUT2+  >   log2(fold)   <   D1R+", y = "corrected p-value", title = 'Immunoprecipitants', family = 'Arial') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

#ggsave(path = '../20210412_blackME_vg2/', filename = 'ip_v_ip_black_volc.tif', 
#       width = 6.7, height = 3, units = 'in', device = 'tiff', dpi = 150)




# lightgreen module for vglut2+

lightgreenME <-modMems_bwm_all %>% filter(Module=='lightgreen')
lightgreenvg2top <- lightgreenME[1] %>% inner_join(vg2_enr1.5x_Gene) %>% inner_join(vg2ip_de1.5x_Gene)

lightgreenvg2top_gost <- gost(query = lightgreenvg2top %>% as.list(), organism = 'mmusculus', ordered_query = FALSE, 
                         evcodes = TRUE, correction_method = 'gSCS', domain_scope = 'annotated', exclude_iea = TRUE,
                         sources = c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'REAC', 'TF', 'WP'))



lightgreenvg2top_gostplot <- gostplot(lightgreenvg2top_gost, interactive = FALSE, capped = TRUE) +
  scale_y_continuous(limits = c(NA,17), breaks = c(0, 4, 8, 12, 16),
                     labels = c(1, 0.0001, 0.00000001, 0.000000000001, 0.0000000000000001)) +
  labs(y = 'corrected p-value', text = element_text(family = "Arial", size = 16)) +
  theme(text = element_text(family = "Arial", size = 16),
        axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'))

lightgreenvg2top_gostplot

publish_gostplot(lightgreenvg2top_gostplot, highlight_terms = lightgreenvg2top_gost$result[c(139,144,86,93,95,106,119,8,22,25,28,47,161,163,164,165,166,189),])



lightgreenvg2top_table <- lightgreenvg2top_gost$result
lightgreenvg2top_table <- lightgreenvg2top_table[,c(10,9,11,5,6,3,16)]




vg2_stats_lightgreen <- vg2_stats %>% rename(Gene=test_id) %>% inner_join(lightgreenME[1])
vg2_stats_minus_lightgreen <- vg2_stats %>% rename(Gene=test_id) %>% anti_join(lightgreenME[1]) %>% filter(status=='OK')
vg2_stats_lightgreenDE <- vg2_stats_lightgreen %>% inner_join(vg2_enr1.5x_Gene)
vg2_stats_lightgreen_nonDE <- vg2_stats_lightgreen %>% anti_join(vg2_stats_lightgreenDE) %>% filter(status=="OK")


ggplot(vg2_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_minus_lightgreen, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='white', shape=21, size=3) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_lightgreen_nonDE,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='lightgreen', shape=21, alpha=0.5, size=3) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_lightgreenDE,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'lightgreen', shape=21, alpha=0.8, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = log2(1.5), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 4.5), breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "Input  >   log2(fold)   <   IP", y = "corrected p-value", title = 'VGLUT2+ Neurons', family = 'Arial') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

ggsave(path = '../20210412_lightgreenME_vg2/', filename = 'vg2_neurons_lightgreen_volc.tif', 
       width = 6.7, height = 3, units = 'in', device = 'tiff', dpi = 150)




ip_stats_lightgreen <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% inner_join(lightgreenME[1])
ip_stats_minus_lightgreen <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% anti_join(lightgreenME[1]) %>% filter(status=='OK')
ip_stats_lightgreenDE_vg2 <- ip_stats_lightgreen %>% filter(log2.fold_change.< - log2(1.5), q_value<0.05)  
ip_stats_lightgreen_nonDE_vg2 <- ip_stats_lightgreen %>% anti_join(ip_stats_lightgreenDE_vg2) %>% filter(status=='OK')


ggplot(vg2ip_v_d1rip_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = ip_stats_minus_lightgreen, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='white', shape=21, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_lightgreen_nonDE_vg2,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='lightgreen', shape=21, alpha=0.5, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_lightgreenDE_vg2,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'lightgreen', shape=21, alpha=0.8, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = -log2(1.5), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 10), breaks = c(-4, -2, 0, 2, 4, 6, 8, 10)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "VGLUT2+  >   log2(fold)   <   D1R+", y = "corrected p-value", title = 'Immunoprecipitants', family = 'Arial') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

ggsave(path = '../20210412_lightgreenME_vg2/', filename = 'ip_v_ip_lightgreen_volc.tif', 
       width = 6.7, height = 3, units = 'in', device = 'tiff', dpi = 150)





# salmon module for vglut2+

salmonME <-modMems_bwm_all %>% filter(Module=='salmon')
salmonvg2top <- salmonME[1] %>% inner_join(vg2_enr1.5x_Gene) %>% inner_join(vg2ip_de1.5x_Gene)

salmonvg2top_gost <- gost(query = salmonvg2top %>% as.list(), organism = 'mmusculus', ordered_query = FALSE, 
                         evcodes = TRUE, correction_method = 'gSCS', domain_scope = 'annotated', exclude_iea = TRUE,
                         sources = c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'REAC', 'TF', 'WP'))



salmonvg2top_gostplot <- gostplot(salmonvg2top_gost, interactive = FALSE, capped = TRUE) +
  scale_y_continuous(limits = c(NA,17), breaks = c(0, 4, 8, 12, 16),
                     labels = c(1, 0.0001, 0.00000001, 0.000000000001, 0.0000000000000001)) +
  labs(y = 'corrected p-value', text = element_text(family = "Arial", size = 16)) +
  theme(text = element_text(family = "Arial", size = 16),
        axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'))

salmonvg2top_gostplot

publish_gostplot(salmonvg2top_gostplot, highlight_terms = salmonvg2top_gost$result[c(139,144,86,93,95,106,119,8,22,25,28,47,161,163,164,165,166,189),])



salmonvg2top_table <- salmonvg2top_gost$result
salmonvg2top_table <- salmonvg2top_table[,c(10,9,11,5,6,3,16)]




vg2_stats_salmon <- vg2_stats %>% rename(Gene=test_id) %>% inner_join(salmonME[1])
vg2_stats_minus_salmon <- vg2_stats %>% rename(Gene=test_id) %>% anti_join(salmonME[1]) %>% filter(status=='OK')
vg2_stats_salmonDE <- vg2_stats_salmon %>% inner_join(vg2_enr1.5x_Gene)
vg2_stats_salmon_nonDE <- vg2_stats_salmon %>% anti_join(vg2_stats_salmonDE) %>% filter(status=="OK")


ggplot(vg2_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_minus_salmon, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='white', shape=21, size=3) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_salmon_nonDE,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='salmon', shape=21, alpha=0.5, size=3) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_salmonDE,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'salmon', shape=21, alpha=0.8, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = log2(1.5), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 4.5), breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "Input  >   log2(fold)   <   IP", y = "corrected p-value", title = 'VGLUT2+ Neurons', family = 'Arial') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

ggsave(path = '../20210412_salmonME_vg2/', filename = 'vg2_neurons_salmon_volc.tif', 
       width = 6.7, height = 3, units = 'in', device = 'tiff', dpi = 150)




ip_stats_salmon <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% inner_join(salmonME[1])
ip_stats_minus_salmon <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% anti_join(salmonME[1]) %>% filter(status=='OK')
ip_stats_salmonDE_vg2 <- ip_stats_salmon %>% filter(log2.fold_change.< - log2(1.5), q_value<0.05)  
ip_stats_salmon_nonDE_vg2 <- ip_stats_salmon %>% anti_join(ip_stats_salmonDE_vg2) %>% filter(status=='OK')


ggplot(vg2ip_v_d1rip_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = ip_stats_minus_salmon, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='white', shape=21, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_salmon_nonDE_vg2,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='salmon', shape=21, alpha=0.5, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_salmonDE_vg2,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'salmon', shape=21, alpha=0.8, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = -log2(1.5), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 10), breaks = c(-4, -2, 0, 2, 4, 6, 8, 10)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "VGLUT2+  >   log2(fold)   <   D1R+", y = "corrected p-value", title = 'Immunoprecipitants', family = 'Arial') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

ggsave(path = '../20210412_salmonME_vg2/', filename = 'ip_v_ip_salmon_volc.tif', 
       width = 6.7, height = 3, units = 'in', device = 'tiff', dpi = 150)


# lightyellow module for vglut2+

lightyellowME <-modMems_bwm_all %>% filter(Module=='lightyellow')
lightyellowvg2top <- lightyellowME[1] %>% inner_join(vg2_enr1.5x_Gene) %>% inner_join(vg2ip_de1.5x_Gene)

# collapse light green, light yellow, and salmon modules manually 

comb_lylgs <- rbind(lightyellowvg2top, lightgreenvg2top, salmonvg2top)
rainbowME <- rbind(salmonME, lightgreenME, lightyellowME)
rainbow_1.5x_enr <- vg2_stats %>% rename(Gene=test_id) %>% inner_join(rainbowME) %>% 
  filter(log2.fold_change. > log2(1.5), q_value < 0.05)
rainbow_1.5x_de <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% inner_join(rainbowME) %>% 
  filter(log2.fold_change. < -log2(1.5), q_value < 0.05)


  comb_lylgs_vg2top_gost <- gost(query = comb_lylgs %>% as.list(), organism = 'mmusculus', ordered_query = FALSE, 
                          evcodes = TRUE, correction_method = 'gSCS', domain_scope = 'annotated', exclude_iea = TRUE,
                          sources = c('GO:MF', 'GO:CC', 'GO:BP', 'KEGG', 'REAC', 'TF', 'WP'))



comb_lylgs_vg2top_gostplot <- gostplot(comb_lylgs_vg2top_gost, interactive = FALSE, capped = TRUE) +
  scale_y_continuous(limits = c(NA,12), breaks = c(0, 4, 8, 12),
                     labels = c(1, 0.0001, 0.00000001, 0.000000000001)) +
  labs(y = 'corrected p-value', text = element_text(family = "Arial", size = 16)) +
  theme(text = element_text(family = "Arial", size = 16),
        axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'))

comb_lylgs_vg2top_gostplot

publish_gostplot(comb_lylgs_vg2top_gostplot, highlight_terms = comb_lylgs_vg2top_table$term_id[c(90,91,94,96,97,54,60,61,63,65,4,12,21,25,29,106,150)])



comb_lylgs_vg2top_table <- comb_lylgs_vg2top_gost$result
comb_lylgs_vg2top_table <- comb_lylgs_vg2top_table[,c(10,9,11,5,6,3,16)]

#write_excel_csv(comb_lylgs_vg2top_table, file = '../20210514_RESULTS/lylgs_vg2_gprofiler.xlsx')



vg2_stats_salmon <- vg2_stats %>% rename(Gene=test_id) %>% inner_join(salmonME[1])
vg2_stats_minus_salmon <- vg2_stats %>% rename(Gene=test_id) %>% anti_join(salmonME[1]) %>% filter(status=='OK')
vg2_stats_salmonDE <- vg2_stats_salmon %>% inner_join(vg2_enr1.5x_Gene)
vg2_stats_salmon_nonDE <- vg2_stats_salmon %>% anti_join(vg2_stats_salmonDE) %>% filter(status=="OK")

vg2_stats_lightyellow <- vg2_stats %>% rename(Gene=test_id) %>% inner_join(lightyellowME[1])
vg2_stats_minus_lightyellow <- vg2_stats %>% rename(Gene=test_id) %>% anti_join(lightyellowME[1]) %>% filter(status=='OK')
vg2_stats_lightyellowDE <- vg2_stats_lightyellow %>% inner_join(vg2_enr1.5x_Gene)
vg2_stats_lightyellow_nonDE <- vg2_stats_lightyellow %>% anti_join(vg2_stats_lightyellowDE) %>% filter(status=="OK")

vg2_stats_lightgreen <- vg2_stats %>% rename(Gene=test_id) %>% inner_join(lightgreenME[1])
vg2_stats_minus_lightgreen <- vg2_stats %>% rename(Gene=test_id) %>% anti_join(lightgreenME[1]) %>% filter(status=='OK')
vg2_stats_lightgreenDE <- vg2_stats_lightgreen %>% inner_join(vg2_enr1.5x_Gene)
vg2_stats_lightgreen_nonDE <- vg2_stats_lightgreen %>% anti_join(vg2_stats_lightgreenDE) %>% filter(status=="OK")

vg2_stats_minus_rainbow <- vg2_stats %>% rename(Gene=test_id) %>% anti_join(vg2_stats_salmon) %>% 
  anti_join(vg2_stats_lightyellow) %>% anti_join(vg2_stats_lightgreen) %>% filter(status=='OK')


ggplot(vg2_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_minus_rainbow, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='white', shape=21, size=3) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_salmon_nonDE,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='salmon', shape=21, alpha=0.7, size=3) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_salmonDE,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'salmon', shape=21, alpha=0.7, size=3) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_lightgreen_nonDE,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='lightgreen', shape=21, alpha=0.7, size=3) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_lightgreenDE,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'lightgreen', shape=21, alpha=0.7, size=3) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_lightyellow_nonDE,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='lightyellow', shape=21, alpha=0.7, size=3) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_lightyellowDE,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'lightyellow', shape=21, alpha=0.7, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = log2(1.5), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 4.5), breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "Input  >   log2(fold)   <   IP", y = "corrected p-value", title = 'VGLUT2+ Neurons', family = 'Arial') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

#ggsave(path = '../20210415_rainbow_vg2/', filename = 'vg2_neurons_rainbow_volc.tif', 
#       width = 6.7, height = 3, units = 'in', device = 'tiff', dpi = 150)




ip_stats_salmon <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% inner_join(salmonME[1])
ip_stats_minus_salmon <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% anti_join(salmonME[1]) %>% filter(status=='OK')
ip_stats_salmonDE_vg2 <- ip_stats_salmon %>% filter(log2.fold_change.< - log2(1.5), q_value<0.05)  
ip_stats_salmon_nonDE_vg2 <- ip_stats_salmon %>% anti_join(ip_stats_salmonDE_vg2) %>% filter(status=='OK')

ip_stats_lightgreen <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% inner_join(lightgreenME[1])
ip_stats_minus_lightgreen <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% anti_join(lightgreenME[1]) %>% filter(status=='OK')
ip_stats_lightgreenDE_vg2 <- ip_stats_lightgreen %>% filter(log2.fold_change.< - log2(1.5), q_value<0.05)  
ip_stats_lightgreen_nonDE_vg2 <- ip_stats_lightgreen %>% anti_join(ip_stats_lightgreenDE_vg2) %>% filter(status=='OK')

ip_stats_lightyellow <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% inner_join(lightyellowME[1])
ip_stats_minus_lightyellow <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% anti_join(lightyellowME[1]) %>% filter(status=='OK')
ip_stats_lightyellowDE_vg2 <- ip_stats_lightyellow %>% filter(log2.fold_change.< - log2(1.5), q_value<0.05)  
ip_stats_lightyellow_nonDE_vg2 <- ip_stats_lightyellow %>% anti_join(ip_stats_lightyellowDE_vg2) %>% filter(status=='OK')

ip_stats_minus_rainbow <- vg2ip_v_d1rip_stats %>% rename(Gene=test_id) %>% anti_join(salmonME) %>% 
  anti_join(lightgreenME) %>% anti_join(lightyellowME)


ggplot(vg2ip_v_d1rip_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = ip_stats_minus_salmon, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='white', shape=21, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_salmon_nonDE_vg2,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='salmon', shape=21, alpha=0.7, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_salmonDE_vg2,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'salmon', shape=21, alpha=0.7, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_lightgreen_nonDE_vg2,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='lightgreen', shape=21, alpha=0.7, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_lightgreenDE_vg2,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'lightgreen', shape=21, alpha=0.7, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_lightyellow_nonDE_vg2,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='lightyellow', shape=21, alpha=0.7, size=3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_lightyellowDE_vg2,
             aes(log2.fold_change., -log10(q_value)),
             color='red', fill = 'lightyellow', shape=21, alpha=0.7, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = -log2(1.5), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 10), breaks = c(-4, -2, 0, 2, 4, 6, 8, 10)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "VGLUT2+  >   log2(fold)   <   D1R+", y = "corrected p-value", title = 'Immunoprecipitants', family = 'Arial') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

#ggsave(path = '../20210415_rainbow_vg2//', filename = 'ip_v_ip_rainbow_volc.tif', 
#       width = 6.7, height = 3, units = 'in', device = 'tiff', dpi = 150)
