# Jonathan Bean, PhD
# 2021-05-20

#libraries
libs <- c("WGCNA", "limma", "gplots",
          "readr", "magrittr", "purrr", 
          "dplyr", "ggplot2", "tidyr", "gprofiler2", 
          "extrafont", "forcats", "limma", "ggfortify")

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


# for PCA analysis 
# using fpkm 
# and giving each sample a unique name
rgt_matrix_fpkm <- rgt[, c(1:3,7)] %>% unite(condition2, c(condition, replicate)) %>% spread(key = condition2, value = FPKM)
rgt_matrix_fpkm2 <- rgt_matrix_fpkm
rgt_matrix_fpkm <- rgt_matrix_fpkm[,-1]
row.names(rgt_matrix_fpkm) <- rgt_matrix_fpkm2[,1]
t_rgt_matrix_fpkm <- t(rgt_matrix_fpkm)
# discard genes without any counts on any of the samples
t_rgt_matrix_fpkm <- t_rgt_matrix_fpkm[, colSums(t_rgt_matrix_fpkm)>0]


# rank. = 100, upto 100 PCs, 24 in this case
pca_fpkm <- prcomp(t_rgt_matrix_fpkm, scale. = TRUE, center = TRUE, rank. = 100)
condition_pca <- c('Drd1_IP','Drd1_IP','Drd1_IP','Drd1_IP','Drd1_IP','Drd1_IP',
                  'Drd1_N','Drd1_N','Drd1_N','Drd1_N','Drd1_N','Drd1_N',
                  'Vglut2_IP','Vglut2_IP','Vglut2_IP','Vglut2_IP','Vglut2_IP','Vglut2_IP',
                  'Vglut2_N','Vglut2_N','Vglut2_N','Vglut2_N','Vglut2_N','Vglut2_N')

pca_fpkm2 <- cbind(condition_pca, t_rgt_matrix_fpkm)



autoplot(pca_fpkm, data = pca_fpkm2, colour = 'condition_pca', frame = FALSE, size = 5, x = 1, y = 2) +
  theme_classic() + scale_color_manual(values = c("#32669D","lightblue","orange","gold")) +
  theme(text = element_text(size = 16, family='Arial'), 
        axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black')) +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# rows 1-6 are D1R IPs
d1r_rgt_matrix_fpkm <- t_rgt_matrix_fpkm[1:6,]
# discard genes without any counts for any of the D1R IPs
d1r_rgt_matrix_fpkm <- d1r_rgt_matrix_fpkm[, colSums(d1r_rgt_matrix_fpkm)>0]


pca_d1r <- prcomp(d1r_rgt_matrix_fpkm, scale. = TRUE, center = TRUE, rank. = 6)
condition_d1r <- c('Drd1_M','Drd1_M','Drd1_M','Drd1_F','Drd1_F','Drd1_F')

pca_d1r2 <- cbind(condition_d1r, d1r_rgt_matrix_fpkm)



autoplot(pca_d1r, data = pca_d1r2, colour = 'condition_d1r', frame = FALSE, size = 5, x = 1, y = 2) +
  theme_classic() + scale_color_manual(values = c("#32669D","salmon")) +
  theme(text = element_text(size = 16, family='Arial'), 
        axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black')) +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())



# rows 13-18 are Vglut2 IPs
vg2_rgt_matrix_fpkm <- t_rgt_matrix_fpkm[13:18,]
#discard genes with 0 counts on any of the Vglut2 IPs
vg2_rgt_matrix_fpkm <- vg2_rgt_matrix_fpkm[, colSums(vg2_rgt_matrix_fpkm)>0]


pca_vg2 <- prcomp(vg2_rgt_matrix_fpkm, scale. = TRUE, center = TRUE, rank. = 6)
condition_vg2 <- c('VG2_M','VG2_M','VG2_M','VG2_F','VG2_F','VG2_F')

pca_vg22 <- cbind(condition_vg2, vg2_rgt_matrix_fpkm)



autoplot(pca_vg2, data = pca_vg22, colour = 'condition_vg2', frame = FALSE, size = 5, x = 1, y = 2) +
  theme_classic() + scale_color_manual(values = c("orange","pink")) +
  theme(text = element_text(size = 16, family='Arial'), 
        axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black')) +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())



# function to rank genes
# rank genes by fold change and q_value
# log2(inf) = 4  or (abs(log2(fc)) > 4) = 4

top_genes <- function(stats_df){
  
  sub_stats_df <- stats_df %>% filter(abs(log2.fold_change.)>1, q_value <0.05)
  
    for(i in 1:length(sub_stats_df[,1])){
      if (abs(sub_stats_df$log2.fold_change.[i]) < 4) {
        sub_stats_df$abs_logFC[i] = abs(sub_stats_df$log2.fold_change.[i])} else {
          sub_stats_df$abs_logFC[i] = 4
        }
      sub_stats_df$rank[i] = sub_stats_df$abs_logFC[i] * -log10(sub_stats_df$q_value[i])
  }
  
  
  
  return(sub_stats_df)
}





top_d1r_n_v_ip <- top_genes(d1r_stats)
# exclude de-enriched genes
top_d1r_n_v_ip <- top_d1r_n_v_ip %>% filter(log2.fold_change.>1)



top_vg2_n_v_ip <- top_genes(vg2_stats)
#exclude de-enriched genes
top_vg2_n_v_ip <- top_vg2_n_v_ip %>% filter(log2.fold_change.>1)



top_d1r_ip_v_ip <- top_genes(vg2ip_v_d1rip_stats)
# select genes differentially expressed in D1R neurons
top_d1r_ip_v_ip <- top_d1r_ip_v_ip %>% filter(log2.fold_change.>1)




top_vg2_ip_v_ip <- top_genes(vg2ip_v_d1rip_stats)
# select genes differentially expressed in Vglut2 neurons
top_vg2_ip_v_ip <- top_vg2_ip_v_ip %>% filter(log2.fold_change. < -1)


# create combined scores for within + between group comparisons
top_d1r_comb <- inner_join(top_d1r_n_v_ip[,c(3,5,6,8:10,13,15,16)], top_d1r_ip_v_ip[,c(3,5,6,8:10,13,15,16)], by = 'gene')
top_d1r_comb <- top_d1r_comb %>% mutate(rank_comb = (rank.x + rank.y))

top_vg2_comb <- inner_join(top_vg2_n_v_ip[,c(3,5,6,8:10,13,15,16)], top_vg2_ip_v_ip[,c(3,5,6,8:10,13,15,16)], by = 'gene')
top_vg2_comb <- top_vg2_comb %>% mutate(rank_comb = (rank.x + rank.y))

top_vg2_comb <- top_vg2_comb[,c(1,2,3,11,4,5,13,6,14,7,15,9,17,18)]
top_d1r_comb <- top_d1r_comb[,c(1,2,3,10,4,5,12,6,14,7,15,9,17,18)]

# write to excell files
write_excel_csv(top_vg2_comb, file = '../20210416_RESULTS/vg2_top_ranked_genes.xlsx')
write_excel_csv(top_d1r_comb, file = '../20210416_RESULTS/d1r_top_ranked_genes.xlsx')





# 2-11 excludes mir411 which has large variance 
# select top 10 genes for each group 
top_genes_for_heatmap <- full_join(top_d1r_comb[c(2:11),], top_vg2_comb[1:10,], by = 'gene') %>% 
  select(gene) %>% rename(tracking_id=gene) %>% as.data.frame()


# create matrix for heatmap of top10 genes
heatmap_fpkms <- inner_join(rgt[, c(1:3,7)], top_genes_for_heatmap) %>% 
  unite(condition2, c(condition, replicate)) %>% spread(key = condition2, value = FPKM)

heatmap_fpkms2 <- heatmap_fpkms
heatmap_fpkms <- heatmap_fpkms[,-1]
row.names(heatmap_fpkms) <- heatmap_fpkms2[,1]
t_heatmap_fpkms <- t(heatmap_fpkms)

heatmap_fpkms <- heatmap_fpkms[,c(19:24,7:18,1:6)]
# log transform due to skewed distribution of d1r ip samples
heatmap_log2fpkm <- log2(heatmap_fpkms+0.01)


# repeating first and last color creates a cap for extreme values
color_scale = colorpanel(50, "Blue", "black", "Yellow")[c(1,1,1,1,1,1,1:15,20,25,30,35:50,50,50,50,50,50,50)]

sizeGrWindow(12,12)

heatmap3::heatmap3(heatmap_log2fpkm, Rowv = NULL, Colv = NA, balanceColor = T, 
                   margins = c(10,10), cexRow = 1, cexCol = 1, 
                   showRowDendro = FALSE, method = 'complete', col = color_scale, useRaster = TRUE)


# parse results to plot in different colors on volcano plot for vglut2 ip vs d1r ip
top_10_vg2 <-top_vg2_comb %>% slice_max(order_by = rank_comb, n=10) %>% select(gene) %>% inner_join(vg2ip_v_d1rip_stats)
top_10_d1r <- top_d1r_comb %>% slice_max(order_by = rank_comb, n=11) %>% select(gene) %>% inner_join(vg2ip_v_d1rip_stats)
ip_stats_2xde <- vg2ip_v_d1rip_stats %>% filter(abs(log2.fold_change.)>1, q_value < 0.05) %>%  
  anti_join(top_10_vg2) %>% anti_join(top_10_d1r)
ip_stats_minus_2xde <- vg2ip_v_d1rip_stats %>% anti_join(ip_stats_2xde) %>% anti_join(top_10_d1r) %>% anti_join(top_10_vg2) %>% 
  filter(status=='OK')


#volcano for 2x and top 10  vglut2 ip vs d1r ip

ggplot(vg2ip_v_d1rip_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = ip_stats_minus_2xde, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='black', shape=21, size=3, alpha = 0.3) +
  geom_point(inherit.aes = FALSE, data = ip_stats_2xde,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='red', shape=21, alpha=0.3, size=3) +
  geom_point(inherit.aes = FALSE, data = top_10_d1r,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill = '#32669D', shape=21, size=3) +
  geom_point(inherit.aes = FALSE, data = top_10_vg2,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill = 'orange', shape=21, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = -log2(2), linetype='dashed', color='black', size=1) +
  geom_vline(xintercept = log2(2), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 10), breaks = c(-4, -2, 0, 2, 4, 6, 8, 10)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "VGLUT2+  >   log2(fold)   <   D1R+", y = "corrected p-value", family = 'Arial', size = 16) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

#ggsave(path = '../20210415_wgcna_heatmap_etc/', filename = 'ip_v_ip_2xde_top10_volc.tif', 
#       width = 6.7, height = 2, units = 'in', device = 'tiff', dpi = 150)


#vglut2 parse results for vglut2 input vs ip
top_10_vg2_vg2_stats <-top_vg2_comb %>% slice_max(order_by = rank_comb, n=10) %>% select(gene) %>% inner_join(vg2_stats)
vg2_stats_2xde <- vg2_stats %>% filter(log2.fold_change.>1, q_value < 0.05) %>%  
  anti_join(top_10_vg2_vg2_stats)
vg2_stats_minus_2xde <- vg2_stats %>% anti_join(vg2_stats_2xde) %>% anti_join(top_10_vg2_vg2_stats) %>% 
  filter(status=='OK')


#volcano for 2x and top 10 for vglut2 input vs ip

ggplot(vg2_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_minus_2xde, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='black', shape=21, size=3, alpha = 0.3) +
  geom_point(inherit.aes = FALSE, data = vg2_stats_2xde,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='red', shape=21, alpha=0.3, size=3) +
  geom_point(inherit.aes = FALSE, data = top_10_vg2_vg2_stats,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill = 'orange', shape=21, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = log2(2), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 4.5), breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "input  >   log2(fold)   <   IP", y = "corrected p-value", family = 'Arial', size = 16) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

#ggsave(path = '../20210415_wgcna_heatmap_etc/', filename = 'vg2_2xde_top10_volc.tif', 
#       width = 6.7, height = 2, units = 'in', device = 'tiff', dpi = 150)


#d1r parse data for d1r input vs ip
top_10_d1r_d1r_stats <-top_d1r_comb %>% slice_max(order_by = rank_comb, n=10) %>% select(gene) %>% inner_join(d1r_stats)
d1r_stats_2xde <- d1r_stats %>% filter(log2.fold_change.>1, q_value < 0.05) %>%  
  anti_join(top_10_d1r_d1r_stats)
d1r_stats_minus_2xde <- d1r_stats %>% anti_join(d1r_stats_2xde) %>% anti_join(top_10_d1r_d1r_stats) %>% 
  filter(status=='OK')


#volcano for 2x and top 10 for d1r input vs ip

ggplot(d1r_stats, aes(log2.fold_change., -log10(q_value))) +
  geom_point(inherit.aes = FALSE, data = d1r_stats_minus_2xde, 
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='black', shape=21, size=3, alpha = 0.3) +
  geom_point(inherit.aes = FALSE, data = d1r_stats_2xde,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill='red', shape=21, alpha=0.3, size=3) +
  geom_point(inherit.aes = FALSE, data = top_10_d1r_d1r_stats,
             aes(log2.fold_change., -log10(q_value)),
             color='black', fill = '#32669D', shape=21, size=3) +
  geom_hline(yintercept = -log10(0.05),linetype = 'dashed', color='black', size=1 ) +
  geom_vline(xintercept = log2(2), linetype='dashed', color='black', size=1) +
  theme_bw() +
  scale_x_continuous(limits = c(-4.5, 4.5), breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0,1,2,3,4), labels = c(1, 0.1, 0.01, 0.001, 0.0001)) +
  labs(x = "input  >   log2(fold)   <   IP", y = "corrected p-value", family = 'Arial', size = 16) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size = 16, family='Arial'), axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
        axis.text.y = element_text(size = 16, family = 'Arial', color = 'black'),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  theme(legend.position="none")

#ggsave(path = '../20210415_wgcna_heatmap_etc/', filename = 'd1r_2xde_top10_volc.tif', 
#       width = 6.7, height = 2, units = 'in', device = 'tiff', dpi = 150)
