# Jonathan Bean, PhD
# 2021-05-20


#libraries
libs <- c("WGCNA", "limma", "gplots",
          "readr", "magrittr", "purrr", "dplyr", "ggplot2", "tidyr", "gprofiler2", "extrafont", "forcats")

lapply(libs, require, character.only = TRUE)

# objects created in 'wgcna_and_gprofiler.R' and 'pca_and_top_genes.R' needed. Please run those scripts first. 



#violin maker function

violin_maker <- function(list_of_genes, read_group_tracking){
  for(val in list_of_genes){

    print(
      
      filter(read_group_tracking, tracking_id == val) %>%
        ggplot(aes(condition,log2(FPKM+0.5),fill=condition)) +
        geom_violin() +
        theme(legend.position="none") +
        scale_x_discrete(limits=c("Vglut2_N","Vglut2_IP","Drd1_N","Drd1_IP"),
                         labels = c("input", "IP", "input", "IP")) +
        scale_fill_manual(values = c("Vglut2_N"="white","Vglut2_IP"="orange","Drd1_N"="white","Drd1_IP"="#32669D")) +
        facet_wrap(~tracking_id,  scales = "free") +
        labs(x="VGLUT2+          D1R+") +
        theme_bw(base_family = "Arial") +
        labs(y = "log2(FPKM + 0.5)", text = element_text(family = "Arial", size = 16)) +
        theme(text = element_text(family = "Arial", size = 16), 
              plot.title = element_text(hjust = 0.5), 
              plot.subtitle = element_text(hjust = 0.5), 
              strip.text.x = element_text(family = "Arial", face = "italic", size = 20), 
              axis.text.x = element_text(size = 16, family = 'Arial', color = 'black'),
              axis.text.y = element_text(size = 16, family = 'Arial', color = 'black')) +
        theme(legend.position="none")   +
        geom_point(size=2)
      
    )
  }
  
}



# needed for cbindX used in following function
library(gdata)

# function to unlist 'intersection' variable from gprofiler results per GO
# then look for genes which correspond to
# enriched or DE genes provided as a list

split_then_top <- function(gp_table, GO_list, ranked_genes, how_many){
  top_gene_go <- data.frame(rank = 1:how_many, stringsAsFactors = TRUE)
  for (GO in GO_list) {
    GO <- filter(gp_table, term_id == GO) %>% select(intersection) %>% as.character() %>% 
      strsplit('[,]') %>% unlist() %>% as.data.frame() %>% dplyr::rename(gene='.') %>% 
      inner_join(ranked_genes) %>% 
      slice_max(rank_comb, n = how_many) %>% select(gene) 
    
    top_gene_go <- cbindX(top_gene_go, GO)
    
  }
  colnames(top_gene_go) <- c('rank', GO_list)
  return(top_gene_go)
  
}



#function to create a single list of genes that are enriched and DE from any of the GOs given as list 

split_then_top_list <- function(gp_table, GO_list, ranked_genes, how_many){
  top_gene_go <- data.frame(gene = character(), stringsAsFactors = TRUE)
  for (GO in GO_list) {
    GO <- filter(gp_table, term_id == GO) %>% select(intersection) %>% as.character() %>% 
      strsplit('[,]') %>% unlist() %>% as.data.frame() %>% rename(gene='.') %>% 
      inner_join(ranked_genes) %>% 
      slice_max(rank_comb, n = how_many) %>% select(gene)

    
    top_gene_go <- full_join(top_gene_go, GO)
    
  }
  
  return(top_gene_go)
  
}



# funtion to calculate median and mean fold-change for the purpose of comparison 

fold_change_calc <- function(gene_list, read_group_tracking, optionABC){

    med_fc_gene <- read_group_tracking %>% rename(gene=tracking_id) %>% 
      dplyr::group_by(condition, gene) %>% 
      dplyr::summarise(med = median(FPKM)) %>% inner_join(gene_list)
  
    mean_fc_gene <- read_group_tracking %>% rename(gene=tracking_id) %>% 
      dplyr::group_by(condition, gene) %>% 
      dplyr::summarise(mean = mean(FPKM)) %>% inner_join(gene_list)
    
    
    if(optionABC == 'A'){
      
      med_fc_table <- med_fc_gene %>% spread(key = "condition", value = "med") %>% 
        mutate(med_d1r_ip_ov_n = Drd1_IP/Drd1_N, med_d1rip_ov_vg2ip = Drd1_IP/Vglut2_IP)
      
      mean_fc_table <- mean_fc_gene %>% spread(key = "condition", value = "mean") %>% 
        mutate(mean_d1r_ip_ov_n = Drd1_IP/Drd1_N, mean_d1rip_ov_vg2ip = Drd1_IP/Vglut2_IP)
      
    } else if(optionABC == 'B') {
      
      med_fc_table <- med_fc_gene %>% spread(key = "condition", value = "med") %>% 
        mutate(med_vg2_ip_ov_n = Vglut2_IP/Vglut2_N, med_d1r_ip_ov_n = Drd1_IP/Drd1_N)
      
      mean_fc_table <- mean_fc_gene %>% spread(key = "condition", value = "mean") %>% 
        mutate(mean_vg2_ip_ov_n = Vglut2_IP/Vglut2_N, mean_d1r_ip_ov_n = Drd1_IP/Drd1_N)
      
    } else if(optionABC == 'C') {
      
      med_fc_table <- med_fc_gene %>% spread(key = "condition", value = "med") %>% 
        mutate(med_vg2_ip_ov_n = Vglut2_IP/Vglut2_N, med_vg2ip_ov_d1rip = Vglut2_IP/Drd1_IP)
      
      mean_fc_table <- mean_fc_gene %>% spread(key = "condition", value = "mean") %>% 
        mutate(mean_vg2_ip_ov_n = Vglut2_IP/Vglut2_N, mean_vg2ip_ov_d1rip = Vglut2_IP/Drd1_IP)
      
    } else {
      optionABC = 'A'
    }
    
    fold_change_table <- inner_join(med_fc_table, mean_fc_table, by = 'gene')
  
  return(fold_change_table[,c(1,6,7,12,13)])

}




#GOs to highlights


#turquoise

turq_go_2_hl <- c('GO:0008094','GO:0043565','GO:0000976','GO:0031981',
                  'GO:0006396','GO:0006259','GO:0000226','WP:WP310')

# top10 2x enriched and 2x DE genes for highlighted terms
top10s_turq_d1r <- split_then_top(turq_top_table, turq_go_2_hl, top_d1r_comb, 10) %>% 
  t() %>% as.data.frame()

top10s_turq_d1r

# top10 2x enriched and 2x DE genes for all terms from gprofiler result
top10s_turq_d1r_all <- split_then_top(turq_top_table,turq_top_table$term_id, top_d1r_comb, 10) %>% 
  t() %>% as.data.frame()

# single list of genes from highlighted terms
top10s_turq_d1r_genes <- split_then_top_list(turq_top_table, turq_go_2_hl, top_d1r_comb, 10)

#median and mean fold-change for highlighted genes / option A for D1R+ enriched and DE
turq_fc <- fold_change_calc(top10s_turq_d1r_genes, rgt, 'A')
turq_fc

# violin plots for genes from highlighted terms
violin_maker(top10s_turq_d1r_genes$gene,rgt)



#midnight blue

mnblue_go_2_hl <- c('GO:0014069','GO:0098978','GO:0043005','GO:0098916',
                    'GO:0048167','GO:0048812','GO:0008306','KEGG:04080')

# top10 2x enriched and 2x DE genes for highlighted terms
top10s_mnblue_d1r <- split_then_top(mnblued1rtop_table, mnblue_go_2_hl, top_d1r_comb, 10) %>% t() %>% 
  as.data.frame()

top10s_mnblue_d1r

# top10 2x enriched and 2x DE genes for all terms from gprofiler result
top10s_mnblue_d1r_all <- split_then_top(mnblued1rtop_table ,mnblued1rtop_table$term_id, top_d1r_comb, 10) %>% t() %>% 
  as.data.frame()

# single list of genes from highlighted terms
top10s_mnblue_d1r_genes <- split_then_top_list(mnblued1rtop_table, mnblue_go_2_hl, top_d1r_comb, 10)

#median and mean fold-change for highlighted genes / option A for D1R+ enriched and DE
mnblue_fc <- fold_change_calc(top10s_mnblue_d1r_genes, rgt, 'A')
mnblue_fc

# violin plots for genes from highlighted terms
violin_maker(top10s_mnblue_d1r_genes$gene,rgt)



### d1r tan 

tan_go_2_hl <- c('GO:0043005','GO:0098978','GO:0007611',
                 'GO:0098916','KEGG:04080','KEGG:04725',
                 'REAC:R-MMU-418594','REAC:R-MMU-416476')

# top10 2x enriched and 2x DE genes for highlighted terms
top10s_tan_d1r <- split_then_top(tan_top_table, tan_go_2_hl, top_d1r_comb, 10) %>% t() %>% as.data.frame()

top10s_tan_d1r

# top10 2x enriched and 2x DE genes for all terms from gprofiler result
top10s_tan_d1r_all <- split_then_top(tan_top_table,tan_top_table$term_id, top_d1r_comb, 10) %>% t() %>% 
  as.data.frame()


# single list of genes from highlighted terms
top10s_tan_d1r_genes <- split_then_top_list(tan_top_table, tan_go_2_hl, top_d1r_comb, 10)

#median and mean fold-change for highlighted genes / option A for D1R+ enriched and DE
tan_fc <- fold_change_calc(top10s_tan_d1r_genes, rgt, 'A')
tan_fc

# violin plots for genes from highlighted terms
violin_maker(top10s_tan_d1r_genes$gene,rgt)



#enriched in both d1r and vglut2
top_both_enr <- inner_join(top_d1r_n_v_ip[, c(3,16)], top_vg2_n_v_ip[, c(3,16)], by = 'gene') %>% 
  mutate(rank_comb = rank.x + rank.y)

# tan module enriched in both

# go to highlight 
tan_both_go_2_hl <- c('GO:0005102','GO:0005515','GO:0043005','GO:0005886',
                      'GO:0071944','GO:0099537','GO:0048666','GO:0022008')

# top10 2x enriched and 2x DE genes for highlighted terms
top10s_tan_both_enr <- split_then_top(tan_top_both_table, tan_both_go_2_hl, top_both_enr, 10) %>% 
  t() %>% as.data.frame()

top10s_tan_both_enr

# top10 2x enriched and 2x DE genes for all terms from gprofiler result
top10s_tan_both_enr_all <- split_then_top(tan_top_both_table, tan_top_both_table$term_id, top_both_enr, 10) %>% 
  t() %>% as.data.frame()

# single list of genes from highlighted terms
top10s_tan_both_genes <- split_then_top_list(tan_top_both_table, tan_both_go_2_hl, top_both_enr, 10)

#median and mean fold-change for highlighted genes / option B for enriched in both
tan_both_fc <- fold_change_calc(top10s_tan_both_genes, rgt, 'B')
tan_both_fc

# violin plots for genes from highlighted terms
violin_maker(top10s_tan_both_genes$gene,rgt)



# vg2 black

black_go_2_hl <- c('GO:0045202','GO:0043005','GO:0006811','GO:0050808',
                   'KEGG:05012','KEGG:05014','KEGG:05020','KEGG:05010')

# top10 2x enriched and 2x DE genes for highlighted terms
top10s_black_vg2 <- split_then_top(blackvg2top_table, black_go_2_hl, top_vg2_comb, 10) %>% 
  t() %>% as.data.frame()

top10s_black_vg2

# top10 2x enriched and 2x DE genes for all terms from gprofiler result
top10s_black_vg2_all <- split_then_top(blackvg2top_table,blackvg2top_table$term_id, top_vg2_comb, 10) %>% t() %>% 
  as.data.frame()

# single list of genes from highlighted terms
top10s_black_vg2_genes <- split_then_top_list(blackvg2top_table, black_go_2_hl, top_vg2_comb, 10)

#median and mean fold-change for highlighted genes / option C for VGLUT2+ enriched and DE
black_vg2_fc <- fold_change_calc(top10s_black_vg2_genes, rgt, 'C')
black_vg2_fc

# violin plots for genes from highlighted terms
violin_maker(top10s_black_vg2_genes$gene,rgt)



#rainbow 

# vg2 rainbow

rainbow_go_2_hl <- c('GO:0005509','GO:0043005','GO:0045202','GO:0031982',
                     'GO:0070887','GO:0006811','GO:0034097','KEGG:04721')

# top10 2x enriched and 2x DE genes for highlighted terms
top10s_rainbow_vg2 <- split_then_top(comb_lylgs_vg2top_table, rainbow_go_2_hl, top_vg2_comb, 10) %>% 
  t() %>% as.data.frame()

top10s_rainbow_vg2

# top10 2x enriched and 2x DE genes for all terms from gprofiler result
top10s_rainbow_vg2_all <- split_then_top(comb_lylgs_vg2top_table,comb_lylgs_vg2top_table$term_id, top_vg2_comb, 10) %>% t() %>% 
  as.data.frame()

# single list of genes from highlighted terms
top10s_rainbow_vg2_genes <- split_then_top_list(comb_lylgs_vg2top_table, rainbow_go_2_hl, top_vg2_comb, 10)

#median and mean fold-change for highlighted genes / option C for VGLUT2+ enriched and DE
rainbow_vg2_fc <- fold_change_calc(top10s_rainbow_vg2_genes, rgt, 'C')
rainbow_vg2_fc

# violin plots for genes from highlighted terms
violin_maker(top10s_rainbow_vg2_genes$gene,rgt)
