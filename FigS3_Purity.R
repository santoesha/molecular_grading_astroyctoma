# load libs ----
library(tidyverse)
library(Hmisc)
library(corrplot) 

source("scripts/clean/R/theme_cellpress.R")

# Correlation methods ------

# ** CATNON -----
meta_catnon <- readRDS("data/clean/catnon/rds/metadata.catnon.Rds") %>% 
               dplyr::filter(RNA_Seq == T) %>% 
               dplyr::filter(predictBrain_12.8_cal_class %in% c("A_IDH_LG", "A_IDH_HG"))
counts.catnon <- readRDS("data/clean/catnon/rds/expression.catnon.vst.Rds") %>% 
                 dplyr::select(meta_catnon$GS_ID)

purities_cat <- readRDS("data/clean/catnon/rds/catnon_purities.Rds") %>% 
            dplyr::filter(BatchID %in% meta_catnon$BatchID)
strom_sig_cat <- readRDS("output/clean/catnon/rds/mckenzie_sign_pca.Rds")

meta_catnon <- meta_catnon %>% 
               dplyr::left_join(., purities_cat, by = c("BatchID" = "BatchID")) %>% 
               dplyr::filter(! is.na(vaf_idh_x2)) %>% 
               dplyr::left_join(., strom_sig_cat) %>%
               dplyr::rename(mic_exp_rna = lts.mic) %>%
               dplyr::rename(oli_exp_rna = lts.oli) %>% 
               dplyr::rename(neu_exp_rna = lts.neu) %>% 
               dplyr::rename(ABSOLUTE_meth = abs_purity) %>% 
               dplyr::rename(ESTIMATE_meth = est_purity) %>% 
               dplyr::rename(InfiniumPurify_meth = getpurity_lgg_meth) %>% 
               dplyr::rename(VAF_IDH_DNA = vaf_idh_x2) 

cor.mat_catnon <- rcorr(as.matrix(meta_catnon %>% dplyr::select(VAF_IDH_DNA, InfiniumPurify_meth, neu_exp_rna, mic_exp_rna, oli_exp_rna),type="spearman"))$r

cor.mat_catnon <- cor.mat_catnon[,c(1:2)]

# pdf(file = "output/figures/paper/sf1/clean/corplot_tp_methods_catnon.pdf", width = 10, height = 10)
corrplot::corrplot(cor.mat_catnon, method="square", diag = F,
         title = paste0("CATNON (n=", nrow(meta_catnon), ")"),
         mar=c(0,0,1,0),
         col = COL2('BrBG'), tl.cex = 1.2, tl.col = "black",
         type = "lower")
# dev.off()

# ** TCGA ----

meta_tcga <- readRDS("data/clean/tcga/rds/metadata_tcga_idhmutnoncodel.Rds") %>% 
                     dplyr::filter(RNA_Seq == T) %>% 
                    dplyr::filter(predictBrain_12.8_cal_class %in% c("A_IDH_LG", "A_IDH_HG"))
counts.tcga <- readRDS("data/clean/tcga/rds/expression.tcga.astr.vst.Rds") %>% 
                       dplyr::select(meta_tcga$Case)

purities_tcga <- readRDS("data/clean/tcga/rds/tcga_purities_updat.Rds") %>% 
                  dplyr::filter(BatchID %in% meta_tcga$BatchID) %>% 
                  dplyr::select(-Case)
strom_sig_tcga <- readRDS("output/clean/tcga/rds/mckenzie_sign_pca.Rds")

meta_tcga <- meta_tcga %>% 
              dplyr::left_join(., purities_tcga, by = c("BatchID" = "BatchID")) %>% 
              dplyr::left_join(., strom_sig_tcga) %>%
              dplyr::rename(mic_exp_rna = lts.mic) %>%
              dplyr::rename(oli_exp_rna = lts.oli) %>% 
              dplyr::rename(neu_exp_rna = lts.neu) %>% 
              dplyr::rename(ABSOLUTE_meth = abs) %>% 
              dplyr::rename(ESTIMATE_meth = est) %>% 
              dplyr::rename(InfiniumPurify_meth = getpurity_lgg_meth) %>% 
              dplyr::rename(VAF_IDH_DNA = vaf_idh_x2) 

cor.mat_tcga <- rcorr(as.matrix(meta_tcga %>% dplyr::select(VAF_IDH_DNA, InfiniumPurify_meth, neu_exp_rna, mic_exp_rna, oli_exp_rna),type="spearman"))$r

cor.mat_tcga <- cor.mat_tcga[,c(1:2)]

#pdf(file = "output/figures/paper/sf1/clean/corplot_tp_methods_tcga.pdf", width = 10, height = 10)
corrplot::corrplot(cor.mat_tcga, method="square", diag = F,
                   title = paste0("TCGA (n=", nrow(meta_tcga), ")"),
                   mar=c(0,0,1,0),
                   col = COL2('BrBG'), tl.cex = 1.2, tl.col = "black",
                   type = "lower")
#dev.off()

# Correlation LFC and methods -----

# ** TCGA ----
getpurity.tcga <- data.frame(cor.t.meth.infiniumpurity = 
                              apply(counts.tcga,1, function(vec) {return( cor.test(meta_tcga$InfiniumPurify_meth, 
                                                                                   as.numeric(vec), method = "pearson")$estimate) })) %>% 
                dplyr::mutate(gene_name = gsub("^.+_","", rownames(.)))

res_ord_tcga <- readRDS("output/clean/tcga/rds/resord_tcga_heidi_num_lncrna.Rds") %>% 
                dplyr::left_join(., getpurity.tcga, by = c("gene_name" = "gene_name")) %>% 
                dplyr::select(log2FoldChange, cor.t.meth.infiniumpurity) %>% 
                tidyr::pivot_longer(!log2FoldChange) %>% 
                dplyr::rename(cor = value) %>% 
                dplyr::rename(method = name) %>% 
                dplyr::filter(! is.na(cor)) %>% 
                dplyr::mutate(Study = "TCGA")

# ** CATNON ----
getpurity.cat <- data.frame(cor.t.meth.infiniumpurity = 
                 apply(counts.catnon,1, function(vec) {return( cor.test(meta_catnon$InfiniumPurify_meth, 
                                                                                     as.numeric(vec), method = "pearson")$estimate) })) %>% 
                 dplyr::mutate(gene_name = gsub("^.+_","", rownames(.)))

res_ord_cat <- readRDS("output/clean/catnon/rds/resord_catnon_heidi_num_lncrna.Rds") %>% 
               dplyr::left_join(., getpurity.cat, by = c("gene_name" = "gene_name")) %>% 
               dplyr::select(log2FoldChange, cor.t.meth.infiniumpurity) %>% 
               tidyr::pivot_longer(!log2FoldChange) %>% 
               dplyr::rename(cor = value) %>% 
               dplyr::rename(method = name) %>% 
               dplyr::filter(! is.na(cor))  %>% 
               dplyr::mutate(Study = "CATNON")
            
plt <- rbind(res_ord_tcga, res_ord_cat)
       
ggplot(plt, aes(x=,log2FoldChange, y=cor)) +
       facet_wrap(~ Study, ncol = 2) +
       xlim(-1,1) +
       theme_classic() +
       ggpubr::stat_cor(method="spearman", aes(label = ..r.label..), size = theme_cellpress_size) +
       #geom_smooth(method = "lm", se = FALSE, formula = y~x,col="#ff2929") +
       labs(x = "log2FoldChange LGC (RNA)", y = "Correlation gene with tumor purity") +
       geom_point(data = plt %>% dplyr::filter(abs(log2FoldChange) > 0.01),
                  pch=21, fill=alpha("white",0.6), cex=2, col = "grey10") +
       theme_cellpress

ggsave("output/figures/paper/sf1/clean/cor_lfc_purity.pdf", width = 10, height = 5)

ggsave("output/figures/paper/rev1/sf2/cor_lfc_purity.pdf", width = 11.2 * 0.95, height = (8.5 * 0.95) * (1/2))

# CDKN2AB -----
cdkn2a_pur_tcga <- meta_tcga %>% 
                   dplyr::select(CDKN2AB_del, InfiniumPurify_meth) %>% 
                   dplyr::mutate(Study = "TCGA")
cdkn2a_pur_cat <- meta_catnon %>% 
                  dplyr::select(CDKN2AB_del, InfiniumPurify_meth) %>% 
                  dplyr::mutate(Study = "CATNON")

plt_cdkn2a <- rbind(cdkn2a_pur_tcga, cdkn2a_pur_cat)

ggplot(plt_cdkn2a, aes(x = CDKN2AB_del, y = InfiniumPurify_meth, 
                col = CDKN2AB_del,fill = CDKN2AB_del)) +
      facet_wrap(~ Study) +
      ggbeeswarm::geom_quasirandom(dodge.width = .8, alpha = .3) +
      geom_boxplot(width = .4, alpha = .33, position = position_dodge(.8), outlier.shape = NA) +
      stat_compare_means(method = "wilcox.test", label = "p.format", size = 6) +
      scale_color_manual(values = c("grey", "black")) +
      scale_fill_manual(values = c("grey", "black")) +
      labs(x = "", y = "Estimated purity (methylation)") +
      stat_n_text(size = 6) +
      theme_cellpress +
      theme(strip.background=element_rect(colour="black",
                                          fill="forestgreen"))
# ggsave("output/figures/paper/sf1/clean/cdkn2ab_purity_boxplot.pdf", width = 10, height = 10)
