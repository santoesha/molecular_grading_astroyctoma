# Load packages ----
library(tidyverse)
library(patchwork)
library(cowplot)
library(EnvStats)
library(ggpubr)
library(ggdensity)

# Source ----

source("scripts/clean/R/job_gg_theme.R")
source("scripts/clean/R/palette.R")
source("scripts/clean/R/sel_genes_dge.R")
source("scripts/clean/R/theme_cellpress.R")

# CATNON ----

meta.catnon <- readRDS("data/clean/catnon/rds/metadata.catnon.Rds") %>% 
               dplyr::filter(RNA_Seq == T) %>% 
               dplyr::filter(predictBrain_12.8_cal_class %in% c("A_IDH_LG", "A_IDH_HG"))

counts.catnon <- readRDS("data/clean/catnon/rds/expression.catnon.vst.Rds") %>% 
                 dplyr::select(meta.catnon$GS_ID) %>% 
                 dplyr::filter(rownames(.) %in% sel_genes)

plt.pca.catnon <- counts.catnon %>%
                  dplyr::mutate(mad =  apply( as.matrix(.), 1, stats::mad) ) %>%  # use MAD instead of SD to find most variable genes
                  dplyr::arrange(-mad) %>% # arrange in asc? order
                  dplyr::slice_head(n = 2000) %>%  # pick top 1000
                  dplyr::mutate(mad = NULL) %>% # remove the sd to obtain original vst matrix
                  t() %>% # transpose, to PCA the genes rather than the patients
                  prcomp(scale = T) 

plt.catnon.scores <-  plt.pca.catnon %>% # PCA
                        purrr::pluck('x') %>%  # take coordinates
                        as.data.frame(stringsAsFactor=F) %>% # transform back from matrix to data.frame 
                        dplyr::select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% 
                        tibble::rownames_to_column('GS_ID') %>% 
                        dplyr::left_join(meta.catnon , by=c('GS_ID'='GS_ID')) 

plt.catnon.loadings <-  plt.pca.catnon %>% 
                      purrr::pluck('rotation') %>%  # take coordinates
                      as.data.frame(stringsAsFactor=F) %>% # transform back from matrix to data.frame 
                      dplyr::select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% 
                      tibble::rownames_to_column('gene') %>% 
                      dplyr::mutate(gene = gsub("^.+_","", gene)) 

plt.catnon.loadings %>% dplyr::slice_min(PC2, n=50) %>% dplyr::pull(gene)

# TCGA ----

meta.tcga <- readRDS("data/clean/tcga/rds/metadata_tcga_idhmutnoncodel.Rds") %>% 
               dplyr::filter(RNA_Seq == T) %>% 
               dplyr::filter(predictBrain_12.8_cal_class %in% c("A_IDH_LG", "A_IDH_HG"))

counts.tcga <- readRDS("data/clean/tcga/rds/expression.tcga.astr.vst.Rds") %>% 
               dplyr::select(meta.tcga$Case) %>% 
               dplyr::filter(rownames(.) %in% sel_genes)

plt.pca.tcga <-   counts.tcga %>%
                  dplyr::mutate(mad =  apply( as.matrix(.), 1, stats::mad) ) %>%  # use MAD instead of SD to find most variable genes
                  dplyr::arrange(-mad) %>% # arrange in asc? order
                  dplyr::slice_head(n = 2000) %>%  # pick top 1000
                  dplyr::mutate(mad = NULL) %>% # remove the sd to obtain original vst matrix
                  t() %>% # transpose, to PCA the genes rather than the patients
                  prcomp(scale = T)

summary(plt.pca.tcga)

plt.tcga.scores <-  plt.pca.tcga %>% # PCA
                        purrr::pluck('x') %>%  # take coordinates
                        as.data.frame(stringsAsFactor=F) %>% # transform back from matrix to data.frame 
                        dplyr::select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% 
                        tibble::rownames_to_column('Case') %>% 
                        dplyr::left_join(meta.tcga , by=c('Case'='Case')) 


plt.tcga.loadings <-  plt.pca.tcga %>% 
                          purrr::pluck('rotation') %>%  # take coordinates
                          as.data.frame(stringsAsFactor=F) %>% # transform back from matrix to data.frame 
                          dplyr::select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% 
                          tibble::rownames_to_column('gene') %>% 
                          dplyr::mutate(gene = gsub("^.+_","", gene)) 

plt.tcga.loadings %>% dplyr::slice_max(PC2, n=50) %>% dplyr::pull(gene)

# GLASS-R ----
meta.glass.r <- readRDS("data/clean/glass/rds/metadata.glass.Rds") %>% 
                dplyr::filter(RNA_Seq == T) %>% 
                dplyr::filter(predictBrain_12.8_cal_class %in% c("A_IDH_LG", "A_IDH_HG")) %>% 
                dplyr::filter(Resectie != 1)
          
counts.glass.r <- readRDS("data/clean/glass/rds/expression.glass.exon.vst.Rds") %>% 
                dplyr::select(meta.glass.r$GS_ID) %>% 
                dplyr::filter(rownames(.) %in% sel_genes)

plt.pca.glass.r <-   counts.glass.r %>%
                     dplyr::mutate(mad =  apply( as.matrix(.), 1, stats::mad) ) %>%  # use MAD instead of SD to find most variable genes
                     dplyr::arrange(-mad) %>% # arrange in asc? order
                     dplyr::slice_head(n = 2000) %>%  # pick top 1000
                     dplyr::mutate(mad = NULL) %>% # remove the sd to obtain original vst matrix
                     t() %>% # transpose, to PCA the genes rather than the patients
                     prcomp(scale = T)

summary(plt.pca.glass)

plt.glass.r.scores <-  plt.pca.glass.r %>% # PCA
                          purrr::pluck('x') %>%  # take coordinates
                          as.data.frame(stringsAsFactor=F) %>% # transform back from matrix to data.frame 
                          dplyr::select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% 
                          dplyr::mutate(PC2 = PC2 * -1) %>% 
                          tibble::rownames_to_column('GS_ID') %>% 
                          dplyr::left_join(meta.glass.r , by=c('GS_ID'='GS_ID')) 


plt.glass.r.loadings <-  plt.pca.glass.r %>% 
                         purrr::pluck('rotation') %>%  # take coordinates
                             as.data.frame(stringsAsFactor=F) %>% # transform back from matrix to data.frame 
                             dplyr::select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% 
                             tibble::rownames_to_column('gene') %>% 
                             dplyr::mutate(gene = gsub("^.+_","", gene)) 


plt.glass.r.loadings %>% dplyr::slice_max(PC2, n=50) %>% dplyr::pull(gene)

# GLASS-P ----
meta.glass.p <- readRDS("data/clean/glass/rds/metadata.glass.Rds") %>% 
                dplyr::filter(RNA_Seq == T) %>% 
                dplyr::filter(predictBrain_12.8_cal_class %in% c("A_IDH_LG", "A_IDH_HG")) %>% 
                dplyr::filter(Resectie == 1)

counts.glass.p <- readRDS("data/clean/glass/rds/expression.glass.exon.vst.Rds") %>% 
                  dplyr::select(meta.glass.p$GS_ID) %>% 
                  dplyr::filter(rownames(.) %in% sel_genes)

plt.pca.glass.p <-   counts.glass.p %>%
                      dplyr::mutate(mad =  apply( as.matrix(.), 1, stats::mad) ) %>%  # use MAD instead of SD to find most variable genes
                      dplyr::arrange(-mad) %>% # arrange in asc? order
                      dplyr::slice_head(n = 2000) %>%  # pick top 1000
                      dplyr::mutate(mad = NULL) %>% # remove the sd to obtain original vst matrix
                      t() %>% # transpose, to PCA the genes rather than the patients
                      prcomp(scale = T)

summary(plt.pca.glass.p)

plt.glass.p.scores <-  plt.pca.glass.p %>% # PCA
                       purrr::pluck('x') %>%  # take coordinates
                       as.data.frame(stringsAsFactor=F) %>% # transform back from matrix to data.frame 
                       dplyr::select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% 
                       dplyr::mutate(PC2 = PC2 * -1) %>% 
                       tibble::rownames_to_column('GS_ID') %>% 
                       dplyr::left_join(meta.glass.p , by=c('GS_ID'='GS_ID'))


plt.glass.p.loadings <-  plt.pca.glass.p %>% 
                         purrr::pluck('rotation') %>%  # take coordinates
                         as.data.frame(stringsAsFactor=F) %>% # transform back from matrix to data.frame 
                         dplyr::select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10) %>% 
                         tibble::rownames_to_column('gene') %>% 
                         dplyr::mutate(gene = gsub("^.+_","", gene)) 
 
rm(sel_genes)

# plt ----

plt.tcga <- plt.tcga.scores %>% dplyr::select(PC1, PC2, PC3, PC4, CDKN2AB_del, log_AIDH_cal_ratio,
                                              CDK4_amp, PDGFRA_amp, RB1_del) %>% 
            dplyr::mutate(Study = paste0("TCGA (n=", nrow(meta.tcga), ")")) 

plt.catnon <- plt.catnon.scores %>% dplyr::select(PC1, PC2, PC3, PC4, CDKN2AB_del, log_AIDH_cal_ratio,
                                                  CDK4_amp, PDGFRA_amp, RB1_del) %>% 
              dplyr::mutate(Study = paste0("CATNON (n=", nrow(meta.catnon), ")")) %>% dplyr::mutate(PC2 = PC2*-1)

plt.glass.p <- plt.glass.p.scores %>% dplyr::select(PC1, PC2, PC3, PC4, CDKN2AB_del, log_AIDH_cal_ratio,
                                                    CDK4_amp, PDGFRA_amp, RB1_del) %>% 
                dplyr::mutate(Study = paste0("GLASS-NL-P (n=", nrow(meta.glass.p), ")")) %>% dplyr::mutate(PC2 = PC2*-1)

plt.glass.r <- plt.glass.r.scores %>% dplyr::select(PC1, PC2, PC3, PC4, CDKN2AB_del, log_AIDH_cal_ratio,
                                                    CDK4_amp, PDGFRA_amp, RB1_del) %>%
                  dplyr::mutate(Study = paste0("GLASS-NL-R (n=", nrow(meta.glass.r), ")"))

plt <- rbind(plt.tcga, plt.catnon, plt.glass.p, plt.glass.r) %>% 
       dplyr::mutate(prog_alter = ifelse(CDKN2AB_del == T | 
                                         CDK4_amp == T |
                                         PDGFRA_amp == T | 
                                         RB1_del == T, T, F)) %>% 
       dplyr::mutate(Study = factor(Study, levels = c("CATNON (n=132)", "TCGA (n=235)", "GLASS-NL-P (n=62)", "GLASS-NL-R (n=87)")))

ggplot(plt, aes(x = PC2, y = log_AIDH_cal_ratio)) +
      facet_wrap(~ Study) +
      geom_point(aes(col = CDKN2AB_del), alpha =0.7) +
      theme_cellpress +
      scale_color_manual(values = c("grey", "#FF0000")) +
      geom_smooth(method = "lm", se=FALSE, col = "grey60") +
      ggpubr::stat_cor(label.y = 21, label.x = -15, col = "black", method="spearman", size = theme_cellpress_size, cor.coef.name = "rho") +
      labs(x = "Principal Component 2 (RNA)", y = "Continuous Grading Coefficient (DNAm)")
    
ggsave("output/figures/paper/rev1/f2/pca_datasets_rna.pdf", width = (11.2 * 0.95) * 0.30, height = (8.5 * 0.95) * 0.40)

# S1: Unsupervised PCA CATNON and TCGA -----


## ** CATNON -----
p1 <- ggplot(plt.pca.catnon, aes(x = PC1, y=PC2, fill = log_AIDH_cal_ratio)) +
      geom_point(size = 4, pch = 21, alpha = 2, stroke = NA) +
      scale_fill_gradient2(midpoint = 0, low = "yellow", mid = "rosybrown",high = "black", space = "Lab" ) +
      theme_classic() +
      labs(x = "CATNON: PC1 (RNA)", y = "PC2 (RNA)", title = paste0("CATNON (n=", 
                                                                    nrow(plt.pca.catnon),
                                                                    ")")) +
      scale_y_continuous(breaks = seq(-20, 20, 20)) +
      theme_cellpress +
      theme(legend.position = "none",
            panel.border = element_rect(colour = "black", fill=NA, size=0.5)) #+
    #labs(title = paste0("CATNON:", " n=", nrow(plt.pca.catnon)))

ydens <- axis_canvas(p1, axis = "y") + 
         geom_density(data = plt.pca.catnon, aes(y = PC2, fill = dkfz_diagnosis, colour = dkfz_diagnosis), alpha = 0.3) +
         scale_fill_manual(values = heidi_color) +
         scale_colour_manual(values = heidi_color) 

p1 <- p1 %>% 
      insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
      cowplot::ggdraw()

rm(ydens)

p2 <- ggplot(plt.pca.catnon, aes(y = PC2, x = log_AIDH_cal_ratio)) +
      geom_point(alpha = 0.8, size = 6, pch = 21, aes(fill =  CDKN2AB_del)) +
      scale_fill_manual(values = c("grey", "#FF0000")) +
      theme_classic() +
      geom_smooth(method = "lm", se=FALSE, col = "black") +
      labs(x = "CATNON: Continuous Grading Coefficient (CGC, methylation)", y = "PC2 (RNA)") +
      ggpubr::stat_cor(label.y = 21, label.x = -15, col = "black", method="spearman", size = theme_cellpress_size, cor.coef.name = "rho") +
      theme(legend.position = "none",
            panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
      scale_y_continuous(breaks = seq(-20, 20, 20))

# ** TCGA -----

p3 <- ggplot(plt.pca.tcga, aes(x = PC1, y=PC2, fill = log_AIDH_cal_ratio)) +
  geom_point(size = 5, pch = 21, alpha = 2, stroke = NA) +
  scale_fill_gradient2(midpoint = 0, low = "yellow", mid = "rosybrown",high = "black", space = "Lab", name = "LGC" ) +
  labs(x = "TCGA: PC1 (RNA)", y = "", title = paste0("TCGA (n=", 
                                                     nrow(plt.pca.tcga),
                                                     ")")) +
  theme_classic() +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) #+
#labs(title = paste0("TCGA:", " n=", nrow(plt.pca.tcga))) 

ydens <- axis_canvas(p3, axis = "y") + 
  geom_density(data = plt.pca.tcga, aes(y = PC2, fill = dkfz_diagnosis, colour = dkfz_diagnosis), alpha = 0.3) +
  scale_fill_manual(values = heidi_color) +
  scale_colour_manual(values = heidi_color) 

p3 <- p3 %>% 
  insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
  cowplot::ggdraw()

rm(ydens)

p4 <- ggplot(plt.pca.tcga, aes(y = PC2, x = log_AIDH_cal_ratio)) +
  geom_point(alpha = 0.8, size = 6, pch = 21, aes(fill =  CDKN2AB_del)) +
  scale_fill_manual(values = c("grey", "#FF0000")) +
  geom_smooth(method = "lm", se=FALSE, col = "black") +
  theme_classic() +
  labs(x = "TCGA: Linear Grading Coefficient (LGC, methylation)", y = "PC2 (RNA)") +
  ggpubr::stat_cor(label.y = 30, label.x = -15, col = "black", method="spearman", size = 10, cor.coef.name = "rho") +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_y_continuous(breaks = seq(-20, 20, 20))


# ** GLASS -----

p5 <- ggplot(plt.pca.glass, aes(x = PC1, y=PC2, fill = log_AIDH_cal_ratio)) +
      geom_point(size = 5, pch = 21, alpha = 2, stroke = NA) +
      scale_fill_gradient2(midpoint = 0, low = "yellow", mid = "rosybrown",high = "black", space = "Lab", name = "LGC" ) +
      theme_classic() +
      labs(x = "GLASS-NL-R: PC1 (RNA)", y = "", title = paste0("GLASS-NL-R (n=", 
                                                               nrow(plt.pca.glass),
                                                               ")")) +
      theme(legend.position = "none",
            panel.border = element_rect(colour = "black", fill=NA, size=0.5)) #+
#labs(title = paste0("TCGA:", " n=", nrow(plt.pca.glass))) 

ydens <- axis_canvas(p5, axis = "y") + 
  geom_density(data = plt.pca.glass, aes(y = PC2, fill = dkfz_diagnosis, colour = dkfz_diagnosis), alpha = 0.3) +
  scale_fill_manual(values = heidi_color) +
  scale_colour_manual(values = heidi_color) 

p5 <- p5 %>% 
  insert_yaxis_grob(ydens, grid::unit(1, "in"), position = "right") %>%
  cowplot::ggdraw()

rm(ydens)

p6 <- ggplot(plt.pca.glass, aes(y = PC2, x = log_AIDH_cal_ratio)) +
  geom_point(alpha = 0.8, size = 6, pch = 21, aes(fill =  CDKN2AB_del)) +
  scale_fill_manual(values = c("grey", "#FF0000")) +
  geom_smooth(method = "lm", se=FALSE, col = "black") +
  theme_classic() +
  labs(x = "GLASS-NL-R: Linear Grading Coefficient (LGC, methylation)", y = "PC2 (RNA)") +
  ggpubr::stat_cor(label.y = 30, label.x = -15, col = "black", method="spearman", size = 10, cor.coef.name = "rho") +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_y_continuous(breaks = seq(-20, 20, 20))

p5+p6

p2 + p4 + p6

remove_y <- theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank()
)

p_1 <- list(p2,p4 + remove_y,p6 + remove_y)

a <- wrap_plots(p_1, nrow = 1) + plot_layout(guides = "collect")

p_2 <- list(
  p1,
  p3 + remove_y,
  p5 + remove_y)

b <- wrap_plots(p_2, nrow = 1) + plot_layout(guides = "collect")

ggpubr::ggarrange(b,a, nrow = 2, align = 'h')

ggsave("output/figures/paper/fig1/pca_all_datasets2.pdf", width = 30, height = 10)

p1

p1 + p3 + p5

p2+p4

p1/p3

p1
# ggsave("output/figures/paper/sf1/unsupervised_pca_catnon.pdf", width = 6, height = 3)

p2
# ggsave("output/figures/paper/sf1/unsupervised_pc2vslgc_catnon.pdf", width = 6, height = 3)


p3
# ggsave("output/figures/paper/sf1/unsupervised_pca_tcga.pdf", width = 8, height = 5)

p4
# ggsave("output/figures/paper/sf1/unsupervised_pc2vslgc_tcga.pdf", width = 6, height = 3)
# 

# WHO vs LGC -----
pc2_cat <- plt.pca.catnon %>% 
           dplyr::rename(Sample_ID = GS_ID) %>%
           dplyr::mutate(Study = "CATNON") %>% 
           dplyr::select(PC1, PC2, Sample_ID, CDKN2AB_del, log_AIDH_cal_ratio, who_2021, Study) 
           
pc2_tcga <- plt.pca.tcga %>% 
            dplyr::rename(Sample_ID = Case) %>% 
            dplyr::mutate(Study = "TCGA") %>% 
            dplyr::mutate(who_2021 = NA) %>% 
            dplyr::select(PC1, PC2, CDKN2AB_del,Sample_ID, log_AIDH_cal_ratio, who_2021, Study) 

pc2_glass <- plt.pca.glass %>% 
             dplyr::rename(Sample_ID = GS_ID) %>% 
             dplyr::mutate(Study = "GLASS-NL-R") %>% 
             dplyr::mutate(who_2021 = ifelse(WHO_Classification2021 == "Astrocytoma, IDH-mutant, WHO grade 4",
                                             "WHO IDHmut grade 4", "WHO IDHmut grade 2/3")) %>% 
             dplyr::select(PC1, PC2, CDKN2AB_del,Sample_ID, log_AIDH_cal_ratio, who_2021, Study) 
          
pc2 <- rbind(pc2_cat, pc2_tcga, pc2_glass) %>% 
       dplyr::mutate(Study = factor(Study, levels = c("CATNON", "TCGA", "GLASS-NL-R")))

# saveRDS(pc2, "output/clean/combined/rds/unsupervised_pc2_rna.Rds")

boxplt1 <- ggplot(pc2, aes(x = Study, y = log_AIDH_cal_ratio,
                                col = CDKN2AB_del,fill = CDKN2AB_del)) +
          ggbeeswarm::geom_quasirandom(dodge.width = .8, alpha = .3, aes(col = CDKN2AB_del)) +
          geom_boxplot(width = .4, alpha = .33, position = position_dodge(.8), outlier.shape = NA) +
          stat_compare_means(method = "wilcox.test", label = "p.format", size = 6) +
          scale_color_manual(values = c("grey", "black")) +
          scale_fill_manual(values = c("grey", "black")) +
          labs(x = "", y = "Linear Grading Coefficient (LGC, methylation)") +
          stat_n_text(size = 6) +
          theme_classic() +
          theme(axis.text = element_text(size = 15),
              axis.title.y = element_text(size = 15))

boxplt2 <- ggplot(pc2, aes(x = Study, y = PC2, 
                           col = CDKN2AB_del,fill = CDKN2AB_del)) +
  ggbeeswarm::geom_quasirandom(dodge.width = .8, alpha = .3, aes(col = CDKN2AB_del)) +
  geom_boxplot(width = .4, alpha = .33, position = position_dodge(.8), outlier.shape = NA) +
  scale_color_manual(values = c("grey", "black")) +
  scale_fill_manual(values = c("grey", "black")) +
  labs(x = "", y = "Unsupervised PC2 (RNA)") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 6) +
  stat_n_text(size = 6) +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 15))

(boxplt1 + boxplt2 & theme(legend.position = "bottom")) + plot_layout(guides = "collect")
# ggsave("output/figures/paper/fig1/who_lgc_boxplot.pdf", width = 12, height = 7)



boxplt3 <- ggplot(pc2 %>% dplyr::filter(! is.na(who_2021)), aes(x = Study, y = PC2, 
                                                                col = who_2021,fill = who_2021)) +
  ggbeeswarm::geom_quasirandom(dodge.width = .8, alpha = .3, aes(col = who_2021)) +
  geom_boxplot(width = .4, alpha = .33, position = position_dodge(.8), outlier.shape = NA) +
  scale_color_manual(values = c("grey", "black")) +
  scale_fill_manual(values = c("grey", "black")) +
  labs(x = "", y = "Unsupervised PC2 (RNA)") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 8) +
  stat_n_text(size = 6) +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.position = "bottom")

boxplt4 <- ggplot(pc2 %>% dplyr::filter(! is.na(who_2021)), aes(x = Study, y = log_AIDH_cal_ratio, 
                                                                col = who_2021,fill = who_2021)) +
  ggbeeswarm::geom_quasirandom(dodge.width = .8, alpha = .3, aes(col = who_2021)) +
  geom_boxplot(width = .4, alpha = .33, position = position_dodge(.8), outlier.shape = NA) +
  scale_color_manual(values = c("grey", "black")) +
  scale_fill_manual(values = c("grey", "black")) +
  labs(x = "", y = "Linear Grading Coefficient (LGC, methylation)") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 8) +
  stat_n_text(size = 6) +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.position = "bottom")

(boxplt3 + boxplt4 & theme(legend.position = "bottom")) + plot_layout(guides = "collect")

# ggsave("output/figures/paper/sf1/clean/who_boxplot.pdf", width = 12, height = 7)

