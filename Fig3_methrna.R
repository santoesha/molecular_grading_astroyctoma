# load libs ----
library(tidyverse)

# load data ----
source("scripts/clean/dge/catnon/load_hclust.R")
source("scripts/clean/R/theme_cellpress.R")

hyp_probes <- readRDS("output/clean/hypermethylated_probes.overlap.Rds")
porath <- read.delim("data/pc2_targets_benporath.txt")

# CATNON ----
meta_cat <- readRDS("data/clean/catnon/rds/metadata.catnon.Rds") %>% 
            dplyr::filter(RNA_Seq == T) %>% 
            dplyr::filter(predictBrain_12.8_cal_class %in% c("A_IDH_LG", "A_IDH_HG")) %>% 
            dplyr::select(- mvalue_median_hyp)

mvalues_cat <- readRDS("data/clean/catnon/rds/mvalues.nomaskedprobes.catnon.Rds") %>% 
               dplyr::filter(probeID %in% hyp_probes) %>% 
               dplyr::select(meta_cat$BatchID) 

plt.cat <- data.frame(BatchID = colnames(mvalues_cat),
                      mvalue_median_hyp = matrixStats::colMedians(as.matrix(mvalues_cat))) %>% 
           dplyr::left_join(., meta_cat, by = c("BatchID" = "BatchID")) %>% 
           dplyr::mutate(Study = "CATNON") %>% 
           dplyr::select(mvalue_median_hyp, lts.emb_up, Study, CDKN2AB_del)

# TCGA ----

meta_tcga <- readRDS("data/clean/tcga/rds/metadata_tcga_idhmutnoncodel.Rds") %>% 
                     dplyr::filter(RNA_Seq == T) %>% 
                     dplyr::filter(predictBrain_12.8_cal_class %in% c("A_IDH_LG", "A_IDH_HG")) %>% 
             dplyr::select(- mvalue_median_hyp)

mvalues_tcga <- readRDS("data/clean/tcga/rds/mvalues.nomaskedprobes.tcga.Rds") %>% 
                       dplyr::filter(probeID %in% hyp_probes) %>% 
                       dplyr::select(meta_tcga$BatchID) 

plt.tcga <- data.frame(BatchID = colnames(mvalues_tcga),
                       mvalue_median_hyp = matrixStats::colMedians(as.matrix(mvalues_tcga))) %>% 
            dplyr::left_join(., meta_tcga) %>% 
            dplyr::mutate(Study = "TCGA") %>% 
            dplyr::select(mvalue_median_hyp, lts.emb_up, Study, CDKN2AB_del)

# GLASS ----
source("scripts/clean/dge/glass/transcriptional_signatures_glass.R")

meta_glass <- readRDS("data/clean/glass/rds/metadata.glass.Rds") %>% 
              dplyr::filter(RNA_Seq == T) %>% 
              dplyr::filter(predictBrain_12.8_cal_class %in% c("A_IDH_LG", "A_IDH_HG")) %>% 
              dplyr::select(BatchID, GLASS_ID, GS_ID, CDKN2AB_del, Resectie) %>% 
              dplyr::left_join(., transcriptional.signatures) %>% 
              dplyr::group_by(GLASS_ID) %>% 
              dplyr::mutate(n_samples = length(GLASS_ID)) 

mvalues_glass <- readRDS("data/clean/glass/rds/mvalues.nomaskedprobes.glass.Rds") %>% 
                 dplyr::filter(probeID %in% hyp_probes) %>% 
                 dplyr::select(meta_glass$BatchID) 

plt.glass <- data.frame(BatchID = colnames(mvalues_glass),
                        mvalue_median_hyp = matrixStats::colMedians(as.matrix(mvalues_glass))) %>% 
             dplyr::left_join(., meta_glass) %>% 
             dplyr::mutate(Study = ifelse(Resectie == 1, "GLASS-NL-P", "GLASS-NL-R")) %>% 
             dplyr::select(mvalue_median_hyp, lts.emb_up, Study, CDKN2AB_del)

# Plot ----

plt <- rbind(plt.cat, plt.tcga, plt.glass) %>% 
       dplyr::mutate(Study = factor(Study, levels = c("CATNON" , "TCGA", "GLASS-NL-P", "GLASS-NL-R")))

ggplot(plt, aes(x = mvalue_median_hyp, y = lts.emb_up)) +
      facet_wrap(~ Study, nrow = 2) +
      geom_point(aes(col = CDKN2AB_del), size = theme_cellpress_size * 2/3, alpha = 0.45) +
      ggpubr::stat_cor(method="spearman", aes(label = ..r.label..), cor.coef.name = "rho", size = theme_cellpress_size) +
      theme_classic() +
      labs(x = "Median M-value hypermethylated probes", y = "Embryonic development RNA signature value (C2)") +
      scale_color_manual(values = c("grey", "#FF0000")) +
      theme_cellpress

ggsave("output/figures/paper/rev1/f3/mvalue_rnasig_emb.pdf", width = (11.2 * 0.95) * 0.37, height = (8.5 * 0.95) * 0.5)

