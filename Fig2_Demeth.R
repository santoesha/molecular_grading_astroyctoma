pc2 <- readRDS("output/clean/combined/rds/unsupervised_pc2_rna.Rds")

meta.catnon <- readRDS("data/clean/catnon/rds/metadata.catnon.Rds") %>% 
               dplyr::filter(RNA_Seq == T) %>% 
               dplyr::rename(Sample_ID = GS_ID) %>% 
               dplyr::select(Sample_ID, demeth_score) 
meta.tcga <- readRDS("data/clean/tcga/rds/metadata_tcga_idhmutnoncodel.Rds") %>% 
             dplyr::filter(RNA_Seq == T) %>% 
             dplyr::rename(Sample_ID = Case) %>% 
             dplyr::select(Sample_ID, demeth_score)
meta.glass <- readRDS("data/clean/glass/rds/metadata.glass.Rds") %>% 
              dplyr::filter(RNA_Seq == T) %>% 
              dplyr::filter(Resectie != 1) %>% 
              dplyr::rename(Sample_ID = GS_ID) %>% 
              dplyr::select(Sample_ID, demeth_score)

meta_rna <- rbind(meta.catnon, meta.tcga, meta.glass) %>% 
            dplyr::left_join(., pc2, by = c("Sample_ID" = "Sample_ID"))

ggplot(data = meta_rna, aes(x = PC2, y = demeth_score)) +
       facet_wrap(~ Study) +
       geom_point(alpha = 0.8, size = 5, aes(fill =  CDKN2AB_del), pch = 21, col = "grey20") +
       scale_fill_manual(values = c("grey", "dodgerblue")) +
       geom_smooth(method = "lm", se=FALSE, col = "black") +
       ggpubr::stat_cor(label.y = 3.5, label.x = -23, col = "black", method="spearman", size = 7, cor.coef.name = "rho") +
       theme_bw() +
       theme(strip.text.x = element_text(size = 20),
             legend.position = "bottom") +
       labs(x = "Unsupervised PC2 (RNA)", 
            y = "Methylation score")
# 
# ggsave("output/figures/paper/fig2/clean/pc2_demeth.pdf", width = 12, height = 8)
