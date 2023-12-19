
# ** CATNON ----
meta_cat <- readRDS("data/clean/catnon/rds/metadata.catnon.Rds") %>% 
            dplyr::filter(RNA_Seq == T) %>% 
            dplyr::select(CDKN2AB_del, lts.cell_cyc_up, lts.col_up, lts.down, lts.emb_up, demeth_score) %>% 
            pivot_longer(cols = c(lts.cell_cyc_up, lts.col_up, lts.down, lts.emb_up), names_to = "cluster", values_to = "rna_score") %>% 
            dplyr::mutate(Study = "CATNON")

meta_tcga <- readRDS("data/clean/tcga/rds/metadata_tcga_idhmutnoncodel.Rds") %>% 
             dplyr::filter(RNA_Seq == T)  %>% 
             dplyr::select(CDKN2AB_del, lts.cell_cyc_up, lts.col_up, lts.down, lts.emb_up, demeth_score) %>% 
             pivot_longer(cols = c(lts.cell_cyc_up, lts.col_up, lts.down, lts.emb_up), names_to = "cluster", values_to = "rna_score") %>% 
             dplyr::mutate(Study = "TCGA")

plt <- rbind(meta_cat, meta_tcga)

ggplot(plt, aes(x = demeth_score, y = rna_score)) +
      facet_wrap(~ cluster) +
      facet_grid(cluster~Study, scales = "free", space = "free") +
      geom_point() +
      ggpubr::stat_cor(method="spearman", aes(label = ..r.label..), cor.coef.name = "rho", size = 8) +
      theme_bw()
      

