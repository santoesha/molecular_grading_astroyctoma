# load libs ----
library(tidyverse)

# load data ----
source("scripts/clean/dge/catnon/load_hclust.R")
source("scripts/clean/R/theme_cellpress.R")

# Down cluster ----
down <- clusters_catnon %>% dplyr::filter(cluster == "down.1")

fetal_ast <- readxl::read_xlsx("data/clean/fpkm_human_celltypes.xlsx", sheet = 1) %>% 
             dplyr::filter(Gene %in% down$gene_name ) %>% 
             t() %>% 
             as.data.frame() %>% 
             janitor::row_to_names(row_number = 1)  %>% 
             dplyr::mutate_if(is.character, as.numeric) %>% 
             dplyr::rowwise() %>%
             dplyr::mutate(mean_fpkm = mean(c_across(where(is.numeric)), na.rm=TRUE)) %>% 
             dplyr::mutate(celltype = "fetal_astrocyte") %>% 
             dplyr::select(mean_fpkm, celltype)

mature_ast <- readxl::read_xlsx("data/clean/fpkm_human_celltypes.xlsx", sheet = 2) %>% 
              dplyr::filter(Gene %in% down$gene_name ) %>% 
              t() %>% 
              as.data.frame() %>% 
              janitor::row_to_names(row_number = 1)  %>% 
              dplyr::mutate_if(is.character, as.numeric) %>% 
              dplyr::rowwise() %>%
              dplyr::mutate(mean_fpkm = mean(c_across(where(is.numeric)), na.rm=TRUE)) %>% 
              dplyr::mutate(celltype = "mature_astrocyte") %>% 
              dplyr::select(mean_fpkm, celltype)

oligo <- readxl::read_xlsx("data/clean/fpkm_human_celltypes.xlsx", sheet = 3) %>% 
         dplyr::filter(Gene %in% down$gene_name ) %>% 
         t() %>% 
         as.data.frame() %>% 
         janitor::row_to_names(row_number = 1)  %>% 
         dplyr::mutate_if(is.character, as.numeric) %>% 
         dplyr::rowwise() %>%
         dplyr::mutate(mean_fpkm = mean(c_across(where(is.numeric)), na.rm=TRUE)) %>% 
         dplyr::mutate(celltype = "oligodendrocyte") %>% 
         dplyr::select(mean_fpkm, celltype)

tam <- readxl::read_xlsx("data/clean/fpkm_human_celltypes.xlsx", sheet = 4) %>% 
       dplyr::filter(Gene %in% down$gene_name ) %>% 
       t() %>% 
       as.data.frame() %>% 
       janitor::row_to_names(row_number = 1)  %>%
       dplyr::mutate_if(is.character, as.numeric) %>% 
       dplyr::rowwise() %>%
       dplyr::mutate(mean_fpkm = mean(c_across(where(is.numeric)), na.rm=TRUE)) %>% 
       dplyr::mutate(celltype = "tam") %>% 
       dplyr::select(mean_fpkm, celltype)

df <- rbind(fetal_ast, mature_ast, oligo, tam)

comp <- list( c("mature_astrocyte", "fetal_astrocyte"), 
              c("mature_astrocyte", "oligodendrocyte"), 
              c("mature_astrocyte", "tam") )


p1 <- ggplot(df, aes(x = celltype, y = mean_fpkm)) +
      ggbeeswarm::geom_quasirandom(dodge.width = .8, alpha = .3) +
      geom_boxplot(width = .4, alpha = .33, position = position_dodge(.8), outlier.shape = NA) +
      labs(x = "Cell type (Zhang et al., 2016)", y = "Mean FPKM (downregulated genes)") +
      ggpubr::stat_compare_means(comparisons = comp, method = "wilcox.test", label = "p.format") +
      EnvStats::stat_n_text() +
      theme_cellpress

rm(df, comp, tam, oligo, fetal_ast, mature_ast, down)

# Up cluster
up <- clusters_catnon %>% dplyr::filter(cluster != "down.1")

fetal_ast <- readxl::read_xlsx("data/clean/fpkm_human_celltypes.xlsx", sheet = 1) %>% 
             dplyr::filter(Gene %in% up$gene_name ) %>%
             dplyr::slice(-c(211,212)) %>% 
             t() %>% 
             as.data.frame() %>% 
             janitor::row_to_names(row_number = 1)  %>% 
             dplyr::mutate_if(is.character, as.numeric) %>% 
             dplyr::rowwise() %>%
             dplyr::mutate(mean_fpkm = mean(c_across(where(is.numeric)), na.rm=TRUE)) %>% 
             dplyr::mutate(celltype = "fetal_astrocyte") %>% 
             dplyr::select(mean_fpkm, celltype)

mature_ast <- readxl::read_xlsx("data/clean/fpkm_human_celltypes.xlsx", sheet = 2) %>% 
              dplyr::filter(Gene %in% up$gene_name ) %>% 
              dplyr::slice(-c(211,212)) %>% 
              t() %>% 
              as.data.frame() %>% 
              janitor::row_to_names(row_number = 1)  %>% 
              dplyr::mutate_if(is.character, as.numeric) %>% 
              dplyr::rowwise() %>%
              dplyr::mutate(mean_fpkm = mean(c_across(where(is.numeric)), na.rm=TRUE)) %>% 
              dplyr::mutate(celltype = "mature_astrocyte") %>% 
              dplyr::select(mean_fpkm, celltype)

oligo <- readxl::read_xlsx("data/clean/fpkm_human_celltypes.xlsx", sheet = 3) %>% 
         dplyr::filter(Gene %in% up$gene_name ) %>% 
         dplyr::slice(-c(211,212)) %>% 
         t() %>% 
         as.data.frame() %>% 
         janitor::row_to_names(row_number = 1)  %>% 
         dplyr::mutate_if(is.character, as.numeric) %>% 
         dplyr::rowwise() %>%
         dplyr::mutate(mean_fpkm = mean(c_across(where(is.numeric)), na.rm=TRUE)) %>% 
         dplyr::mutate(celltype = "oligodendrocyte") %>% 
         dplyr::select(mean_fpkm, celltype)

tam <- readxl::read_xlsx("data/clean/fpkm_human_celltypes.xlsx", sheet = 4) %>% 
       dplyr::filter(Gene %in% up$gene_name ) %>% 
       dplyr::slice(-c(211,212)) %>% 
       t() %>% 
       as.data.frame() %>% 
       janitor::row_to_names(row_number = 1)  %>%
       dplyr::mutate_if(is.character, as.numeric) %>% 
       dplyr::rowwise() %>%
       dplyr::mutate(mean_fpkm = mean(c_across(where(is.numeric)), na.rm=TRUE)) %>% 
       dplyr::mutate(celltype = "tam") %>% 
       dplyr::select(mean_fpkm, celltype)

df <- rbind(fetal_ast, mature_ast, oligo, tam)

comp <- list( c("fetal_astrocyte", "mature_astrocyte"), 
              c("fetal_astrocyte", "oligodendrocyte"), 
              c("fetal_astrocyte", "tam") )


p2 <- ggplot(df, aes(x = celltype, y = mean_fpkm)) +
      ggbeeswarm::geom_quasirandom(dodge.width = .8, alpha = .3) +
      geom_boxplot(width = .4, alpha = .33, position = position_dodge(.8), outlier.shape = NA) +
      labs(x = "Cell type (Zhang et al., 2016)", y = "Mean FPKM (upregulated genes)") +
      ggpubr::stat_compare_means(comparisons = comp, method = "wilcox.test", label = "p.format") +
      EnvStats::stat_n_text() +
      theme_cellpress

p1+p2

ggsave("output/figures/paper/sf3/celltype_healthy_brain_clusters.pdf", width = (8.5 * 0.95), height = 4)

