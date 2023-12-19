# load libs ----
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
source("scripts/clean/R/theme_cellpress.R")

# load data -----

ecm_markers <- read.csv("data/GO:0005201_ECM.csv") %>% 
               dplyr::rename(gene_name = name ) %>% 
               dplyr::pull(gene_name)

his_markers <- read.csv("data/GO:0030527_HIS.csv") %>% 
               dplyr::rename(gene_name = name ) %>% 
               dplyr::pull(gene_name)

emb_markers <- read.csv("data/GO:0048598_EMB.csv") %>% 
               dplyr::rename(gene_name = name ) %>% 
               dplyr::pull(gene_name)

source("scripts/clean/dge/catnon/load_hclust.R")

source("scripts/clean/dge/signature_genes.R")

# ** CATNON -----
meta.catnon <- readRDS("data/clean/catnon/rds/metadata.catnon.Rds") %>% 
               dplyr::filter(RNA_Seq == T)
counts.catnon <- readRDS("data/clean/catnon/rds/expression.catnon.vst.Rds") %>% 
                 dplyr::select(meta.catnon$GS_ID)

dtable_cat <- readRDS("")
  
# ** TCGA -----
meta.tcga <- readRDS("data/clean/tcga/rds/metadata_tcga_idhmutnoncodel.Rds") %>% 
              dplyr::filter(RNA_Seq == T)
counts.tcga <- readRDS("data/clean/tcga/rds/expression.tcga.astr.vst.Rds") %>% 
               dplyr::select(meta.tcga$Case)
# Corplot -----

# ** CATNON -----

exp_de_catnon <- counts.catnon %>%
                 dplyr::filter( rownames(.) %in% dtable_cat$gene_uid) %>% 
                 `rownames<-`(., gsub("^.+_","",rownames(.)))

labels_cat <- dtable_cat %>% 
              remove_rownames(.) %>% 
              dplyr::select(gene_name,up_gene)  %>% 
              dplyr::mutate(emb_gene = ifelse(gene_name %in% emb_markers, T, F))  %>% 
              dplyr::mutate(ecm_markers = ifelse(gene_name %in% ecm_markers, T, F)) %>% 
              tibble::column_to_rownames('gene_name') %>% 
              dplyr::mutate(gene = NULL) 

h <- recursiveCorPlot::recursiveCorPlot(exp_de_catnon,labels_cat, 
                                        5, 5,
                                        caption =  paste0("Genes: n=",       nrow(dtable_cat), ", ",
                                                          "CATNON samples: n=",ncol(exp_de_catnon)))
h

ggsave("output/figures/paper/fig2/clean/recursiveCorPlot_catnon.png", height=7, width=7, device = "png", dpi = 1000)

# i <- recursiveCorPlot::recursiveCorPlot(exp_de_catnon,labels_cat,3.5,2,return_h_object = T)
# write_rds(i, "output/clean/catnon/rds/corplot_clusters_numeric_cat.Rds")


# ** TCGA -----

exp_de_tcga <- counts.tcga %>%
  dplyr::filter( rownames(.) %in% dtable_tcga$gene_uid) %>% 
  `rownames<-`(., gsub("^.+_","",rownames(.)))

labels_tcga <- dtable_tcga %>% 
  remove_rownames(.) %>% 
  dplyr::select(gene_name,up_gene)  %>% 
  dplyr::mutate(emb_gene = ifelse(gene_name %in% emb_markers, T, F))  %>% 
  dplyr::mutate(ecm_markers = ifelse(gene_name %in% ecm_markers, T, F)) %>% 
  tibble::column_to_rownames('gene_name') %>% 
  dplyr::mutate(gene = NULL) 

j <- recursiveCorPlot::recursiveCorPlot(exp_de_tcga,labels_tcga, 
                                        5, 5,
                                        caption =  paste0("Genes: n=",       nrow(dtable_tcga), ", ",
                                                          "TCGA samples: n=",ncol(exp_de_tcga)))
j

c
# 
# k <- recursiveCorPlot::recursiveCorPlot(exp_de_tcga,labels_tcga,3.5,2,return_h_object = T)
# write_rds(k, "output/clean/tcga/rds/corplot_clusters_numeric_tcga.Rds")

# Pathway -----

source("scripts/clean/dge/catnon/load_hclust.R")
source("scripts/clean/dge/tcga/load_hclust.R")

msigdbr_c5 <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
              dplyr::distinct(gs_name, gene_symbol) %>%
              as.data.frame()

# ** ECM -----
msig_ecm_cat <- enricher(gene = clusters_catnon %>% dplyr::filter(cluster == "col_up") %>% dplyr::pull(gene_name),
                     TERM2GENE = msigdbr_c5) %>%
                    as_tibble() %>% 
                    slice_min(qvalue, n=5) %>% 
                    dplyr::mutate(qvalue = log2(qvalue)) %>% 
                    dplyr::mutate(cluster = "C3")
msig_ecm_tcga <- enricher(gene = clusters_tcga %>% dplyr::filter(cluster == "col_up") %>% dplyr::pull(gene_name),
                          TERM2GENE = msigdbr_c5) %>%
                as_tibble() %>% 
                slice_min(qvalue, n=5) %>% 
                dplyr::mutate(qvalue = log2(qvalue))

p1 <- ggpubr::ggbarplot(msig_ecm_cat, x = "ID", y = "qvalue",
                        col = "black",
                        fill = "red",
                        sort.val = "desc",          
                        #rotate = TRUE,
                        sort.by.groups = FALSE,     
                        x.text.angle = 45,         
                        ylab = "log2(q-value)",
                        xlab = "",
                        lab.size = theme_cellpress_size) +
       theme_cellpress +
       coord_flip()


p2 <- ggpubr::ggbarplot(msig_ecm_tcga, x = "ID", y = "qvalue",
                        col = "black",
                        fill = "red",
                        sort.val = "desc",          
                        rotate = TRUE,
                        sort.by.groups = FALSE,     
                        x.text.angle = 45,         
                        ylab = "log2(q-value)",
                        xlab = "",
                        lab.size =theme_cellpress_size) +
  theme_cellpress

# ** EMB -----

msig_emb_cat <- enricher(gene = clusters_catnon %>% dplyr::filter(cluster == "emb_up") %>% dplyr::pull(gene_name),
                         TERM2GENE = msigdbr_c5) %>%
                as_tibble() %>% 
                slice_min(qvalue, n=5) %>% 
                dplyr::mutate(qvalue = log2(qvalue)) %>% 
                dplyr::mutate(cluster = "C2")
msig_emb_tcga <- enricher(gene = clusters_tcga %>% dplyr::filter(cluster == "emb_up") %>% dplyr::pull(gene_name),
                          TERM2GENE = msigdbr_c5) %>%
                  as_tibble() %>% 
                  slice_min(qvalue, n=5) %>% 
                  dplyr::mutate(qvalue = log2(qvalue))

p3 <- ggpubr::ggbarplot(msig_emb_cat, x = "ID", y = "qvalue",
                        col = "black",
                        fill = "#ffd600",
                        sort.val = "desc",          
                        rotate = TRUE,
                        sort.by.groups = FALSE,     
                        x.text.angle = 45,         
                        ylab = "log2(q-value)",
                        xlab = "",
                        lab.size =theme_cellpress_size) +
  theme_cellpress


p4 <- ggpubr::ggbarplot(msig_emb_tcga, x = "ID", y = "qvalue",
                        col = "black",
                        fill = "#ffd600",
                        sort.val = "desc",          
                        rotate = TRUE,
                        sort.by.groups = FALSE,     
                        x.text.angle = 45,         
                        ylab = "log2(q-value)",
                        xlab = "",
                        lab.size =theme_cellpress_size) +
  theme_cellpress


# ** EMB -----

msig_cyc_cat <- enricher(gene = clusters_catnon %>% dplyr::filter(cluster == "cell_cyc_up") %>% dplyr::pull(gene_name),
                         TERM2GENE = msigdbr_c5) %>%
  as_tibble() %>% 
  slice_min(qvalue, n=5) %>% 
  dplyr::mutate(qvalue = log2(qvalue)) %>% 
  dplyr::mutate(cluster = "C1")
msig_cyc_tcga <- enricher(gene = clusters_tcga %>% dplyr::filter(cluster == "cell_cyc_up") %>% dplyr::pull(gene_name),
                          TERM2GENE = msigdbr_c5) %>%
  as_tibble() %>% 
  slice_min(qvalue, n=5) %>% 
  dplyr::mutate(qvalue = log2(qvalue))

p5 <- ggpubr::ggbarplot(msig_cyc_cat, x = "ID", y = "qvalue",
                        col = "black",
                        fill = "#00ff80",
                        sort.val = "desc",          
                        rotate = TRUE,
                        sort.by.groups = FALSE,     
                        x.text.angle = 45,         
                        ylab = "log2(q-value)",
                        xlab = "",
                        lab.size =theme_cellpress_size) +
  theme_cellpress


p6 <- ggpubr::ggbarplot(msig_cyc_tcga, x = "ID", y = "qvalue",
                        col = "black",
                        fill = "#00ff80",
                        sort.val = "desc",          
                        rotate = TRUE,
                        sort.by.groups = FALSE,     
                        x.text.angle = 45,         
                        ylab = "log2(q-value)",
                        xlab = "",
                        lab.size =theme_cellpress_size) +
  theme_cellpress

p5+p6

# plt ----

plt <- rbind(msig_cyc_cat, msig_ecm_cat, msig_emb_cat) 

ggplot(data = plt) +
      facet_wrap(~ cluster, scale="free") +
      geom_bar(aes(x = reorder(ID,qvalue), y = qvalue, color = cluster, fill = cluster), stat = "identity") +
      scale_color_manual(values = c("#00ff80",  "#ffd600", "red")) +
      scale_fill_manual(values = c("#00ff80", "#ffd600", "red")) +
      theme_cellpress +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.position = "none") +
      labs(x = "", y = "q-value")
  
ggsave("output/figures/paper/rev1/f2/barplot_methylation_dmp.pdf", width = 30, height = 11)


ggsave("output/figures/paper/rev1/f2/barplot_pathway.pdf", width = 20, height = 13)
