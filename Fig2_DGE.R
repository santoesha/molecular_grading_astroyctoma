# load libs ----
library(tidyverse)
library(msigdb)
library(msigdbr)
library(clusterProfiler)

source("scripts/clean/R/theme_cellpress.R")

# load data -----

source("scripts/clean/dge/catnon/load_hclust.R")

meta_cat <- readRDS("data/clean/catnon/rds/metadata.catnon.Rds") %>% 
            #dplyr::filter(RNA_Seq == T) %>% 
            dplyr::filter(predictBrain_12.8_cal_class %in% c("A_IDH_LG", "A_IDH_HG"))
meta_tcga <- readRDS("data/clean/tcga/rds/metadata_tcga_idhmutnoncodel.Rds") %>% 
             dplyr::filter(RNA_Seq == T) %>% 
            dplyr::filter(predictBrain_12.8_cal_class %in% c("A_IDH_LG", "A_IDH_HG"))

res_ord_cat <- readRDS("output/clean/catnon/rds/resord_catnon_heidi_num_lncrna.Rds")
dtable <- readRDS("output/clean/catnon/rds/dtable_catnon_heidi_num_lncrna.Rds")

e_vst_cat <- readRDS("data/clean/catnon/rds/expression.catnon.vst.Rds")
e_vst_tcga <- readRDS("data/clean/tcga/rds/expression.tcga.astr.vst.Rds") %>% 
              dplyr::select(meta_tcga$Case)

# Corplot ----

# ** CATNON -----
exp_de <- e_vst_cat %>%
          dplyr::select(meta_cat$GS_ID) %>% 
          dplyr::filter( gsub("^.+_","",rownames(.)) %in% dtable$gene_name) %>% 
          `rownames<-`(ifelse(rownames(.) == "ENSG00000285077_ARHGAP11B", "ENSG00000285077_ARHGAP11B.1", rownames(.) ))%>%
          `rownames<-`(., gsub("^.+_","",rownames(.))) 

labels <- clusters_catnon %>% 
           dplyr::mutate(value = T) %>% 
           tidyr::pivot_wider(names_from = cluster, values_from = value) %>% 
           tibble::column_to_rownames('gene_name') %>% 
           base::replace(is.na(.), FALSE)

h <- recursiveCorPlot::recursiveCorPlot(exp_de,labels, 
                                        5, 5,
                                        caption =  paste0("Genes: n=",       nrow(labels), ", ",
                                                          "Samples: n=",ncol(exp_de)))
h

# ggsave("output/figures/paper/fig2/clean/corplot_catnon.png", width=8,height=8)

# ** TCGA ----
sig_res_tcga <- readRDS("output/clean/tcga/rds/resord_tcga_heidi_num_lncrna.Rds") %>% 
                dplyr::filter(abs(log2FoldChange) > 0.5, padj < 0.01 )

exp_de_tcga <- e_vst_tcga %>%
               dplyr::filter( gsub("^.+_","",rownames(.)) %in% sig_res_tcga$gene_name) %>% 
              `rownames<-`(ifelse(rownames(.) == "ENSG00000285077_ARHGAP11B", "ENSG00000285077_ARHGAP11B.1", rownames(.) ))%>%
              `rownames<-`(., gsub("^.+_","",rownames(.))) 

labels_tcga <- data.frame(gene_name = rownames(exp_de_tcga)) %>%
           dplyr::mutate(de_catnon = ifelse(gene_name %in% clusters_catnon$gene_name, T, F)) %>% 
           column_to_rownames('gene_name')
 
l <- recursiveCorPlot::recursiveCorPlot(exp_de_tcga,labels_tcga, 
                                        5, 5,
                                        caption =  paste0("Genes: n=",       nrow(labels_tcga), ", ",
                                                          "Samples: n=",ncol(exp_de_tcga)))
l

#ggsave("output/figures/paper/sf3/corplot_tcga.png", width=8,height=8)



# Pathway -----

msigdbr_c5 <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
              dplyr::distinct(gs_name, gene_symbol) %>%
              as.data.frame()

msig_ecm <- enricher(gene = clusters_catnon %>% dplyr::filter(cluster == "col_up") %>% dplyr::pull(gene_name),
                         TERM2GENE = msigdbr_c5) %>%
            as_tibble() %>% 
            slice_min(qvalue, n=5) %>% 
            dplyr::mutate(qvalue = log2(qvalue))

msig_emb <- enricher(gene = clusters_catnon %>% dplyr::filter(cluster == "emb_up") %>% dplyr::pull(gene_name),
                     TERM2GENE = msigdbr_c5) %>%
            as_tibble() %>% 
            slice_min(qvalue, n=5) %>% 
            dplyr::mutate(qvalue = log2(qvalue))

msig_cc <- enricher(gene = clusters_catnon %>% dplyr::filter(cluster == "cell_cyc_up") %>% dplyr::pull(gene_name),
                     TERM2GENE = msigdbr_c5) %>%
           as_tibble() %>% 
           slice_min(qvalue, n=5) %>% 
           dplyr::mutate(qvalue = log2(qvalue))
        

p1 <- ggpubr::ggbarplot(msig_ecm, x = "ID", y = "qvalue",
                        col = "black",
                        fill = "red",
                        sort.val = "desc",          
                        rotate = TRUE,
                        sort.by.groups = FALSE,     
                        x.text.angle = 45,         
                        ylab = "log2(q-value)",
                        xlab = "",
                        lab.size =1) +
     coord_flip() +
     scale_y_continuous() +
     scale_x_discrete(position = "top") +
     theme(axis.text.y = element_text(size = 20, hjust = 0),  # Adjust hjust to align y-axis labels
          axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust =1 ),
          axis.title = element_text(size = 20)) 
  

p2 <- ggpubr::ggbarplot(msig_emb, x = "ID", y = "qvalue",
                        col = "black",
                        fill = "#ffd600",
                        sort.val = "desc",          
                        rotate = TRUE,
                        sort.by.groups = FALSE,     
                        x.text.angle = 45,         
                        ylab = "log2(q-value)",
                        xlab = "",
                        lab.size =1) +
     coord_flip() +
     scale_y_continuous() +
     scale_x_discrete(position = "top") +
     theme(axis.text.y = element_text(size = 20, hjust = 0),  # Adjust hjust to align y-axis labels
          axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust =1 ),
          axis.title = element_text(size = 20)) 

p3 <- ggpubr::ggbarplot(msig_cc, x = "ID", y = "qvalue",
                        col = "black",
                        fill = "#00ff80",
                        sort.val = "desc",          
                        rotate = TRUE,
                        sort.by.groups = FALSE,     
                        x.text.angle = 45,         
                        ylab = "log2(q-value)",
                        xlab = "",
                        lab.size =1) +
     coord_flip() +
     scale_y_continuous() +
     scale_x_discrete(position = "top") +
     theme(axis.text.y = element_text(size = 20, hjust = 0),  # Adjust hjust to align y-axis labels
          axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust =1 ),
          axis.title = element_text(size = 20)) 

p3/p2/p1
# 
# ggsave("output/figures/paper/fig2/clean/gsea_go_dge_rna.pdf", width=20,height=8)


# Volcano ----

# ** CATNON ----
res_ord_cat <- res_ord_cat %>% 
               dplyr::left_join(., clusters_catnon) %>% 
               dplyr::rename(stat_catnon = stat)

p4 <- ggplot(res_ord_cat, aes(x = log2FoldChange , y = -log10(padj))) +
      geom_point(data = subset(res_ord_cat, is.na(cluster)), size = 1.5, fill = "lightgrey", alpha = 0.8,pch=21, color = "NA") +
      geom_point(data = subset(res_ord_cat, cluster == "down.1"), size = 1.5, aes(fill = "C0"), alpha = 0.8,pch=21, color = "NA") +
      geom_point(data = subset(res_ord_cat, cluster == "cell_cyc_up"), size = 1.5, aes(fill = "C1"), alpha = 0.8,pch=21, color = "NA") +
      geom_point(data = subset(res_ord_cat, cluster == "emb_up"), size = 1.5, aes(fill = "C2"), alpha = 0.8,pch=21, color = "NA") +
      geom_point(data = subset(res_ord_cat, cluster == "col_up"), size = 1.5, aes(fill = "C3"), alpha = 0.8,pch=21, color = "NA") +
      geom_point(data = subset(res_ord_cat, cluster %in% c("up.4", "up.5")), size = 1.5, aes(fill = ""), alpha = 0.8,pch=21, color = "NA") +
      scale_fill_manual(values = c("grey30","#996632",  "#00ff80", "#ffd600", "red" )) +
      geom_hline(yintercept = -log10(0.01), linetype = "longdash") +
      geom_vline(xintercept = .5, linetype = "longdash") +
      geom_vline(xintercept = -.5, linetype = "longdash") +
      theme_classic() +
      xlim(-3.5, 3.5) +
      ylim(0, 31) +
      labs(x = "log2(Fold Change RNA) CGC", y = "-log10(FDR adjusted p-value)", title = "CATNON (n=132)") +
      theme_cellpress +
      theme(legend.position = "top", legend.title = element_blank())

# ** TCGA ----

res_ord_tcga <- readRDS("output/clean/tcga/rds/resord_tcga_heidi_num_lncrna.Rds") %>% 
                dplyr::left_join(., clusters_catnon) %>% 
                dplyr::rename(stat_tcga = stat)

p5 <- ggplot(res_ord_tcga, aes(x = log2FoldChange , y = -log10(padj))) +
      geom_point(data = subset(res_ord_tcga, is.na(cluster)), size = 1.5, fill = "lightgrey", alpha = 0.8,pch=21, color = "NA") +
      geom_point(data = subset(res_ord_tcga, cluster == "down.1"), size = 1.5, aes(fill = "C0"), alpha = 0.8,pch=21, color = "NA") +
      geom_point(data = subset(res_ord_tcga, cluster == "cell_cyc_up"), size = 1.5, aes(fill = "C1"), alpha = 0.8,pch=21, color = "NA") +
      geom_point(data = subset(res_ord_tcga, cluster == "emb_up"), size = 1.5, aes(fill = "C2"), alpha = 0.8,pch=21, color = "NA") +
      geom_point(data = subset(res_ord_tcga, cluster == "col_up"), size = 1.5, aes(fill = "C3"), alpha = 0.8,pch=21, color = "NA") +
      geom_point(data = subset(res_ord_tcga, cluster %in% c("up.4", "up.5")), size = 1.5, aes(fill = ""), alpha = 0.8,pch=21, color = "NA") +
      scale_fill_manual(values = c("grey30","#996632",  "#00ff80", "#ffd600", "red" )) +
      geom_hline(yintercept = -log10(0.01), linetype = "longdash") +
      geom_vline(xintercept = .5, linetype = "longdash") +
      geom_vline(xintercept = -.5, linetype = "longdash") +
      theme_classic() +
      xlim(-3.5, 3.5) +
      ylim(0, 50.5) +
      labs(x = "log2(Fold Change RNA) CGC", y = "-log10(FDR adjusted p-value)", title = "TCGA (n=235)") +
      theme_cellpress +
      theme(legend.position = "none", legend.title = element_blank())

p4+p5

# Correlation tests -----

res_ord_cat_cor <- res_ord_cat %>% 
                   dplyr::select(stat_catnon, gene_name, cluster)
res_ord_tcga_cor <- res_ord_tcga %>% 
                    dplyr::select(stat_tcga, gene_name)

plt <- res_ord_cat_cor %>% 
         dplyr::left_join(., res_ord_tcga_cor, by = c("gene_name" = "gene_name")) 

p6 <- ggplot(plt, aes(x = stat_catnon, y=stat_tcga)) +
      ggpubr::stat_cor(method="spearman", aes(label = ..r.label..), cor.coef.name = "rho", size = theme_cellpress_size) +
      geom_point(data = subset(plt, is.na(cluster)), size = 1.5, col = "grey70", alpha = 2,pch=21) +
      geom_point(data = subset(plt, cluster %in% c("emb_up", "down.1", "col_up", "cell_cyc_up")), aes(fill = cluster, col = cluster),
                 size = 1.5, alpha = 0.8,pch=21) +
      scale_fill_manual(values = c("#996632", "#00ff80","red",  "#ffd600", "grey70")) +
      scale_color_manual(values = c("#996632", "#00ff80","red",  "#ffd600", "grey70")) +
      labs(x = "Wald statistics CATNON RNA: CGC", y = "Wald statistics TCGA RNA: CGC") +
      #ggrepel::geom_text_repel(data = subset(plt, embryo_gene == T & cluster == "tf_up"), aes(label = gene_name), box.padding = 1, max.overlaps = 50, size = 2.5, col = "black") +
      theme_classic() +
      theme_cellpress +
      theme(legend.title = element_blank(),
            legend.position = "top") 

p4+p5+p6 + plot_layout(widths = c(2, 2, 1))

ggsave("output/figures/paper/rev1/f2/volcano_lgc.pdf", width = (11.2 * 0.95) * 0.7, height = (8.5 * 0.95) * 0.3) 


# ggsave("output/figures/paper/fig2/clean/correlation_de_tcgacatnon.pdf", width = 5, height = 7)

# Venn diagram -----

# ** CATNON ----

pie_catnon <- clusters_catnon %>% 
              dplyr::mutate(value = dplyr::case_when(is.na(cluster) ~ "other",
                                         cluster %in% c("up.4", "up.5") ~ "other",
                                         cluster == "down.1" ~ "C0",
                                         cluster == "cell_cyc_up" ~ "C1",
                                         cluster == "emb_up" ~ "C2",
                                         cluster == "col_up" ~ "C3"))
labels_catnon <- c("C0", "C1 (cell cycling)", "C2 (embryonic development)", "C3 (extracellular matrix)", "Other")
values_catnon <- c(length(which(pie_catnon$value == "C0")), 
                 length(which(pie_catnon$value == "C1")),
                 length(which(pie_catnon$value == "C2")),
                 length(which(pie_catnon$value == "C3")),
                 length(which(pie_catnon$value == "other")))

pdf(file = "output/figures/paper/fig2/clean/pie_catnon.pdf",
    width = 7,
    height = 5)
pie(values_catnon, labels = labels_catnon, col = c("#996632", "#00ff80","#ffd600" , "red" , "grey70"))
dev.off()

# ** TCGA ----

dtable_tcga <- readRDS("output/clean/tcga/rds/dtable_tcga_heidi_num_lncrna.Rds")

pie_tcga <- dtable_tcga %>% 
       dplyr::left_join(., clusters_catnon) %>% 
       dplyr::mutate(value = dplyr::case_when(is.na(cluster) ~ "other",
                                              cluster %in% c("up.4", "up.5") ~ "other",
                                              cluster == "down.1" ~ "C0",
                                              cluster == "cell_cyc_up" ~ "C1",
                                              cluster == "emb_up" ~ "C2",
                                              cluster == "col_up" ~ "C3"))

labels_tcga <- c("C0", "C1 (cell cycling)", "C2 (embryonic development)", "C3 (extracellular matrix)", "Other")
values_tcga <- c(length(which(pie_tcga$value == "C0")), 
            length(which(pie_tcga$value == "C1")),
            length(which(pie_tcga$value == "C2")),
            length(which(pie_tcga$value == "C3")),
            length(which(pie_tcga$value == "other")))

pdf(file = "output/figures/paper/fig2/clean/pie_tcga.pdf",
    width = 7,
    height = 5)
pie(values_tcga, labels = labels_tcga, col = c("#996632", "#00ff80","#ffd600" , "red" , "grey70"))
dev.off()
