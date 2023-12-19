# load libs ----
library(tidyverse)

# load data -----
source("scripts/R/probe_to_gene.R")
source("scripts/clean/R/sel_genes_dge.R")
source("scripts/clean/dge/catnon/load_hclust.R")
source("scripts/clean/R/theme_cellpress.R")

porath <- read.delim("data/pc2_targets_benporath.txt")
g_cimp_low_probes <- readRDS("data/rds/gcimp_low_probes.Rds")

res_ord_met_cat <- readRDS("output/clean/catnon/rds/methylation_res_cat.Rds") %>% 
                   dplyr::mutate(significant = ifelse(abs(logFC) > 0.5 & adj.P.Val < 0.01, T, F)) %>% 
                   dplyr::mutate(sign = dplyr::case_when(significant == T & logFC < -0.5 ~ "hypomethylated",
                                                         significant == T & logFC > 0.5 ~ "hypermethylated",
                                                         T ~ NA)) 

sig_met_cat <- res_ord_met_cat %>% 
               dplyr::filter(abs(logFC) > 0.5, adj.P.Val < 0.01) %>% 
               dplyr::mutate(sign = ifelse(logFC < 0, "hypomethylated", "hypermethylated"))
res_ord_cat <- readRDS("output/clean/catnon/rds/resord_catnon_heidi_num_lncrna.Rds") 

res_ord_met_tcga <- readRDS("output/clean/tcga/rds/methylation_res_tcga.Rds") %>% 
                    dplyr::mutate(significant = ifelse(abs(logFC) > 0.5 & adj.P.Val < 0.01, T, F)) %>% 
                    dplyr::mutate(sign = ifelse(significant == T & logFC < 0, "hypomethylated", "hypermethylated"))
sig_met_tcga <- res_ord_met_tcga %>% 
                dplyr::filter(abs(logFC) > 0.5, adj.P.Val < 0.01) %>% 
                dplyr::mutate(sign = ifelse(logFC < 0, "hypomethylated", "hypermethylated"))

res_ord_tcga <- readRDS("output/clean/tcga/rds/resord_tcga_heidi_num_lncrna.Rds")

annoEPIC <- read.csv("data/methylation/annoEPIC.csv", row.names = c(1))  %>% 
            dplyr::filter(! chr %in% c("chrX", "chrY")) %>% 
            dplyr::filter(UCSC_RefGene_Name != "")

sig_rna_to_met <- gene_to_probes %>% 
                  dplyr::filter(value %in% clusters_catnon$gene_name) %>% 
                  dplyr::left_join(., annoEPIC %>% dplyr::select(Name, Relation_to_Island))

input_probe_overlap <- intersect(res_ord_met_cat$probeID, res_ord_met_tcga$probeID)


ggplot(res_ord_met_cat, aes(x = logFC , y = -log10(adj.P.Val))) +
      geom_point(data = subset(res_ord_met_cat, significant == F), size = 2, fill = "lightgrey", alpha = 0.8,pch=21, color = "NA") +
      geom_point(data = subset(res_ord_met_cat, sign == "hypermethylated"), size = 2, aes(fill = "hypermethylated"), alpha = 0.8,pch=21, color = "NA") +
      geom_point(data = subset(res_ord_met_cat, sign == "hypomethylated"), size = 2, aes(fill = "hypomethylated"), alpha = 0.8,pch=21, color = "NA") +
      geom_hline(yintercept = -log10(0.01), linetype = "longdash") +
      geom_vline(xintercept = .5, linetype = "longdash") +
      geom_vline(xintercept = -.5, linetype = "longdash") +
      theme_cellpress +
      scale_fill_manual(values = c("orange","powderblue"))  +
      labs(x = "log2(Fold Change CGC)", y = "-log10(FDR adjusted p-value") 
# 
# ggsave("output/figures/volcano_dmp.pdf", width = (8.5 * 0.95), height = 5)

# Barplot ----


bar_cat_hypo <- sig_met_cat %>% 
                  dplyr::select(probeID, sign) %>% 
                  dplyr::group_by(sign) %>% 
                  dplyr::summarise(count = n()) %>% 
                  dplyr::ungroup() %>% 
                  dplyr::mutate(fraction = count/sum(count)) %>% 
                  dplyr::mutate(Study = "CATNON") 

bar_tcga_hypo <- sig_met_tcga %>% 
                 dplyr::select(probeID, sign) %>% 
                 dplyr::group_by(sign) %>% 
                 dplyr::summarise(count = n()) %>% 
                 dplyr::ungroup() %>% 
                 dplyr::mutate(fraction = count/sum(count)) %>% 
                 dplyr::mutate(Study = "TCGA") 

plt_hypo <- rbind(bar_cat_hypo, bar_tcga_hypo)                

p1 <- ggplot(plt_hypo, aes(fill=sign, x=Study, y = fraction)) + 
      geom_bar(stat="identity") +
      theme_classic() +
      labs(x = "", y = "Fraction of significant probes") +
      geom_text(aes(label=round(fraction,2), group=sign),size=theme_cellpress_size, position = position_fill(vjust=0.5), col = "black") +
      scale_fill_manual(values = c("#0074D9", "#FFD700")) +
      theme_cellpress +
      theme(legend.position = "bottom",
            legend.title = element_blank()) 


bar_expected <-   annoEPIC %>% 
                  dplyr::select(Relation_to_Island) %>% 
                  dplyr::group_by(Relation_to_Island) %>% 
                  dplyr::summarise(count = n()) %>% 
                  dplyr::ungroup() %>% 
                  dplyr::mutate(fraction = count/nrow(annoEPIC)) %>% 
                  dplyr::mutate(sign = "expected") %>% 
                  dplyr::mutate(label = "Distribution across genome")
                
bar_cat_island <- sig_met_cat %>% 
                  dplyr::select(probeID, sign) %>% 
                  dplyr::left_join(., annoEPIC, by = c("probeID" = "Name")) %>% 
                  dplyr::group_by(Relation_to_Island, sign) %>% 
                  dplyr::summarise(count = n()) %>% 
                  dplyr::ungroup() %>% 
                  dplyr::group_by(sign) %>% 
                  dplyr::mutate(fraction = count/sum(count)) %>% 
                  dplyr::mutate(label = "Distribution across DMP")
bar_cat_island <- rbind(bar_cat_island, bar_expected) 

bar_tcga_island <- sig_met_tcga %>% 
                   dplyr::select(probeID, sign) %>% 
                   dplyr::left_join(., annoEPIC, by = c("probeID" = "Name")) %>% 
                   dplyr::group_by(Relation_to_Island, sign) %>% 
                   dplyr::summarise(count = n()) %>% 
                   dplyr::ungroup() %>% 
                   dplyr::group_by(sign) %>% 
                   dplyr::mutate(fraction = count/sum(count)) %>% 
                   dplyr::mutate(label = "Distribution across DMP")
bar_tcga_island <- rbind(bar_tcga_island, bar_expected) 

plt_island <- rbind(bar_cat_island, bar_tcga_island)

p2 <- ggplot(bar_expected , aes(fill=Relation_to_Island, x=sign, y = fraction)) + 
      geom_bar(stat="identity") +
      theme_classic() +
      labs(x = "", y = "") +
      geom_text(aes(label=round(fraction,3), group=Relation_to_Island),size=theme_cellpress_size, position = position_fill(vjust=0.5)) +
      scale_fill_brewer(palette = "Dark2") +
      theme_cellpress +
      theme_void() +
      theme(legend.position = "none",
            legend.title = element_blank()) 

p3 <- ggplot(bar_cat_island %>% dplyr::filter(sign != "expected"), aes(fill=Relation_to_Island, x=sign, y = fraction)) + 
      geom_bar(stat="identity") +
      theme_classic() +
      labs(x = "", y = "") +
      geom_text(aes(label=round(fraction,3), group=Relation_to_Island),size=theme_cellpress_size, position = position_fill(vjust=0.5)) +
      scale_fill_brewer(palette = "Dark2") +
      theme_cellpress +
      theme(legend.position = "none",
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 45, vjust =1, hjust=1)) 

p4 <- ggplot(bar_tcga_island %>% dplyr::filter(sign != "expected"), aes(fill=Relation_to_Island, x=sign, y = fraction)) + 
      geom_bar(stat="identity") +
      theme_classic() +
      geom_text(aes(label=round(fraction,3), group=Relation_to_Island),size=theme_cellpress_size, position = position_fill(vjust=0.5)) +
      scale_fill_brewer(palette = "Dark2") +
      theme_cellpress +
      theme(legend.position = "right",
            legend.title = element_blank(),
            axis.text.x = element_text(angle = 45, vjust =1, hjust=1)) +
      labs(x = "", y = "") 

plt2 <- (p2)+(p3+p4) + plot_layout(widths = c(1, 3))

plt2
# 
# ggsave("output/figures/paper/rev1/f2/barplot_methylation_dmp.pdf", width = (11.2 * 0.95) * 1/3, height = (8.5 * 0.95) * (1/3))

# Statistical test ------

sig_cat_island <- sig_met_cat %>% 
                  dplyr::select(probeID, sign) %>% 
                  dplyr::left_join(., annoEPIC, by = c("probeID" = "Name")) 
                  dplyr::mutate(cond = ifelse(Relation_to_Island == "Island", "Island", "No Island"))


fisher.test(table(sig_cat_island$sign, sig_cat_island$cond), alternative = "greater")

obs_cat <- as.data.frame(table(sig_cat_island$Relation_to_Island)) %>% dplyr::rename(n_obs = Freq)
exp_cat <- as.data.frame(table(annoEPIC$Relation_to_Island)) %>% dplyr::rename(n_exp = Freq)
test_cat <- obs_cat %>% 
            dplyr::left_join(., exp_cat) %>% 
            dplyr::mutate(frac_obs = n_obs / nrow(sig_met_cat)) %>% 
            dplyr::mutate(frac_exp = (n_exp / nrow(annoEPIC))) 

chisq.test(test_cat$n_obs, p=test_cat$frac_exp)

sig_tcga_island <- sig_met_tcga %>% 
                   dplyr::select(probeID, sign) %>% 
                   dplyr::left_join(., annoEPIC, by = c("probeID" = "Name")) %>% 
                   dplyr::mutate(cond = ifelse(Relation_to_Island == "Island", "Island", "No Island"))

fisher.test(table(sig_tcga_island$sign, sig_tcga_island$cond), alternative = "greater")


# Diff probes to genes -----

dmp_to_gene_cat <- sig_met_cat %>% 
                   dplyr::left_join(., gene_to_probes, by = c("probeID" = "Name"))
length(unique(dmp_to_gene_cat$value))

dmp_to_gene_tcga <- sig_met_tcga %>% 
                    dplyr::left_join(., gene_to_probes, by = c("probeID" = "Name"))
length(unique(dmp_to_gene_tcga$value))

length(intersect(unique(dmp_to_gene_cat$value), unique(dmp_to_gene_tcga$value)))

a <- as.data.frame(unique(dmp_to_gene_cat$value))

# Hypermethylated -----
hyper_cat  <- sig_met_cat %>% 
              dplyr::filter(sign == "hypermethylated")

hyper_tcga  <- sig_met_tcga %>% 
               dplyr::filter(sign == "hypermethylated") 

hypo_cat  <- sig_met_cat %>% 
             dplyr::filter(sign == "hypomethylated")

hypo_tcga  <- sig_met_tcga %>% 
              dplyr::filter(sign == "hypomethylated") 

# saveRDS(intersect(hyper_cat$probeID, hyper_tcga$probeID), "output/clean/hypermethylated_probes.overlap.Rds")

# Hypermethylated + upregulated -----

hyper_cat  <- sig_met_cat %>% 
              dplyr::filter(sign == "hypermethylated") %>% 
              dplyr::filter(probeID %in% sig_rna_to_met$Name) %>% 
              dplyr::left_join(sig_rna_to_met, by = c("probeID" = "Name")) %>% 
              dplyr::filter(value %in% (clusters_catnon %>% 
                                         dplyr::filter(! cluster == "down.1") %>% 
                                         dplyr::pull(gene_name))) %>% 
              dplyr::mutate(Study = "CATNON") %>% 
              dplyr::left_join(., res_ord_cat %>% dplyr::select(gene_name, log2FoldChange), by = c("value" = "gene_name"))

length(unique(hyper_cat$value))

hyper_tcga  <- sig_met_tcga %>% 
               dplyr::filter(sign == "hypermethylated") %>% 
               dplyr::filter(probeID %in% sig_rna_to_met$Name) %>% 
               dplyr::left_join(sig_rna_to_met, by = c("probeID" = "Name")) %>% 
               dplyr::filter(value %in% (clusters_catnon %>% 
                                        dplyr::filter(! cluster == "down.1") %>% 
                                        dplyr::pull(gene_name))) %>% 
               dplyr::mutate(Study = "TCGA") %>% 
               dplyr::left_join(., res_ord_tcga %>% dplyr::select(gene_name, log2FoldChange), by = c("value" = "gene_name"))

length(unique(hyper_tcga$value))

length(intersect(hyper_cat$value, hyper_tcga$value))

tab_cat <- rbind(hyper_cat, hyper_tcga) %>% 
           dplyr::filter(value %in% intersect(hyper_cat$value, hyper_tcga$value)) %>%
           dplyr::mutate(label = paste0(probeID, " (", Relation_to_Island, ")")) %>% 
           dplyr::select(Study, logFC, log2FoldChange, value, label) %>% 
           dplyr::group_by(Study, value) %>% 
           dplyr::mutate(median_lfc_met = median(logFC)) 

# saveRDS(tab_cat, "output/clean/hypermeth_upreg_values.Rds")

label_tab <- tab_cat %>% 
             dplyr::group_by(value, Study) %>%
             dplyr::summarise(label = paste(label, collapse = ",")) %>% 
             dplyr::left_join(., tab_cat %>% dplyr::select(-label, -logFC), by = c("value", "Study")) %>% 
             distinct() %>% 
             dplyr::mutate(log2FoldChange = round(log2FoldChange, 2)) %>% 
             dplyr::mutate(median_lfc_met = round(median_lfc_met, 2)) %>% 
             dplyr::filter(log2FoldChange > 0.25) %>% 
             dplyr::group_by(value) %>% 
             dplyr::filter(n() > 1)

# write.csv(label_tab, "output/clean/hypermet_upreg.table.csv", row.names = F)

msigdbr_c2 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP") %>% 
              dplyr::distinct(gs_name, gene_symbol) %>%
              as.data.frame()

msig_hyp <- enricher(gene = unique(label_tab$value),
                     TERM2GENE = msigdbr_c2, 
                     universe = res_ord_cat$gene_name) %>%
            as_tibble() %>% 
            slice_min(qvalue, n=25) %>% 
            dplyr::mutate(qvalue = log2(qvalue))

a <- enricher(gene = clusters_catnon %>% dplyr::filter(cluster == "emb_up") %>% dplyr::pull(gene_name),
         TERM2GENE = msigdbr_c2, 
         universe = res_ord_cat$gene_name) %>%
  as_tibble() %>% 
  slice_min(qvalue, n=50) %>% 
  dplyr::mutate(qvalue = log2(qvalue))

# ** Porath -----

test_porath <- data.frame(gene = unique(c(porath$prc2_target, label_tab$value))) %>% 
               dplyr::mutate(porath = ifelse(gene %in% porath$prc2_target, "prc2_target", "no_prc2_target")) %>% 
               dplyr::mutate(hm_up = ifelse(gene %in% label_tab$value, "hm_up", "no_hm_up"))

fisher.test(table(test_porath$porath, test_porath$hm_up))


# ** DMP plot -----

dmp_cat_plt <- res_ord_met_cat %>% 
               dplyr::left_join(., annoEPIC, by = c("probeID" = "Name")) %>% 
               dplyr::left_join(., gene_to_probes, by = c("probeID" = "Name")) %>% 
               dplyr::filter(value %in% res_ord_cat$gene_name) %>% 
               dplyr::group_by(value) %>% 
               dplyr::summarise(median_stat_met = median(stat_meth)) %>% 
               dplyr::left_join(., res_ord_cat, by = c("value" = "gene_name")) %>% 
               dplyr::left_join(., clusters_catnon, by = c("value" = "gene_name"))  %>% 
               dplyr::mutate(cluster_name = dplyr::case_when(is.na(cluster) ~ NA,
                                                     cluster %in% c("up.4", "up.5") ~ "other",
                                                     cluster == "down.1" ~ "C0",
                                                     cluster == "cell_cyc_up" ~ "C1",
                                                     cluster == "emb_up" ~ "C2",
                                                     cluster == "col_up" ~ "C3"))

p3 <- ggplot(dmp_cat_plt, aes(x=median_stat_met, y=stat, col = cluster_name)) +
      geom_point(data = subset(dmp_cat_plt, is.na(cluster_name)), alpha = 0.5, col = "grey90") +
      geom_point(data = subset(dmp_cat_plt, ! is.na(cluster_name)), alpha = 0.5) +
      scale_color_manual(values = c("#996632", "#00ff80",  "#ffd600", "red","grey90", "grey90")) +
      geom_vline(xintercept = 0, col = "black", linewidth = theme_cellpress_lwd) +
      geom_hline(yintercept = 0, col = "black", linewidth = theme_cellpress_lwd) +
      ggrepel::geom_text_repel(data = subset(dmp_cat_plt, cluster_name == "C2" & median_stat_met >  1 & stat > 0 ), aes(label = value),
                              max.overlaps = Inf, size = theme_cellpress_size, col = "black") +
      theme_classic() +
      labs(x = "logFC / LFCse (methylation)", y = "logFC / LFCse (RNA)", title = "CATNON") +
      theme_cellpress +
      theme(legend.title = element_blank(),
            legend.position = "top")

dmp_tcga_plt <- res_ord_met_tcga %>% 
                dplyr::left_join(., annoEPIC, by = c("probeID" = "Name")) %>% 
                dplyr::left_join(., gene_to_probes, by = c("probeID" = "Name")) %>% 
                dplyr::filter(value %in% res_ord_tcga$gene_name) %>% 
                dplyr::group_by(value) %>% 
                dplyr::summarise(median_stat_met = median(stat_meth)) %>% 
                dplyr::left_join(., res_ord_tcga, by = c("value" = "gene_name")) %>% 
                dplyr::left_join(., clusters_catnon, by = c("value" = "gene_name"))  %>% 
                dplyr::mutate(cluster_name = dplyr::case_when(is.na(cluster) ~ NA,
                                                              cluster %in% c("up.4", "up.5") ~ "other",
                                                              cluster == "down.1" ~ "C0",
                                                              cluster == "cell_cyc_up" ~ "C1",
                                                              cluster == "emb_up" ~ "C2",
                                                              cluster == "col_up" ~ "C3"))
              

p4 <- ggplot(dmp_tcga_plt, aes(x=median_stat_met, y=stat, col = cluster_name)) +
      geom_point(data = subset(dmp_tcga_plt, is.na(cluster_name)), alpha = 0.5, col = "grey90") +
      geom_point(data = subset(dmp_tcga_plt, ! is.na(cluster_name)), alpha = 0.5) +
      scale_color_manual(values = c("#996632", "#00ff80",  "#ffd600", "red","grey90", "grey90")) +
      geom_vline(xintercept = 0, col = "black", linewidth = theme_cellpress_lwd) +
      geom_hline(yintercept = 0, col = "black", linewidth = theme_cellpress_lwd) +
      ggrepel::geom_text_repel(data = subset(dmp_tcga_plt, cluster_name == "C2" & median_stat_met >  1 & stat > 0 ), aes(label = value),
                              max.overlaps = Inf, size = theme_cellpress_size, col = "black") +
      theme_classic() +
      labs(x = "logFC / LFCse (methylation)", y = "", title = "TCGA") +
      theme_cellpress + 
      theme(legend.title = element_blank(),
            legend.position = "none")

p3+p4

# ggsave("output/figures/paper/rev1/f3/rna_vs_meth_lfc.pdf", width = (11.2 * 0.95) * 0.75, height = (8.5 * 0.95) * 0.4)
