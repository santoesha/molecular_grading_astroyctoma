library(msigdbr)

source("scripts/R/probe_to_gene.R")
source("scripts/clean/dge/catnon/load_hclust.R")
source("scripts/clean/R/theme_cellpress.R")


porath <- read.delim("data/pc2_targets_benporath.txt")

annoEPIC <- read.csv("data/methylation/annoEPIC.csv", row.names = c(1))  %>% 
            dplyr::filter(! chr %in% c("chrX", "chrY")) %>% 
            dplyr::filter(UCSC_RefGene_Name != "")

res_ord_cat <- readRDS("output/clean/catnon/rds/resord_catnon_heidi_num_lncrna.Rds") 
res_ord_tcga <- readRDS("output/clean/tcga/rds/resord_tcga_heidi_num_lncrna.Rds")

sig_met_cat <- readRDS("output/clean/catnon/rds/methylation_res_cat.Rds") %>% 
                       dplyr::mutate(significant = ifelse(abs(logFC) > 0.5 & adj.P.Val < 0.01, T, F)) %>% 
                       dplyr::mutate(sign = dplyr::case_when(significant == T & logFC < -0.5 ~ "hypomethylated",
                                                             significant == T & logFC > 0.5 ~ "hypermethylated",
                                                             T ~ NA)) %>% 
                       dplyr::filter(sign == 'hypermethylated')

sig_met_tcga <- readRDS("output/clean/tcga/rds/methylation_res_tcga.Rds") %>% 
                        dplyr::mutate(significant = ifelse(abs(logFC) > 0.5 & adj.P.Val < 0.01, T, F)) %>% 
                        dplyr::mutate(sign = dplyr::case_when(significant == T & logFC < -0.5 ~ "hypomethylated",
                                      significant == T & logFC > 0.5 ~ "hypermethylated",
                                      T ~ NA)) %>% 
                        dplyr::filter(sign == 'hypermethylated')

rna_to_met <- gene_to_probes %>% 
              #dplyr::filter(value %in% clusters_catnon$gene_name) %>% 
              dplyr::left_join(., annoEPIC %>% dplyr::select(Name, Relation_to_Island)) 

hyp_probes <- data.frame(probeID = intersect(sig_met_cat$probeID, sig_met_tcga$probeID)) %>% 
              dplyr::left_join(., rna_to_met, by = c("probeID" = "Name")) %>% 
              dplyr::filter(value %in% res_ord_cat$gene_name) %>% 
              dplyr::left_join(., clusters_catnon, by = c("value" = "gene_name")) %>%
              dplyr::mutate(porath_gene = ifelse(value %in% porath$prc2_target, T ,F)) %>% 
              distinct()
write.csv(hyp_probes, 'output/tables/rds/hyp_probes.csv')

length(unique(hyp_probes$value))

length(intersect(clusters_catnon %>% dplyr::filter(cluster != "down.1") %>% dplyr::pull(gene_name), hyp_probes$value))
length(intersect(clusters_catnon %>% dplyr::filter(cluster == "emb_up") %>% dplyr::pull(gene_name), hyp_probes$value))

k <- as.data.frame(unique(intersect(hyp_probes$value, porath$prc2_target)))

hyp_up <- intersect(clusters_catnon %>% dplyr::filter(cluster == "emb_up") %>% dplyr::pull(gene_name), hyp_probes$value)


msigdbr_c2 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP") %>% 
              dplyr::distinct(gs_name, gene_symbol) %>%
              as.data.frame()

msig_hyp <- enricher(gene = unique(hyp_probes$value),
                     TERM2GENE = msigdbr_c2, 
                     universe = res_ord_cat$gene_name) %>%
            as_tibble() %>% 
            slice_min(qvalue, n=25) %>% 
            dplyr::mutate(qvalue = log2(qvalue))


# Hypermethylated probes vs C2 -----

hyp_probes <- readRDS("output/clean/hypermethylated_probes.overlap.Rds")

mval_glass_hyper <- readRDS("data/clean/glass/rds/mvalues.nomaskedprobes.glass.Rds") %>% 
  dplyr::select(probeID, meta_glass$BatchID) %>% 
  dplyr::filter(probeID %in% hyp_probes) %>% 
  tibble::column_to_rownames('probeID') 





hyper_up_cat  <- sig_met_cat %>% 
                 dplyr::filter(sign == "hypermethylated") %>% 
                 dplyr::filter(probeID %in% hyp_probes) %>% 
                 dplyr::filter(probeID %in% rna_to_met$Name) %>% 
                 dplyr::left_join(rna_to_met, by = c("probeID" = "Name")) %>% 
                 dplyr::left_join(., res_ord_cat, by = c("value" = "gene_name")) %>% 
                 dplyr::left_join(., clusters_catnon, by = c("value" = "gene_name")) 
                 dplyr::filter(value %in% (clusters_catnon %>% 
                                            dplyr::filter(! cluster == "down.1") %>% 
                                            dplyr::pull(gene_name))) %>% 
                dplyr::mutate(Study = "CATNON") 

hyper_up_tcga  <- sig_met_tcga %>% 
                dplyr::filter(sign == "hypermethylated") %>% 
                dplyr::filter(probeID %in% sig_rna_to_met$Name) %>% 
                dplyr::left_join(sig_rna_to_met, by = c("probeID" = "Name")) %>% 
                dplyr::filter(value %in% (clusters_catnon %>% 
                                            dplyr::filter(! cluster == "down.1") %>% 
                                            dplyr::pull(gene_name))) %>% 
                dplyr::mutate(Study = "TCGA") 

a <- unique(intersect(hyper_up_cat$value, hyper_tcga$value))

length(unique(hyper_tcga$value))

length(intersect(hyper_cat$value, hyper_tcga$value))

tab_cat <- rbind(hyper_up_cat, hyper_up_tcga) %>% 
            dplyr::filter(value %in% intersect(hyper_up_cat$value, hyper_up_tcga$value)) %>%
            dplyr::mutate(label = paste0(probeID, " (", Relation_to_Island, ")")) %>% 
            dplyr::select(Study, logFC, log2FoldChange, value, label) %>% 
            dplyr::group_by(Study, value) %>% 
            dplyr::mutate(median_lfc_met = median(logFC)) 

# saveRDS(tab_cat, "output/clean/hypermeth_upreg_values.Rds")\



a <- enricher(gene = clusters_catnon %>% dplyr::filter(cluster == "emb_up") %>% dplyr::pull(gene_name),
              TERM2GENE = msigdbr_c2, 
              universe = res_ord_cat$gene_name) %>%
  as_tibble() %>% 
  slice_min(qvalue, n=50) %>% 
  dplyr::mutate(qvalue = log2(qvalue))
