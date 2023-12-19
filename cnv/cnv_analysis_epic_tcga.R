# load libs ----
library(tidyverse)
library(biovizBase)

# load data ----
source('scripts/R/load_chrom_sizes.R')
source("scripts/R/chrom_length_hg19.R")

data(hg19IdeogramCyto)
data(genesymbol)

hg19_bands <- as.data.frame(hg19IdeogramCyto) %>% 
  dplyr::mutate(start = start/10^6) %>% 
  dplyr::mutate(end = end/10^6) 

meta_tcga <- readRDS("data/clean/tcga/rds/metadata_tcga_idhmutnoncodel.Rds") %>% 
             dplyr::filter(predictBrain_12.8_cal_class %in% c('A_IDH_LG', 'A_IDH_HG'))

cnv.metadata <- readRDS("data/clean/tcga/rds/heidelberg_binvals.Rds") %>% 
                dplyr::filter(! cnv.segments.chrom %in% c("chrX", "chrY")) %>% 
                dplyr::select(cnv.segments.id, cnv.segments.chrom, cnv.segments.start, cnv.segments.end, cnv.segments.length) %>% 
                dplyr::mutate(cnv.segments.id = paste0('cnv.segment.',cnv.segments.id,':',cnv.segments.start,'-',cnv.segments.end)) 

cnv <- readRDS("data/clean/tcga/rds/heidelberg_binvals.Rds") %>% 
        dplyr::filter(! cnv.segments.chrom %in% c("chrX", "chrY")) %>% 
        dplyr::mutate(cnv.segments.id = paste0('cnv.segment.',cnv.segments.id,':',cnv.segments.start,'-',cnv.segments.end)) %>%
        dplyr::select(-cnv.segments.chrom, -cnv.segments.start, -cnv.segments.end, -cnv.segments.length) %>%
        tibble::column_to_rownames('cnv.segments.id') %>%
        t() %>%
        as.data.frame %>%
        tibble::rownames_to_column('BatchID') %>% 
        dplyr::mutate(BatchID = gsub("X", "", BatchID)) %>% 
        dplyr::filter(BatchID %in% meta_tcga$BatchID)

oncovals <- readRDS("data/clean/tcga/rds/heidelberg_genevals.Rds") 

# CNV plot Heidelberg methylation bins -----

# ** numeric -----
data <- cnv %>%
        tibble::column_to_rownames('BatchID')  

identical(rownames(data),meta_tcga$BatchID)

cnv.cor.lgc <- data.frame(cor.t.cnv.lgc = apply(data,2, function(vec) {return( cor.test(meta_tcga$log_AIDH_cal_ratio, 
                                                                                        as.numeric(vec), method = "spearman")$estimate) })) %>% 
  dplyr::mutate(cnv.segments.id = rownames(.)) %>% 
  dplyr::left_join(.,cnv.metadata )

plt.cnv.cor.lgc <- data.frame(chrs_hg19_s) %>% 
  tibble::rownames_to_column('chr') %>% 
  dplyr::filter(! chr %in% c('chrX', 'chrY')) %>% 
  dplyr::left_join(., cnv.cor.lgc, by = c("chr" = "cnv.segments.chrom")) %>% 
  dplyr::mutate(pos = chrs_hg19_s + cnv.segments.start) %>% 
  dplyr::mutate(chr = factor(chr, levels=gtools::mixedsort(unique(as.character(chr))) ))


ggplot(plt.cnv.cor.lgc, aes(x=pos, y=cor.t.cnv.lgc, col=chr)) +
  labs(x = "Position", y = "Correlation CNA locus with linear grading coefficient",x=NULL) +
  geom_point() +
  theme_classic() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = chrs_hg19_s) 

# ** discrete -----

data_amp <- cnv %>%
  tibble::column_to_rownames('BatchID')  %>%
  dplyr::mutate_all(~ case_when(. > 0.35 ~ "AMP",
                                TRUE ~ "NO-AMP")) %>% 
  tibble::rownames_to_column('BatchID')

data_del <- cnv %>%
  tibble::column_to_rownames('BatchID')  %>%
  dplyr::mutate_all(~ case_when(. < -0.415 ~ "DEL",
                                TRUE ~ "NO-DEL")) %>% 
  tibble::rownames_to_column('BatchID')

disc.cnv.out.del <- data_del %>% 
  tidyr::pivot_longer(!BatchID, names_to = "alteration", values_to = "status") %>% 
  dplyr::left_join(., meta_tcga) %>% 
  dplyr::select(log_AIDH_cal_ratio, alteration, status, BatchID) %>% 
  dplyr::mutate(group = ifelse(status == "NO-DEL", 1, 2)) %>% 
  dplyr::group_by(alteration) %>%  
  dplyr::mutate(n_levels = n_distinct(group)) %>% 
  dplyr::filter(n_levels > 1) %>% 
  dplyr::do(w = wilcox.test(log_AIDH_cal_ratio ~ group, data= .)) %>% 
  dplyr::summarise(alteration, wilc.pval = w$p.value) %>% 
  dplyr::mutate(wilc.padj = p.adjust(wilc.pval, method = "BH")) %>% 
  dplyr::mutate(p.sign = ifelse(wilc.padj < 0.05, T, F)) %>% 
  dplyr::left_join(.,cnv.metadata,by = c("alteration" = "cnv.segments.id") ) %>% 
  dplyr::mutate(log_padj = log10(wilc.padj) * 1) %>% 
  dplyr::mutate(test = "HD")

disc.cnv.out.amp <- data_amp %>% 
  tidyr::pivot_longer(!BatchID, names_to = "alteration", values_to = "status") %>% 
  dplyr::left_join(., meta_tcga) %>% 
  dplyr::select(log_AIDH_cal_ratio, alteration, status, BatchID) %>% 
  dplyr::mutate(group = ifelse(status == "NO-AMP", 1, 2)) %>% 
  dplyr::group_by(alteration) %>%  
  dplyr::mutate(n_levels = n_distinct(group)) %>% 
  dplyr::filter(n_levels > 1) %>% 
  dplyr::do(w = wilcox.test(log_AIDH_cal_ratio ~ group, data= .)) %>% 
  dplyr::summarise(alteration, wilc.pval = w$p.value) %>% 
  dplyr::mutate(wilc.padj = p.adjust(wilc.pval, method = "BH")) %>% 
  dplyr::mutate(p.sign = ifelse(wilc.padj < 0.05, T, F)) %>% 
  dplyr::left_join(.,cnv.metadata,by = c("alteration" = "cnv.segments.id") ) %>% 
  dplyr::mutate(log_padj = log10(wilc.padj) * -1) %>% 
  dplyr::mutate(test = "AMP")

disc.cnv.out <- rbind(disc.cnv.out.del, disc.cnv.out.amp)

plt.disc.cnv.out <- data.frame(chrs_hg19_s) %>% 
  tibble::rownames_to_column('chr') %>% 
  dplyr::filter(! chr %in% c('chrX', 'chrY')) %>% 
  dplyr::left_join(., disc.cnv.out, by = c("chr" = "cnv.segments.chrom")) %>% 
  dplyr::mutate(pos = chrs_hg19_s + cnv.segments.start) %>% 
  dplyr::mutate(chr = factor(chr, levels=gtools::mixedsort(unique(as.character(chr))) )) %>% 
  dplyr::mutate(start = cnv.segments.start/1000000) %>% 
  dplyr::mutate(end = cnv.segments.end/1000000) 

c <- plt.disc.cnv.out %>%
     dplyr::mutate(label = case_when(chr == "chr4" & start > 55.09 & end < 55.16 & wilc.padj < 0.05  & 
                                    test == "AMP" ~ "PDGFRA",
                                  chr == "chr4" & start > 57.49 & end < 57.56 & wilc.padj < 0.05 &
                                    test == "AMP" ~ "HOPX",
                                  chr == "chr12" & start > 58.09 & end < 58.16 & wilc.padj < 0.05 &
                                    test == "AMP" ~ "CDK4",
                                  chr == "chr9" & start > 21.40 & end < 22.01 & wilc.padj < 0.05 &
                                    test == "HD" ~ "CDKN2AB",
                                  chr == "chr9" & start > 2.00 & end < 2.06 & wilc.padj < 0.05 &
                                    test == "HD" ~ "SMARCA2",
                                  chr == "chr2" & start > 240.24 & end < 240.31 & wilc.padj < 0.05 &
                                    test == "HD" ~ "HDAC4",
                                  chr == "chr13" & start > 48.90 & end < 49.16 & wilc.padj < 0.05 &
                                    test == "HD" ~ "RB1",
                                  chr == "chr19" & start > 54.19 & end < 54.26 & wilc.padj < 0.05 &
                                    test == "HD" ~ "C19MC",
                                  chr == "chr10" & start > 89.65 & end < 89.76 & wilc.padj < 0.05 &
                                    test == "HD" ~ "PTEN")) %>% 
      dplyr::left_join(., chr_sizes) %>% 
      dplyr::mutate(Study = 'TCGA')

c.sig <- c %>% 
         dplyr::filter(p.sign == T) %>% 
         dplyr::group_by(chr) %>% 
         dplyr::mutate(freq = n()/length * 100) 

# saveRDS(c, 'output/tables/rds/cnv_tcga_tab.Rds')


ggplot(c, aes(x=pos, y=log_padj, label = label, col = p.sign)) +
       geom_hline(yintercept = log10(0.05)) +
       geom_hline(yintercept = -log10(0.05)) +
       facet_wrap(~ chr, nrow = 1) +
       labs(x = "Chromosomal position", 
            y = "Association CNA locus with linear grading coefficient",x=NULL,
            title = paste0("TCGA (n=", nrow(cnv), ")")) +
       geom_point() +
       theme_bw() +
       scale_color_manual(values = c("grey70", "#1F77B4")) +
       ggrepel::geom_text_repel(data = subset(c, log_padj > 0 & !is.na(label)),cex=2.8,
                           col="black",  direction = "y", nudge_y = 0.5, segment.size=0.25)  +
       ggrepel::geom_text_repel(data = subset(c, log_padj < 0 & !is.na(label)),cex=2.8,
                           col="black",  direction = "y", nudge_y = -0.5, segment.size=0.25) +
       theme(axis.text.x = element_blank()) 


ggsave("output/figures/paper/sf4/tcga_cnv_lgc.pdf", width = 17, height = 5)


subsetByOverlaps(hg19IdeogramCyto, genesymbol["SMARCA2",])

# Oncovals Heidelberg -----

# ** numeric ----
data_num <- meta_tcga %>% 
  dplyr::select(BatchID) %>% 
  dplyr::left_join(., oncovals, by = c("BatchID" = "BatchID")) %>% 
  tibble::column_to_rownames('BatchID')  

identical(rownames(data_num),meta_tcga$BatchID)

cor.onco.lgc <- data.frame(cor.t.onco.lgc = apply(data_num,2, function(vec) {return( cor.test(meta_tcga$log_AIDH_cal_ratio, 
                                                                                              as.numeric(vec), method = "spearman")$estimate) })) %>% 
  dplyr::mutate(cnv.segments.id = rownames(.)) %>% 
  dplyr::left_join(.,cnv.metadata )


# Surv ----

# cnv.segment.chr9-0120:21950001-22000000 - CDKN2AB -
# cnv.segment.chr13-0314:49050001-49100000 - RB1 -
# cnv.segment.chr12-0598:58100001-58150000 - CDK4 -
# cnv.segment.chr4-0410:55100001-55150000 - PDGFRA -
# cnv.segment.chr9-0016:2050001-2200000 - SMARCA2 -

surv_meta <- meta_tcga %>% 
  dplyr::left_join(., data_del) %>% 
  dplyr::mutate(`cnv.segment.chr9-0035:21400001-22000000` = factor(`cnv.segment.chr9-0035:21400001-22000000`, 
                                                                   levels = c("NO-DEL", "DEL"))) %>% 
  dplyr::mutate(`cnv.segment.chr13-0194:48900001-49150000` = factor(`cnv.segment.chr13-0194:48900001-49150000`, 
                                                                    levels = c("NO-DEL", "DEL"))) %>% 
  dplyr::mutate(`cnv.segment.chr9-0006:2000001-2050000` = factor(`cnv.segment.chr9-0006:2000001-2050000`, 
                                                                 levels = c("NO-DEL", "DEL")))

fit_cdkn2ab <- coxph(Surv(time = Survival, event = as.numeric(Deceased)) ~ `cnv.segment.chr9-0035:21400001-22000000`, data = surv_meta)
fit_rb1 <- coxph(Surv(time = Survival, event = as.numeric(Deceased)) ~ `cnv.segment.chr13-0194:48900001-49150000`, data = surv_meta)
fit_smarca2 <- coxph(Surv(time = Survival, event = as.numeric(Deceased)) ~ `cnv.segment.chr9-0006:2000001-2050000`, data = surv_meta)

summary(fit_cdkn2ab)$conf.int
summary(fit_rb1)$conf.int
summary(fit_smarca2)$conf.int

rm(surv_meta)

surv_meta <- meta_tcga %>% 
              dplyr::left_join(., data_amp) %>% 
              dplyr::mutate(`cnv.segment.chr12-0388:58100001-58150000` = factor(`cnv.segment.chr12-0388:58100001-58150000`, 
                                                                                levels = c("NO-AMP", "AMP"))) %>% 
              dplyr::mutate(`cnv.segment.chr4-0271:55100001-55150000` = factor(`cnv.segment.chr4-0271:55100001-55150000`, 
                                                                               levels = c("NO-AMP", "AMP")))

fit_cdk4 <- coxph(Surv(time = Survival, event = as.numeric(Deceased)) ~ `cnv.segment.chr12-0388:58100001-58150000`, data = surv_meta)
fit_pdgfra <- coxph(Surv(time = Survival, event = as.numeric(Deceased)) ~ `cnv.segment.chr4-0271:55100001-55150000`, data = surv_meta)

summary(fit_cdk4)$conf.int
summary(fit_pdgfra)$conf.int