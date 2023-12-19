library(patchwork)
library(RColorBrewer)

source("scripts/clean/R/theme_cellpress.R")
source("scrips/R/palette.R")

# load data -----
meta_catnon <- readRDS("data/clean/catnon/rds/metadata.catnon.Rds")
meta_tcga <- readRDS("data/clean/tcga/rds/metadata_tcga_idhmutnoncodel.Rds")
meta_glass_r1 <- readRDS("data/clean/glass/rds/metadata.glass.Rds") %>% 
                 dplyr::filter(Resectie == 1)
meta_glass_rec <- readRDS("data/clean/glass/rds/metadata.glass.Rds") %>% 
                  dplyr::filter(Resectie != 1) %>% 
                  dplyr::group_by(GLASS_ID) %>% 
                  dplyr::mutate(rec_samp = n()) 

# grade -----

catnon.tmp <- meta_catnon %>% 
              dplyr::mutate(who_2016 = ifelse(`Necrosis and/or microvascular proliferation` == "Present", "G4", "G2/G3"))

glass.tmp <- read.csv("~/mnt/neuro-genomic-ro/glass/Metadata/Samples/WHOclassification_03052022.csv") %>% 
                      dplyr::filter(GLASS_ID %in% meta_glass_r1$GLASS_ID) %>%
                      dplyr::filter(Sample_Type == "I") %>% 
                      dplyr::mutate(who_2016 = ifelse(WHO_Classification2016 == "Glioblastoma, IDH-mutant, WHO grade IV", "G4", "G2/G3")) 

tcga.tmp <-  xlsx::read.xlsx('data/TCGA methylation subtypes.xlsx', sheetName = 'S1A. TCGA discovery dataset', startRow = 2) %>% 
             dplyr::filter(Case %in% meta_tcga$Case) %>% 
             dplyr::mutate(Grade = ifelse(Grade == "NA", NA, Grade)) %>% 
             dplyr::mutate(who_2016 = ifelse(Grade == "G4", "G4", "G2/G3"))

chisq.test(glass.tmp$who_2016, catnon.tmp$who_2016)

dat <- matrix(c(table(catnon.tmp$who_2016), table(glass.tmp$who_2016)), nrow = 2)

fisher.test(dat)

# age -----

# CATNON
round(median(meta_catnon$age), 0)
round(min(meta_catnon$age), 0)
round(max(meta_catnon$age), 0)

# TCGA
round(median(meta_tcga$age_at_diagnosis, na.rm = T), 0)
round(min(meta_tcga$age_at_diagnosis, na.rm = T), 0)
round(max(meta_tcga$age_at_diagnosis, na.rm = T), 0)

# GLASS
round(median(meta_glass_r1$age_at_diagnosis, na.rm = T), 0)
round(min(meta_glass_r1$age_at_diagnosis, na.rm = T), 0)
round(max(meta_glass_r1$age_at_diagnosis, na.rm = T), 0)

wilcox.test(meta_glass_r1$age_at_diagnosis, meta_catnon$age)
wilcox.test(meta_glass_r1$age_at_diagnosis, meta_tcga$age)

# WHO ----

tmp.cat <- meta_catnon %>% dplyr::filter(! is.na(who_2021)) 
table(tmp.cat$who_2021)

tmp1.glass <- meta_glass_r1 %>% dplyr::filter(! is.na(WHO_Classification2021)) 
table(tmp1.glass$WHO_Classification2021)

tmp2.glass <- meta_glass_rec %>% dplyr::filter(! is.na(WHO_Classification2021)) 
table(tmp2.glass$WHO_Classification2021)

rm(tmp.cat, tmp1.glass, tmp2.glass)

# Methylation diagnosis -----

tmp <- meta_catnon %>% 
       dplyr::filter(! predictBrain_12.8_cal_class %in% c("A_IDH_LG", "A_IDH_HG")) %>% 
       dplyr::select(who_2021, predictBrain_12.8_cal_class, BatchID, 
                     `G-CIMP subtype`, predictBrain_12.8_cal_O_IDH, predictBrain_12.8_cal_OLIGOSARC_IDH,
                     predictBrain_12.8_cal_pedHGG_RTK1A, predictBrain_12.8_cal_INFLAM_ENV, predictBrain_12.8_cal_A_IDH_HG,
                     predictBrain_12.8_cal_A_IDH_LG) %>% 
       rowwise() %>% 
       dplyr::mutate(max_score = base::max(predictBrain_12.8_cal_O_IDH, predictBrain_12.8_cal_OLIGOSARC_IDH,
                                           predictBrain_12.8_cal_pedHGG_RTK1A, predictBrain_12.8_cal_INFLAM_ENV)) 

calc_ratio_detP <- t(detP) %>%
                   as.data.frame(stringsAsFactors=F) %>%
                   dplyr::mutate(failed_probes = rowSums(. > 0.01)/ncol(.)) %>%
                   dplyr::mutate(no_astro = ifelse(rownames(.) %in% tmp$BatchID, T , F)) %>%
                   dplyr::filter(failed_probes > 0.05)


calc_ratio_detP <- detP %>%
                   as.data.frame(stringsAsFactors=F) %>%
                   dplyr::mutate(failed_probes = rowSums(. > 0.01)/ncol(.)) %>%
  dplyr::mutate(no_astro = ifelse(rownames(.) %in% tmp$BatchID, T , F)) %>%
  dplyr::filter(failed_probes > 0.025)


ggplot(calc_ratio_detP, aes(x = reorder(rownames(calc_ratio_detP),failed_probes), y = failed_probes)) +
      geom_point()







meth_samp_cat <- readRDS("data/clean/catnon/rds/methylation.targets.catnon.Rds")


w <- function(fn, prefix) {
  
  a = read.csv(fn) %>%
    tibble::column_to_rownames('X') %>%
    `colnames<-`('pval') %>%
    dplyr::arrange(-pval) 
  
  top <- a %>%
    tibble::rownames_to_column('class') %>%
    dplyr::slice_head(n=1) %>%
    dplyr::pull(class)
  
  a<- a %>%
    t() %>%
    as.data.frame() %>%
    dplyr::mutate(class = top) %>%
    dplyr::rename_with( ~ paste0(prefix,"cal_", .x))
  
  return(a)
}


catnon <- data.frame(predictBrain.scores.file = Sys.glob("tmp/cache/CATNON/v12.8/*/predictBrain_v12.8/*_scores_cal.csv")) %>%
          dplyr::mutate(BatchID = gsub("^.+_v12.8/([^/]+)_scores_cal.csv$","\\1", predictBrain.scores.file))  %>%
          dplyr::filter(BatchID %in% meth_samp_cat$BatchID) %>% 
          dplyr::mutate(version = gsub("^.+predictBrain_([^_\\/]+)[_/].+$","\\1", predictBrain.scores.file)) %>%
          tidyr::pivot_wider(id_cols = BatchID, names_from = version, values_from = c(predictBrain.scores.file),
                             names_prefix = "heidelberg_reportBrain_") %>%
          dplyr::rowwise() %>%
          dplyr::mutate(tmp = w(heidelberg_reportBrain_v12.8, "predictBrain_12.8_")) %>%
          dplyr::ungroup() %>%
          tidyr::unnest(tmp) %>%
          dplyr::select(BatchID, predictBrain_12.8_cal_class, predictBrain_12.8_cal_O_IDH,
                        predictBrain_12.8_cal_A_IDH_LG, predictBrain_12.8_cal_A_IDH_HG,
                        predictBrain_12.8_cal_CTRL_HEMI, predictBrain_12.8_cal_OLIGOSARC_IDH,
                        predictBrain_12.8_cal_INFLAM_ENV, predictBrain_12.8_cal_pedHGG_RTK1A,
                        predictBrain_12.8_cal_GBM_RTK2, predictBrain_12.8_cal_EPN_ST_ZFTA_RELA_A,
                        predictBrain_12.8_cal_DLBCL, predictBrain_12.8_cal_CPH_ADM) %>% 
          dplyr::mutate(Study = "CATNON")

glass <- data.frame(predictBrain.scores.file = Sys.glob("~/mnt/neuro-genomic-ro/glass/Methylation/Heidelberg/brain_classifier_v12.8_sample_report__v1.1__131/*/predictBrain_v12.8/*_scores_cal.csv")) %>%
         dplyr::mutate(BatchID = gsub("^.+_v12.8/([^/]+)_scores_cal.csv$","\\1", predictBrain.scores.file))  %>%
         dplyr::mutate(version = gsub("^.+predictBrain_([^_\\/]+)[_/].+$","\\1", predictBrain.scores.file)) %>%
         tidyr::pivot_wider(id_cols = BatchID, names_from = version, values_from = c(predictBrain.scores.file),
                           names_prefix = "heidelberg_reportBrain_") %>%
         dplyr::rowwise() %>%
         dplyr::mutate(tmp = w(heidelberg_reportBrain_v12.8, "predictBrain_12.8_")) %>%
         dplyr::ungroup() %>%
         tidyr::unnest(tmp) %>%
         dplyr::select(BatchID, predictBrain_12.8_cal_class, predictBrain_12.8_cal_O_IDH,
                      predictBrain_12.8_cal_A_IDH_LG, predictBrain_12.8_cal_A_IDH_HG,
                      predictBrain_12.8_cal_CTRL_HEMI, predictBrain_12.8_cal_OLIGOSARC_IDH,
                      predictBrain_12.8_cal_INFLAM_ENV, predictBrain_12.8_cal_pedHGG_RTK1A,
                      predictBrain_12.8_cal_GBM_RTK2, predictBrain_12.8_cal_EPN_ST_ZFTA_RELA_A,
                      predictBrain_12.8_cal_DLBCL, predictBrain_12.8_cal_CPH_ADM) %>% 
         dplyr::mutate(Study = "GLASS-NL")

tcga <- data.frame(predictBrain.scores.file = Sys.glob("tmp/cache/TCGA/v12.8/*/predictBrain_v12.8/*_scores_cal.csv")) %>%
        dplyr::mutate(Case = gsub("\\__executed.*","\\1", predictBrain.scores.file))  %>% 
        dplyr::mutate(Case = gsub(".*(TCGA)", "\\1", Case)) %>% 
        dplyr::mutate(version = gsub("^.+predictBrain_([^_\\/]+)[_/].+$","\\1", predictBrain.scores.file)) %>%
        tidyr::pivot_wider(id_cols = Case, names_from = version, values_from = c(predictBrain.scores.file),
                           names_prefix = "heidelberg_reportBrain_") %>%
        dplyr::rowwise() %>%
        dplyr::mutate(tmp = w(heidelberg_reportBrain_v12.8, "predictBrain_12.8_")) %>%
        dplyr::ungroup() %>%
        tidyr::unnest(tmp) %>%
        dplyr::rename(BatchID = Case) %>% 
        dplyr::select(BatchID, predictBrain_12.8_cal_class, predictBrain_12.8_cal_O_IDH,
                      predictBrain_12.8_cal_A_IDH_LG, predictBrain_12.8_cal_A_IDH_HG,
                      predictBrain_12.8_cal_CTRL_HEMI, predictBrain_12.8_cal_OLIGOSARC_IDH,
                      predictBrain_12.8_cal_INFLAM_ENV, predictBrain_12.8_cal_pedHGG_RTK1A,
                      predictBrain_12.8_cal_GBM_RTK2, predictBrain_12.8_cal_EPN_ST_ZFTA_RELA_A,
                      predictBrain_12.8_cal_DLBCL, predictBrain_12.8_cal_CPH_ADM) %>% 
        dplyr::mutate(Study = "TCGA")

class_heidi <- rbind(catnon, glass,tcga) %>% 
               dplyr::mutate(log_AIDH_cal_ratio = log(predictBrain_12.8_cal_A_IDH_HG/predictBrain_12.8_cal_A_IDH_LG)) 

ord <- class_heidi %>%
       dplyr::select(BatchID,log_AIDH_cal_ratio, Study, 
                     predictBrain_12.8_cal_A_IDH_HG, predictBrain_12.8_cal_A_IDH_LG,
                     predictBrain_12.8_cal_O_IDH, predictBrain_12.8_cal_OLIGOSARC_IDH)

dat_mut <- class_heidi %>%
           dplyr::select(BatchID,predictBrain_12.8_cal_class) %>% 
           reshape2::melt(., "BatchID") %>% 
           dplyr::left_join(., ord) %>% 
           dplyr::mutate(astro = ifelse(value %in% c("A_IDH_LG", "A_IDH_HG"), "IDHmut-astro", "Other diagnosis")) %>% 
           dplyr::group_by(Study, astro) %>% 
           dplyr::mutate(count = n()) %>% 
           dplyr::ungroup() %>% 
           dplyr::mutate(label = paste0(astro, " (n=", count, ")")) %>% 
           dplyr::mutate(other_cal = 1 - predictBrain_12.8_cal_A_IDH_HG - predictBrain_12.8_cal_A_IDH_LG -
                                         predictBrain_12.8_cal_O_IDH - predictBrain_12.8_cal_OLIGOSARC_IDH)

a <- dat_mut %>% 
     dplyr::filter(astro == "Other diagnosis") %>% 
     dplyr::rename(diagnosis = value) %>% 
     dplyr::select(predictBrain_12.8_cal_A_IDH_HG, predictBrain_12.8_cal_A_IDH_LG,
                   predictBrain_12.8_cal_O_IDH, predictBrain_12.8_cal_OLIGOSARC_IDH,
                   other_cal, diagnosis, Study, BatchID) %>% 
     tidyr::pivot_longer(cols = c("predictBrain_12.8_cal_A_IDH_HG", 
                                  "predictBrain_12.8_cal_A_IDH_LG",
                                  "predictBrain_12.8_cal_O_IDH", 
                                  "predictBrain_12.8_cal_OLIGOSARC_IDH", 
                                  "other_cal")) %>% 
     dplyr::mutate(label = ifelse(! diagnosis %in% c("OLIGOSARC_IDH", "O_IDH"), "Other", diagnosis))

 
p1 <- ggplot(dat_mut, aes(x = reorder(BatchID, value), y = variable, fill = value)) +
      facet_grid(Study ~ astro) +
      geom_tile(size = 3) +
      theme_cellpress +
      theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y =element_blank(),
        axis.ticks.x = element_blank()) +
      scale_fill_manual(values = diag_col) 

ggsave("output/figures/paper/sf5/classes_heidi.pdf", width = 11, height = 8.5*1/3)

p2 <- ggplot(data = a , aes(x = BatchID, y = value, fill = name)) +
           facet_grid(Study ~ label) +
           geom_bar(position="stack", stat="identity") +
           theme_cellpress +
           theme(axis.title.x=element_blank(),
                 axis.text.x = element_blank(),
                 axis.title.y=element_blank(),
                 axis.ticks.y =element_blank(),
                 axis.ticks.x = element_blank()) +
           scale_fill_manual(values = prob_col)
