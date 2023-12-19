# load libs ----
library(tidvyerse)
library(patchwork)
library(ggpubr)

# load data ----
source("scripts/clean/R/metadata_all.R")
source("scripts/clean/R/theme_cellpress.R")

cols <- c("A_IDH_HG" = "#E41A1C",
          "A_IDH_LG" = "#377EB8",
          "CPH_ADM" = "#4DAF4A",
          "CTRL_CBM" = "#984EA3",
          "CTRL_CORPCAL" = "#FF7F00",
          "CTRL_HEMI" = "#FFFF33",
          "CTRL_HYPOTHAL" = "#A65628",
          "CTRL_REACTIVE" = "black",
          "GBM_MES_TYP" = "#999999",
          "INFLAM_ENV" = "#66C2A5",
          "NB_TMM_NEG" = "#FC8D62",
          "O_IDH" = "#8DA0CB",
          "OLIGOSARC_IDH" = "#E78AC3",
          "pedHGG_RTK1A" = "#A6D854")

# LGC ----

# ** CATNON ----

lgc_cat <- ggplot(meta_catnon) +
                   geom_point(aes(x = reorder(BatchID, log_AIDH_cal_ratio),
                                  y = log_AIDH_cal_ratio,
                                  col = predictBrain_12.8_cal_class)) +
                   coord_equal() +
                   labs(y = "CGC") +
                   theme_cellpress +
                   theme(axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank(),
                         legend.title = element_blank(),
                         legend.position="right") +
                   scale_color_manual(values = cols) 

# ** TCGA ----

lgc_tcga <- ggplot(meta_tcga) +
                  geom_point(aes(x = reorder(Case, log_AIDH_cal_ratio),
                                 y = log_AIDH_cal_ratio,
                                 col = predictBrain_12.8_cal_class)) +
                  coord_equal() +
                  labs(y = "CGC") +
                  theme_cellpress +
                  theme(axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        legend.title = element_blank(),
                        legend.position="none") +
                  scale_color_manual(values = cols) +
                  scale_y_continuous(limits = c(-20, 20), breaks = c(-10, 0, 10)) 

# ** GLASS-P ----

lgc_glass_p <- ggplot(meta_glass_r1) +
                      geom_point(aes(x = reorder(BatchID, log_AIDH_cal_ratio),
                                     y = log_AIDH_cal_ratio,
                                     col = predictBrain_12.8_cal_class)) +
                      coord_equal() +
                      labs(y = "CGC") +
                      theme_cellpress +
                      theme(axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank(),
                            legend.title = element_blank(),
                            legend.position="none") +
                      scale_color_manual(values = cols) +
                      scale_y_continuous(limits = c(-20, 20), breaks = c(-10, 0, 10))

# ** GLASS-R ----

lgc_glass_r <- ggplot(meta_glass_r2) +
              geom_point(aes(x = reorder(BatchID, log_AIDH_cal_ratio),
                             y = log_AIDH_cal_ratio,
                             col = predictBrain_12.8_cal_class)) +
              coord_equal() +
              labs(y = "CGC") +
              theme_cellpress +
              theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    legend.title = element_blank(),
                    legend.position="none
                    ") +
              scale_color_manual(values = cols)  +
              scale_y_continuous(limits = c(-20, 20), breaks = c(-10, 0, 10))


# Surv ----

# ** CATNON ----

surv_cat <- ggplot(meta_catnon %>% dplyr::mutate(Deceased = ifelse(Deceased == 0, "Alive/Censored", "Deceased")) , 
                   aes(x=reorder(BatchID,log_AIDH_cal_ratio),
                       y=Survival, 
                       fill=Deceased)) +
            geom_bar(position = "stack", stat='identity') +
            labs(y = "Survival (years)") +
            theme_cellpress +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  legend.title = element_blank(),
                  legend.position="none") +
            scale_fill_brewer(palette = "Reds") +
            coord_equal() +
            scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5))


# ** TCGA ----

surv_tcga <- ggplot(meta_tcga %>% dplyr::mutate(Deceased = ifelse(Deceased == 0, "Alive/Censored", "Deceased")) , 
                    aes(x=reorder(Case,log_AIDH_cal_ratio),
                        y=Survival, 
                        fill=Deceased)) +
            geom_bar(position = "stack", stat='identity') +
            labs(y = "Survival (years)") +
            theme_cellpress +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  legend.title = element_blank(),
                  legend.position="none") +
            scale_fill_brewer(palette = "Reds") +
            coord_equal() +
            scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5))

# ** GLASS-P ----

surv_glass_p <- ggplot(meta_glass_r1 %>% dplyr::mutate(Deceased = ifelse(Deceased == 0, "Alive/Censored", "Deceased")) , 
                       aes(x=reorder(BatchID,log_AIDH_cal_ratio),
                           y=Survival, 
                           fill=Deceased)) +
                geom_bar(position = "stack", stat='identity') +
                labs(y = "Survival (years)") +
                theme_cellpress +
                theme(axis.title.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.x=element_blank(),
                      legend.title = element_blank(),
                      legend.position="none") +
                scale_fill_brewer(palette = "Reds") +
                coord_equal() +
                scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 5))

# Mut ----

# ** CATNON ----

ord <- meta_catnon %>%
       dplyr::select(BatchID,log_AIDH_cal_ratio)

dat_mut <- meta_catnon %>%
          dplyr::select(BatchID, CDKN2AB_del, CDK4_amp, PDGFRA_amp, RB1_del, CDK6_amp) %>% 
          reshape2::melt(., "BatchID") %>% 
          dplyr::left_join(., ord) 

mut_cat <- ggplot(dat_mut, aes(x = reorder(BatchID,log_AIDH_cal_ratio), y = variable, fill = value)) +
            geom_tile(colour = "grey70", size = 0.05) +
            theme(axis.title.x=element_blank(),
                  axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1),
                  axis.title.y=element_blank(),
                  axis.ticks.y =element_blank(),
                  legend.title = element_blank(),
                  legend.position="right") +
            scale_fill_manual(values = c("white","seagreen")) +
            coord_equal() +
            theme_cellpress +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  legend.title = element_blank()) +
            labs(y = "")

rm(dat_mut, ord)


# ** TCGA -----


ord <- meta_tcga %>%
      dplyr::select(Case,log_AIDH_cal_ratio)

dat_mut <- meta_tcga %>%
          dplyr::select(Case, CDKN2AB_del, CDK4_amp, PDGFRA_amp, RB1_del, CDK6_amp) %>% 
          reshape2::melt(., "Case") %>% 
          dplyr::left_join(., ord) 

mut_tcga <- ggplot(dat_mut, aes(x = reorder(Case,log_AIDH_cal_ratio), y = variable, fill = value)) +
            geom_tile(colour = "grey70", size = 0.05) +
            scale_fill_manual(values = c("white","seagreen")) +
            coord_equal() +
            theme_cellpress +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  legend.title = element_blank(),
                  legend.position = "none") +
            labs(y = "")

rm(dat_mut, ord)


# ** GLASS-P ----

ord <- meta_glass_r1 %>%
       dplyr::select(BatchID,log_AIDH_cal_ratio)

dat_mut <- meta_glass_r1 %>%
            dplyr::select(BatchID, CDKN2AB_del, CDK4_amp, PDGFRA_amp, RB1_del, CDK6_amp) %>% 
            reshape2::melt(., "BatchID") %>% 
            dplyr::left_join(., ord) 

mut_glass_p <- ggplot(dat_mut, aes(x = reorder(BatchID,log_AIDH_cal_ratio), y = variable, fill = value)) +
                      geom_tile(colour = "grey70", size = 0.05) +
                      scale_fill_manual(values = c("white","seagreen")) +
                      coord_equal() +
                      theme_cellpress +
                      theme(axis.title.x=element_blank(),
                            axis.text.x = element_blank(),
                            axis.title.y=element_blank(),
                            axis.ticks.y =element_blank(),
                            legend.title = element_blank(),
                            legend.position="none")


rm(dat_mut, ord)

# ** GLASS-R ----

ord <- meta_glass_r2 %>%
        dplyr::select(BatchID,log_AIDH_cal_ratio)

dat_mut <- meta_glass_r2 %>%
            dplyr::select(BatchID, CDKN2AB_del, CDK4_amp, PDGFRA_amp, RB1_del, CDK6_amp) %>% 
            reshape2::melt(., "BatchID") %>% 
            dplyr::left_join(., ord) %>% 
            dplyr::mutate(value = ifelse(value %in% c(TRUE, "IDH1_R132H"), T, F))

mut_glass_r <- ggplot(dat_mut, aes(x = reorder(BatchID,log_AIDH_cal_ratio), y = variable, fill = value)) +
                  geom_tile(colour = "grey70", size = 0.05) +
                  scale_fill_manual(values = c("white","seagreen")) +
                  coord_equal() +
                  theme_cellpress +
                  theme(axis.title.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        legend.title = element_blank(),
                        legend.position = "none") +
                  labs(y = "")


# Resection GLASS -----
res_glass <- readRDS("data/clean/glass/rds/metadata.glass.Rds") %>%
              dplyr::select(BatchID, Resectie, log_AIDH_cal_ratio) %>%
              dplyr::mutate(r2 = ifelse(Resectie == 2, T, F)) %>%
              dplyr::mutate(r3 = ifelse(Resectie == 3, T, F)) %>% 
              dplyr::mutate(r4 = ifelse(Resectie == 4, T, F)) %>% 
              dplyr::filter(Resectie != 1) %>% 
              dplyr::select(BatchID, r2, r3, r4, log_AIDH_cal_ratio) %>% 
              tidyr::pivot_longer(cols = starts_with("r"))

res_glass_r <- ggplot(res_glass, aes(x = reorder(BatchID,log_AIDH_cal_ratio), y = name, fill = value)) +
                      geom_tile(colour = "grey70", size = 0.05) +
                      scale_fill_manual(values = c("white","seagreen")) +
                      coord_equal() +
                      theme_cellpress +
                      theme(axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank(),
                            legend.title = element_blank(),
                            legend.position="none") +
                      labs(y = "")


rm(dat_mut, ord)

# Plot ----
p1 <- mut_cat/lgc_cat/surv_cat  + plot_layout(guides = "collect") & theme(legend.position = 'top')
p1
# ggsave("output/figures/paper/rev1/f1/lgc_f1_catnon.pdf", width = 11.2, height = (8.5 * 0.95) * (1/3))

p2 <- (mut_tcga/lgc_tcga/surv_tcga) +
      theme_cellpress & 
      theme(legend.position = 'none',
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            legend.title = element_blank())
p2
# ggsave("output/figures/paper/rev1/f1/lgc_f1_tcga.pdf", width = 11.2 * 1/3, height = (8.5 * 0.95) * (0.4))

p3 <- (mut_glass_p/lgc_glass_p/surv_glass_p) +
      theme_cellpress & 
      theme(legend.position = 'none',
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            legend.title = element_blank())


mut_glass_p
ggsave("output/figures/paper/rev1/f1/mut_f1_glassp.pdf", width = 11.2 * 1/2, height = (6 * 0.95) * (1/3))

ggsave("output/figures/paper/rev1/f1/lgc_f1_glass.pdf", width = 11.2 * 1/2, height = (6 * 0.95) * (1/3))

p4 <- mut_glass_r/lgc_glass_r/res_glass_r +
      theme_cellpress & 
      theme(legend.position = 'none',
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            legend.title = element_blank())
p4

ggsave("output/figures/paper/rev1/f1/lgc_f1_glass_r2.pdf", width = 11.2 * 1/2, height = (6 * 0.95) * (1/3))


ggsave("output/figures/paper/rev1/f1/lgc_f1_glasstcga.pdf", width = 11.2, height = (8.5 * 0.95) * (1/3))

ggarrange(p3, p4, nrow = 1, align = "h")
ggsave("output/figures/paper/rev1/f1/lgc_f1_glasstcga.pdf", width = 5, height = (8.5 * 0.95) * (1/3))

ggarrange(p3, p4, ncol = 3,align = 'hv')

# ggsave("output/figures/paper/rev1/f1/lgc_f1_tcgaglass.pdf", units = "in", expand = c(0.1, 0.1))


a = cowplot::plot_grid(mut_tcga, lgc_tcga, surv_tcga, ncol = 1,  align = 'v', axis = 'lr', rel_heights = c(3,2,2))
b = cowplot::plot_grid(mut_glass_p, lgc_glass_p, surv_glass_p, ncol = 1,  align = 'v', axis = 'lr', rel_heights = c(3,2,2))
c = cowplot::plot_grid(mut_glass_r, lgc_glass_r, res_glass_r, ncol = 1,  align = 'v', axis = 'lr', rel_heights = c(3,2,2))

cowplot::plot_grid(a,b,c, ncol = 3,  align = 'hv', axis = 'tblr')

cowplot::plot_grid(mut_tcga,  mut_glass_p,  mut_glass_r, 
                   lgc_tcga,  lgc_glass_p,  lgc_glass_r,
                   surv_tcga, surv_glass_p, res_glass_r, 
                   ncol = 3,  align = 'hv', axis = 'tblr')





