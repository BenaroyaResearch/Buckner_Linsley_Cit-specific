#### Set up #### 
setwd("~/Box/RA Cliff Metaclustering")
library(tidyverse)
library(openxlsx)
library(geepack)
library(gridExtra)
library(grid)
library(broom)
library(visreg)

# Set up the ggplot default params
theme_set(theme_classic(10) + 
            theme(axis.title = element_text(size=10, face="bold"),
                  axis.text.x = element_text(size=10, face="bold", color = "black"),
                  axis.text.y = element_text(size=8, color = "black")))

#### Load data #### 
# mc_data_path = file.path("Output", "190320", "heatmaps", 
#                          "190320metacluster_counts_for_modeling.csv")
# 
# # Clean up data
# mc_count_data = read_csv(mc_data_path)
# 
# das28_relabel_vector = c(`Near remission` = "< 2.6",
#                          `Low` = "2.6 \u2264 3.2",
#                          `Moderate` = "3.2  \u2264 5.1",
#                          `High` = "> 5.1")
# 
# mc_count_data = mc_count_data %>%
#   mutate(disease_status = factor(disease_status, levels = c("HC", "RA")),
#          # duration = factor(duration, levels = c("Recent onset", "Long standing")),
#          duration = case_when(duration == "Recent onset" ~ "\u2264 5",
#                               duration == "Long standing" ~ ">5") %>%
#            factor(levels = c("\u2264 5", ">5")),
#          # rapid3_cat = factor(rapid3_cat, levels = c("Near remission", "Low", "Moderate", "High"), ordered = T),
#          rapid3_cat =  case_when(rapid3_cat == "Near remission" ~ "0 - 1",
#                                  rapid3_cat == "Low" ~ "1.3 - 2",
#                                  rapid3_cat == "Moderate" ~ "2.3 - 4",
#                                  rapid3_cat == "High" ~ "4.3 - 10",
#                                  is.na(rapid3_cat) ~ NA_character_) %>%
#            factor(levels = c("0 - 1", "1.3 - 2", "2.3 - 4", "4.3 - 10"), ordered = T),
#          das28_cat = factor(das28_cat, levels = c("Near remission", "Low", "Moderate", "High"), ordered = T) %>%
#            recode_factor(!!!das28_relabel_vector),
#          das28_esr4_cat = factor(das28_esr4_cat, levels = c("Near remission", "Low", "Moderate", "High"), ordered = T) %>%
#            recode_factor(!!!das28_relabel_vector),
#          das28_crp3_cat = factor(das28_crp3_cat, levels = c("Near remission", "Low", "Moderate", "High"), ordered = T) %>%
#            recode_factor(!!!das28_relabel_vector),
#          das28_crp4_cat = factor(das28_crp4_cat, levels = c("Near remission", "Low", "Moderate", "High"), ordered = T) %>%
#            recode_factor(!!!das28_relabel_vector),
#          crp_cat = factor(crp_cat, levels = c("Normal", "High")),
#          TNFi = ifelse(medClass == "TNFblock", "TNFi",
#                        ifelse(medClass == "DMARD", "DMARD", NA))  %>%
#            factor(levels = c("DMARD", "TNFi")),
#          age_cat = factor(age_cat, levels = c("Young", "Old")),
#          sex = factor(sex),
#          Group = factor(Group, levels = c("HC", "RO", "LS"))) %>%
#   dplyr::rename(duration_cat = duration, duration = duration_yrs, group = Group)

mc_data_path = file.path("Output", "191211", "heatmaps", 
                         "191211metacluster_counts_for_modeling.csv")
mc_count_data = read_csv(mc_data_path)

load("20191210_combined_phenos.RData")
phenos <- abbrev_phenos

mc_count_data = mc_count_data %>%
  left_join(phenos, by = c("subject" = "subj")) %>%
  mutate(rapid3_cat =  case_when(rapid3_cat == "Near remission" ~ "0 - 1",
                                 rapid3_cat == "Low" ~ "1.3 - 2",
                                 rapid3_cat == "Moderate" ~ "2.3 - 4",
                                 rapid3_cat == "High" ~ "4.3 - 10",
                                 is.na(rapid3_cat) ~ NA_character_) %>%
           factor(levels = c("0 - 1", "1.3 - 2", "2.3 - 4", "4.3 - 10"), ordered = T),
         duration_cat = case_when(duration_cat == "Recent onset" ~ "\u2264 5",
                                  duration_cat == "Long standing" ~ ">5") %>%
           factor(levels = c("\u2264 5", ">5")))


#### Age & sex-controlled models - Wald Test ####
## CD4 in MC3 with Disease
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 3, !is.na(disease_status), Influenza_tot >= 1)

total_cell_model = glm(Cells_in_MC/Cells_tot ~ 1 + disease_status + age + sex,
                       weights = Cells_tot,
                       family = quasibinomial(),
                       data = analysis_data) 

fig2A_pt1 = visreg(total_cell_model, 
                   xvar = "disease_status", 
                   ylim = range(analysis_data$Cells_in_MC/analysis_data$Cells_tot), 
                   rug = F, 
                   line=list(col="red"),
                   scale = "response",
                   gg = T) +
  geom_point(data = analysis_data, 
             aes(x = ifelse(disease_status == "HC", 0.225, 0.775), 
                 y = Cells_in_MC/Cells_tot, 
                 shape = disease_status), 
             size = 1.5, position = position_jitter(width = 0.1)) +
  # annotate('text', -Inf, Inf,
  #          label=paste("italic('p') == 0.0399"), parse=TRUE,
  #          hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% CD4+ cells in AC3", x = "") + 
  scale_shape_manual(values = setNames(c(19, 1), c("HC", "RA")), guide = F) +
  scale_y_continuous(labels = function(x) (x*100))

## CD4 in MC2 with Rapid3cat
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 2, !is.na(rapid3_cat), Influenza_tot >= 1)

total_cell_model = glm(Cells_in_MC/Cells_tot ~ 1 + rapid3_cat + age + sex,
                       weights = Cells_tot,
                       family = quasibinomial(),
                       data = analysis_data) 

fig2A_pt2 = visreg(total_cell_model, 
                   xvar = "rapid3_cat", 
                   ylim = range(analysis_data$Cells_in_MC/analysis_data$Cells_tot), 
                   rug = F, 
                   line=list(col="red"),
                   scale = "response",
                   gg = T) +
  geom_point(data = analysis_data, 
             aes(x = ifelse(rapid3_cat == "0 - 1", 0.11,
                            ifelse(rapid3_cat == "1.3 - 2", 0.36,
                                   ifelse(rapid3_cat == "2.3 - 4", 0.63, 0.89))), 
                 y =Cells_in_MC/Cells_tot), 
             size = 1.5, shape = 1, position = position_jitter(width = 0.07)) +
  # annotate('text', -Inf, Inf,
  #          label=paste("italic('p') == 0.0005"), parse=TRUE,
  #          hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% CD4+ cells in AC2", x = "Disease Activity Score\n(Rapid3)") + 
  scale_y_continuous(labels = function(x) (x*100))

# Output Fig 2A
fig2A = grid.draw(rbind(ggplotGrob(fig2A_pt1),
              ggplotGrob(fig2A_pt2),
              size = "last"))


## Agg in MC5 with Duration
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 5, !is.na(duration_cat), Aggrecan_tot >= 8)

antigenSpec_model = glm(Aggrecan_in_MC/Aggrecan_tot ~ 1 + duration_cat + age + sex,
                        weights = Aggrecan_tot,
                        family = quasibinomial(),
                        data = analysis_data) 

fig2B_pt1 = visreg(antigenSpec_model, 
                   xvar = "duration_cat", 
                   ylim = range(analysis_data$Aggrecan_in_MC/analysis_data$Aggrecan_tot), 
                   rug = F, 
                   line=list(col="red"),
                   scale = "response",
                   gg = T) +
  geom_point(data = analysis_data, 
             aes(x = ifelse(duration_cat == "\u2264 5", 0.225, 0.775), 
                 y = Aggrecan_in_MC/Aggrecan_tot), 
             size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
  # annotate('text', -Inf, Inf,
  #          label=paste("italic('p') == 0.0026"), parse=TRUE,
  #          hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% Aggrecan-specific\nT cells in AC5", x = "") + 
  scale_y_continuous(labels = function(x) (x*100))

## Cilp in MC1 with Duration
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 1, !is.na(duration_cont), Cilp_tot >= 8)

antigenSpec_model = glm(Cilp_in_MC/Cilp_tot ~ 1 + duration_cont + age + sex,
                        weights = Cilp_tot,
                        family = quasibinomial(),
                        data = analysis_data)

fig2B_pt2 = visreg(antigenSpec_model, 
                   xvar = "duration_cont", 
                   ylim = range(analysis_data$Cilp_in_MC/analysis_data$Cilp_tot), 
                   rug = F, 
                   line=list(col="red"), 
                   scale = "response",
                   gg = T) +
  geom_point(data = analysis_data, 
             aes(x = duration_cont, y = Cilp_in_MC/Cilp_tot), 
             size = 1.5, shape = 1, position = position_jitter(width = 0.2)) + 
  # annotate('text', -Inf, Inf,
  #          label=paste("italic('p') == 0.0342"), parse=TRUE,
  #          hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% CILP-specific\nT cells in AC1", x = "Disease Duration\n(Years)") + 
  scale_y_continuous(labels = function(x) (x*100)) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25))

# Output Fig 2B
fig2B = grid.draw(rbind(ggplotGrob(fig2B_pt1),
              ggplotGrob(fig2B_pt2),
              size = "last"))


## VF in MC5 with rapid3
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 5, !is.na(rapid3_cont), Vimentin_Fibrinogen_tot >= 8)

antigenSpec_model = glm(Vimentin_Fibrinogen_in_MC/Vimentin_Fibrinogen_tot ~ 1 + rapid3_cont + age + sex,
                        weights = Vimentin_Fibrinogen_tot,
                        family = quasibinomial(),
                        data = analysis_data)

fig2C_pt1 = visreg(antigenSpec_model, 
                   xvar = "rapid3_cont", 
                   ylim = range(analysis_data$Vimentin_Fibrinogen_in_MC/analysis_data$Vimentin_Fibrinogen_tot), 
                   rug = F, 
                   line=list(col="red"), 
                   scale = "response",
                   gg = T) +
  geom_point(data = analysis_data, 
             aes(x = rapid3_cont, y = Vimentin_Fibrinogen_in_MC/Vimentin_Fibrinogen_tot), 
             size = 1.5, shape = 1, position = position_jitter(width = 0.1)) + 
  # annotate('text', -Inf, Inf,
  #          label=paste("italic('p') == 0.0442"), parse=TRUE,
  #          hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% Vim+Fib-specific\nT cells in AC5", x = "") + 
  scale_y_continuous(labels = function(x) (x*100))


## VF in MC6 with Rapid3_cat
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 6, !is.na(rapid3_cat), Vimentin_Fibrinogen_tot >= 8)

antigenSpec_model = glm(Vimentin_Fibrinogen_in_MC/Vimentin_Fibrinogen_tot ~ 1 + rapid3_cat + age + sex,
                        weights = Vimentin_Fibrinogen_tot,
                        family = quasibinomial(),
                        data = analysis_data) 

fig2C_pt2 = visreg(antigenSpec_model, 
                   xvar = "rapid3_cat", 
                   ylim = range(analysis_data$Vimentin_Fibrinogen_in_MC/analysis_data$Vimentin_Fibrinogen_tot), 
                   rug = F, 
                   line=list(col="red"), 
                   scale = "response",
                   gg = T) +
  geom_point(data = analysis_data, aes(x = ifelse(rapid3_cat == "0 - 1", 0.11,
                                                  ifelse(rapid3_cat == "1.3 - 2", 0.36,
                                                         ifelse(rapid3_cat == "2.3 - 4", 0.63, 0.89))), 
                                       y = Vimentin_Fibrinogen_in_MC/Vimentin_Fibrinogen_tot), 
             size = 1.5, shape = 1, position = position_jitter(width = 0.08)) +
  # annotate('text', -Inf, Inf,
  #          label=paste("italic('p') == 0.0437"), parse=TRUE,
  #          hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% Vim+Fib-specific\nT cells in AC6", x = "Disease Activity Score\n(Rapid3)") + 
  scale_y_continuous(labels = function(y) (y*100), 
                     breaks = c(0, .10, .20, .30, .40, .50, .60, .70),
                     limits = c(NA, .65)) 

# Output Fig 2C
fig2C = grid.draw(rbind(ggplotGrob(fig2C_pt1),
              ggplotGrob(fig2C_pt2),
              size = "last"))


# ## Cilp in MC6 with CRP
# analysis_data = mc_count_data %>%
#   dplyr::filter(metacluster == 6, !is.na(crp), Cilp_tot >= 8)
# 
# antigenSpec_model = glm(Cilp_in_MC/Cilp_tot ~ 1 + crp + age + sex,
#                         weights = Cilp_tot,
#                         family = quasibinomial(),
#                         data = analysis_data)
# 
# fig2D = visreg(antigenSpec_model, 
#                xvar = "crp", 
#                ylim = range(analysis_data$Cilp_in_MC/analysis_data$Cilp_tot), 
#                rug = F, 
#                line=list(col="red"),  #, lty = 2
#                scale = "response", 
#                gg = T) +
#   geom_point(data = analysis_data, 
#              aes(x = crp, y = Cilp_in_MC/Cilp_tot), 
#              size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
#   # annotate('text', -Inf, Inf,
#   #          label=paste("italic('p') == 0.0397"), parse=TRUE,
#   #          hjust=-0.3, vjust = 2, size = 3) +
#   labs(y = "% CILP-specific\nT cells in AC6", x = "CRP") + 
#   scale_y_continuous(labels = function(x) (x*100))


## Cilp in MC3 with CCP2
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 3, !is.na(ccp2_cat), Cilp_tot >= 8)

antigenSpec_model = glm(Cilp_in_MC/Cilp_tot ~ 1 + ccp2_cat + age + sex,
                        weights = Cilp_tot,
                        family = quasibinomial(),
                        data = analysis_data)

fig2D = visreg(antigenSpec_model,
               xvar = "ccp2_cat",
               ylim = range(analysis_data$Cilp_in_MC/analysis_data$Cilp_tot),
               rug = F,
               line=list(col="red"),  #, lty = 2
               scale = "response",
               gg = T) +
  geom_point(data = analysis_data,
             aes(x = ifelse(ccp2_cat == "Negative", 0.225, 0.775), 
                 y = Cilp_in_MC/Cilp_tot),
             size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
  # annotate('text', -Inf, Inf,
  #          label=paste("italic('p') == 0.0397"), parse=TRUE,
  #          hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% CILP-specific\nT cells in AC3", x = "CCP2") +
  scale_y_continuous(labels = function(x) (x*100))


## Cilp in MC8 with TNF_org
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 8, !is.na(tnfi_org_cat), Cilp_tot >= 8)

antigenSpec_model = glm(Cilp_in_MC/Cilp_tot ~ 1 + tnfi_org_cat + age + sex,
                        weights = Cilp_tot,
                        family = quasibinomial(),
                        data = analysis_data)

fig2E = visreg(antigenSpec_model,
               xvar = "tnfi_org_cat",
               ylim = range(analysis_data$Cilp_in_MC/analysis_data$Cilp_tot),
               rug = F,
               line=list(col="red"),  #, lty = 2
               scale = "response",
               gg = T) +
  geom_point(data = analysis_data,
             aes(x = ifelse(tnfi_org_cat == "dmard", 0.225, 0.775), 
                 y = Cilp_in_MC/Cilp_tot),
             size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
  # annotate('text', -Inf, Inf,
  #          label=paste("italic('p') == 0.0397"), parse=TRUE,
  #          hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% CILP-specific\nT cells in AC8", x = "Therapy") +
  scale_y_continuous(labels = function(x) (x*100))


## Agg in MC8 with CCP3
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 8, !is.na(ccp3_cat), Aggrecan_tot >= 8) %>%
  mutate(ccp3_cat = ifelse(ccp3_cat == "Negative", "Negative", "Positive") %>%
           factor(levels = c("Negative", "Positive")))

antigenSpec_model = glm(Aggrecan_in_MC/Aggrecan_tot ~ 1 + ccp3_cat + age + sex,
                        weights = Aggrecan_tot,
                        family = quasibinomial(),
                        data = analysis_data)

fig2F = visreg(antigenSpec_model,
               xvar = "ccp3_cat",
               ylim = range(analysis_data$Aggrecan_in_MC/analysis_data$Aggrecan_tot),
               rug = F,
               line=list(col="red"),  #, lty = 2
               scale = "response",
               gg = T) +
  geom_point(data = analysis_data,
             aes(x = ifelse(ccp3_cat == "Negative", 0.225, 0.775), 
             # aes(x = ifelse(ccp3_cat == "Negative", 0.11,
             #                ifelse(ccp3_cat == "Weak Positive", 0.36,
             #                       ifelse(ccp3_cat == "Moderate Positive", 0.63, 0.89))), 
                     y = Aggrecan_in_MC/Aggrecan_tot),
             size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
  # annotate('text', -Inf, Inf,
  #          label=paste("italic('p') == 0.0397"), parse=TRUE,
  #          hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% Aggrecan-specific\nT cells in AC8", x = "CCP3") +
  scale_y_continuous(labels = function(x) (x*100))


quartz(file = "NEW_Figure_2_no_pvals.pdf", type = "pdf", width = 12, height = 5, dpi = 300)
cowplot::plot_grid(fig2A_pt1, fig2B_pt1, fig2C_pt1, fig2F, fig2D, fig2E, fig2A_pt2, fig2B_pt2, fig2C_pt2, align = "hv", ncol = 5)
dev.off()

# # Wouldn't italicize the "p"
# cairo_pdf("Figure_2.pdf", width = 12, height = 5)
# cowplot::plot_grid(fig2A_pt1, fig2B_pt1, fig2C_pt1, fig2D, fig2A_pt2, fig2B_pt2, fig2C_pt2, align = "hv", ncol = 4)
# dev.off()
# 
# # Wouldn't print unicode characters
# cowplot::save_plot("Figure_2.pdf", test, ncol = 4, nrow = 2, base_height = 2.25, base_width = 3)




# 
# #### Age & sex-controlled models - LRT TEST ####
# ## CD4 in MC3 with Disease
# analysis_data = mc_count_data %>%
#   dplyr::filter(metacluster == 3, !is.na(disease_status), Influenza_tot >= 1)
# 
# total_cell_model = glm(Cells_in_MC/Cells_tot ~ 1 + disease_status + age + sex,
#                        weights = Cells_tot,
#                        family = quasibinomial(),
#                        data = analysis_data) 
# 
# fig1 = visreg(total_cell_model,
#               xvar = "disease_status", 
#               ylim = range(analysis_data$Cells_in_MC/analysis_data$Cells_tot), 
#               rug = F, 
#               alpha = 1,
#               line=list(col="red"),
#               scale = "response",
#               gg = T) +
#   geom_point(data = analysis_data, 
#              aes(x = ifelse(disease_status == "HC", 0.225, 0.775), 
#                  y = Cells_in_MC/Cells_tot, 
#                  shape = disease_status), 
#              size = 1.5, position = position_jitter(width = 0.1)) +
#   annotate('text', -Inf, Inf,
#            label=paste("italic('p') == 0.026"), parse=TRUE,
#            hjust=-0.3, vjust = 2, size = 3) +
#   labs(y = "% CD4+ cells in AC3", x = "") + 
#   scale_shape_manual(values = setNames(c(19, 1), c("HC", "RA")), guide = F) +
#   scale_y_continuous(labels = function(x) (x*100))
# 
# ## CD4 in MC2 with Rapid3cat
# analysis_data = mc_count_data %>%
#   dplyr::filter(metacluster == 2, !is.na(rapid3), Influenza_tot >= 1)
# 
# total_cell_model = glm(Cells_in_MC/Cells_tot ~ 1 + rapid3 + age + sex,
#                        weights = Cells_tot,
#                        family = quasibinomial(),
#                        data = analysis_data) 
# 
# fig2 = visreg(total_cell_model, 
#               xvar = "rapid3", 
#               ylim = range(analysis_data$Cells_in_MC/analysis_data$Cells_tot), 
#               rug = F, 
#               alpha = 1,
#               line=list(col="red"),
#               scale = "response",
#               gg = T) +
#   geom_point(data = analysis_data, 
#              aes(x = rapid3, 
#                  y = Cells_in_MC/Cells_tot), 
#              size = 1.5, shape = 1, position = position_jitter(width = 0.07)) +
#   annotate('text', -Inf, Inf,
#            label=paste("italic('p') == 0.011"), parse=TRUE,
#            hjust=-0.3, vjust = 2, size = 3) +
#   labs(y = "% CD4+ cells in AC2", x = "Disease Activity Score\n(Rapid3)") + 
#   scale_y_continuous(labels = function(x) (x*100))
# 
# 
# ## CD4 in MC6 with DAS28_CRP3_cat
# analysis_data = mc_count_data %>%
#   dplyr::filter(metacluster == 6, !is.na(das28_crp3_cat), Influenza_tot >= 1)
# 
# total_cell_model = glm(Cells_in_MC/Cells_tot ~ 1 + das28_crp3_cat + age + sex,
#                        weights = Cells_tot,
#                        family = quasibinomial(),
#                        data = analysis_data)
# 
# fig3 = visreg(total_cell_model,
#               xvar = "das28_crp3_cat",
#               ylim = range(analysis_data$Cells_in_MC/analysis_data$Cells_tot),
#               rug = F,
#               alpha = 1,
#               line=list(col="red"),
#               scale = "response",
#               gg = T) +
#   geom_point(data = analysis_data,
#              aes(x = ifelse(das28_crp3_cat == "< 2.6", 0.11,
#                         ifelse(das28_crp3_cat == "2.6 \u2264 3.2", 0.36,
#                                ifelse(das28_crp3_cat == "3.2  \u2264 5.1", 0.63, 0.89))),
#                  y = Cells_in_MC/Cells_tot),
#              size = 1.5, shape = 1, position = position_jitter(width = 0.07)) +
#   annotate('text', -Inf, Inf,
#            label=paste("italic('p') == 0.037"), parse=TRUE,
#            hjust=-0.3, vjust = 2, size = 3) +
#   labs(y = "% CD4+ cells in AC6", x = "Disease Activity Score\n(DAS28 - CRP3)") +
#   scale_y_continuous(labels = function(x) (x*100))
# 
# 
# ## Agg in MC1 with Duration
# analysis_data = mc_count_data %>%
#   dplyr::filter(metacluster == 1, !is.na(duration_cat), Aggrecan_tot >= 8)
# 
# antigenSpec_model = glm(Aggrecan_in_MC/Aggrecan_tot ~ 1 + duration_cat + age + sex,
#                         weights = Aggrecan_tot,
#                         family = quasibinomial(),
#                         data = analysis_data) 
# 
# fig4 = visreg(antigenSpec_model, 
#               xvar = "duration_cat", 
#               ylim = range(analysis_data$Aggrecan_in_MC/analysis_data$Aggrecan_tot), 
#               rug = F, 
#               alpha = 1,
#               line=list(col="red"),
#               scale = "response",
#               gg = T) +
#   geom_point(data = analysis_data, 
#              aes(x = ifelse(duration_cat == "\u2264 5", 0.225, 0.775), 
#                  y = Aggrecan_in_MC/Aggrecan_tot), 
#              size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
#   annotate('text', -Inf, Inf,
#            label=paste("italic('p') == 0.018"), parse=TRUE,
#            hjust=-0.3, vjust = 2, size = 3) +
#   labs(y = "% Aggrecan-specific\nT cells in AC1", x = "") + 
#   scale_y_continuous(labels = function(x) (x*100))
# 
# ## Agg in MC2 with Duration
# analysis_data = mc_count_data %>%
#   dplyr::filter(metacluster == 2, !is.na(duration_cat), Aggrecan_tot >= 8)
# 
# antigenSpec_model = glm(Aggrecan_in_MC/Aggrecan_tot ~ 1 + duration_cat + age + sex,
#                         weights = Aggrecan_tot,
#                         family = quasibinomial(),
#                         data = analysis_data) 
# 
# fig5 = visreg(antigenSpec_model, 
#               xvar = "duration_cat", 
#               ylim = range(analysis_data$Aggrecan_in_MC/analysis_data$Aggrecan_tot), 
#               rug = F, 
#               alpha = 1,
#               line=list(col="red"),
#               scale = "response",
#               gg = T) +
#   geom_point(data = analysis_data, 
#              aes(x = ifelse(duration_cat == "\u2264 5", 0.225, 0.775), 
#                  y = Aggrecan_in_MC/Aggrecan_tot), 
#              size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
#   annotate('text', -Inf, Inf,
#            label=paste("italic('p') == 0.022"), parse=TRUE,
#            hjust=-0.3, vjust = 2, size = 3) +
#   labs(y = "% Aggrecan-specific\nT cells in AC2", x = "") + 
#   scale_y_continuous(labels = function(x) (x*100))
# 
# 
# ## Agg in MC5 with Duration
# analysis_data = mc_count_data %>%
#   dplyr::filter(metacluster == 5, !is.na(duration_cat), Aggrecan_tot >= 8)
# 
# antigenSpec_model = glm(Aggrecan_in_MC/Aggrecan_tot ~ 1 + duration_cat + age + sex,
#                         weights = Aggrecan_tot,
#                         family = quasibinomial(),
#                         data = analysis_data) 
# 
# fig6 = visreg(antigenSpec_model, 
#               xvar = "duration_cat", 
#               ylim = range(analysis_data$Aggrecan_in_MC/analysis_data$Aggrecan_tot), 
#               alpha = 1,
#               rug = F, 
#               line=list(col="red"),
#               scale = "response",
#               gg = T) +
#   geom_point(data = analysis_data, 
#              aes(x = ifelse(duration_cat == "\u2264 5", 0.225, 0.775), 
#                  y = Aggrecan_in_MC/Aggrecan_tot), 
#              size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
#   annotate('text', -Inf, Inf,
#            label=paste("italic('p') == '4.0e-5'"), parse=TRUE,
#            hjust=-0.3, vjust = 2, size = 3) +
#   labs(y = "% Aggrecan-specific\nT cells in AC5", x = "") + 
#   scale_y_continuous(labels = function(x) (x*100))
# 
# 
# ## Agg in MC7 with DAS28_ESR3_cat
# analysis_data = mc_count_data %>%
#   dplyr::filter(metacluster == 7, !is.na(das28_cat), Aggrecan_tot >= 8)
# 
# antigenSpec_model = glm(Aggrecan_in_MC/Aggrecan_tot ~ 1 + das28_cat + age + sex,
#                        weights = Aggrecan_tot,
#                        family = quasibinomial(),
#                        data = analysis_data)
# 
# fig7 = visreg(antigenSpec_model,
#               xvar = "das28_cat",
#               ylim = range(analysis_data$Aggrecan_in_MC/analysis_data$Aggrecan_tot),
#               rug = F,
#               alpha = 1,
#               line=list(col="red"),
#               scale = "response",
#               gg = T) +
#   geom_point(data = analysis_data,
#              aes(x = ifelse(das28_cat == "< 2.6", 0.135,
#                             ifelse(das28_cat == "2.6 \u2264 3.2", 0.5, 0.855)),
#                  y = Aggrecan_in_MC/Aggrecan_tot),
#              size = 1.5, shape = 1, position = position_jitter(width = 0.07)) +
#   annotate('text', -Inf, Inf,
#            label=paste("italic('p') == 0.012"), parse=TRUE,
#            hjust=-0.3, vjust = 2, size = 3) +
#   labs(y = "% Aggrecan-specific\nT cells in AC7", x = "Disease Activity Score\n(DAS28 - ESR3)") +
#   scale_y_continuous(labels = function(x) (x*100))
# 
# 
# ## Agg in MC7 with TNFi
# analysis_data = mc_count_data %>%
#   dplyr::filter(metacluster == 7, !is.na(TNFi), Aggrecan_tot >= 8)
# 
# antigenSpec_model = glm(Aggrecan_in_MC/Aggrecan_tot ~ 1 + TNFi + age + sex,
#                         weights = Aggrecan_tot,
#                         family = quasibinomial(),
#                         data = analysis_data) 
# 
# fig8 = visreg(antigenSpec_model, 
#               xvar = "TNFi", 
#               ylim = range(analysis_data$Aggrecan_in_MC/analysis_data$Aggrecan_tot), 
#               rug = F, 
#               alpha = 1,
#               line=list(col="red"),
#               scale = "response",
#               gg = T) +
#   geom_point(data = analysis_data, 
#              aes(x = ifelse(TNFi == "DMARD", 0.225, 0.775), 
#                  y = Aggrecan_in_MC/Aggrecan_tot), 
#              size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
#   annotate('text', -Inf, Inf,
#            label=paste("italic('p') == 0.026"), parse=TRUE,
#            hjust=-0.3, vjust = 2, size = 3) +
#   labs(y = "% Aggrecan-specific\nT cells in AC7", x = "") + 
#   scale_y_continuous(labels = function(x) (x*100))
# 
# 
# ## Cilp in MC1 with Duration
# analysis_data = mc_count_data %>%
#   dplyr::filter(metacluster == 1, !is.na(duration), Cilp_tot >= 8)
# 
# antigenSpec_model = glm(Cilp_in_MC/Cilp_tot ~ 1 + duration + age + sex,
#                         weights = Cilp_tot,
#                         family = quasibinomial(),
#                         data = analysis_data)
# 
# fig9 = visreg(antigenSpec_model, 
#               xvar = "duration", 
#               ylim = range(analysis_data$Cilp_in_MC/analysis_data$Cilp_tot), 
#               rug = F, 
#               alpha = 1,
#               line=list(col="red"), 
#               scale = "response",
#               gg = T) +
#   geom_point(data = analysis_data, 
#              aes(x = duration, y = Cilp_in_MC/Cilp_tot), 
#              size = 1.5, shape = 1, position = position_jitter(width = 0.2)) + 
#   annotate('text', -Inf, Inf,
#            label=paste("italic('p') == 0.044"), parse=TRUE,
#            hjust=-0.3, vjust = 2, size = 3) +
#   labs(y = "% CILP-specific\nT cells in AC1", x = "Disease Duration\n(Years)") + 
#   scale_y_continuous(labels = function(x) (x*100)) +
#   scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25))
# 
# 
# ## Cilp in MC3 with Das28_ESR4 categorical
# analysis_data = mc_count_data %>%
#   dplyr::filter(metacluster == 3, !is.na(das28_esr4_cat), Cilp_tot >= 8)
# 
# antigenSpec_model = glm(Cilp_in_MC/Cilp_tot ~ 1 + das28_esr4_cat + age + sex,
#                         weights = Cilp_tot,
#                         family = quasibinomial(),
#                         data = analysis_data)
# 
# fig10 = visreg(antigenSpec_model,
#                xvar = "das28_esr4_cat",
#                ylim = range(analysis_data$Cilp_in_MC/analysis_data$Cilp_tot),
#                rug = F,
#                alpha = 1,
#                line=list(col="red"),  #, lty = 2
#                scale = "response",
#                gg = T) +
#   geom_point(data = analysis_data,
#              aes(x = ifelse(das28_esr4_cat == "< 2.6", 0.135,
#                             ifelse(das28_esr4_cat == "2.6 \u2264 3.2", 0.5, 0.855)),
#                  y = Cilp_in_MC/Cilp_tot),
#              size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
#   annotate('text', -Inf, Inf,
#            label=paste("italic('p') == 0.0024"), parse=TRUE,
#            hjust=-0.3, vjust = 2, size = 3) +
#   labs(y = "% CILP-specific\nT cells in AC3", x = "Disease Activity Score\n(DAS28 - ESR4)") +
#   scale_y_continuous(labels = function(x) (x*100))
# 
# 
# ## Cilp in MC5 with DAS28_ESR3
# analysis_data = mc_count_data %>%
#   dplyr::filter(metacluster == 5, !is.na(das28_cat), Cilp_tot >= 8)
# 
# antigenSpec_model = glm(Cilp_in_MC/Cilp_tot ~ 1 + das28_cat + age + sex,
#                         weights = Cilp_tot,
#                         family = quasibinomial(),
#                         data = analysis_data) 
# 
# fig11 = visreg(antigenSpec_model, 
#                xvar = "das28_cat", 
#                ylim = range(analysis_data$Cilp_in_MC/analysis_data$Cilp_tot), 
#                rug = F, 
#                alpha = 1,
#                line=list(col="red"),  #, lty = 2
#                scale = "response", 
#                gg = T) +
#   geom_point(data = analysis_data, 
#              aes(x = ifelse(das28_cat == "< 2.6", 0.135,
#                             ifelse(das28_cat == "2.6 \u2264 3.2", 0.5, 0.855)), 
#                             y = Cilp_in_MC/Cilp_tot), 
#              size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
#   annotate('text', -Inf, Inf,
#            label=paste("italic('p') == 0.047"), parse=TRUE,
#            hjust=-0.3, vjust = 2, size = 3) +
#   labs(y = "% CILP-specific\nT cells in AC5", x = "Disease Activity Score\n(DAS28 - ESR3)") + 
#   scale_y_continuous(labels = function(x) (x*100))
# 
# 
# ## Cilp in MC6 with CRP
# analysis_data = mc_count_data %>%
#   dplyr::filter(metacluster == 6, !is.na(crp), Cilp_tot >= 8)
# 
# antigenSpec_model = glm(Cilp_in_MC/Cilp_tot ~ 1 + crp + age + sex,
#                         weights = Cilp_tot,
#                         family = quasibinomial(),
#                         data = analysis_data)
# 
# fig12 = visreg(antigenSpec_model, 
#                xvar = "crp", 
#                ylim = range(analysis_data$Cilp_in_MC/analysis_data$Cilp_tot), 
#                rug = F, 
#                alpha = 1,
#                line=list(col="red"),  #, lty = 2
#                scale = "response", 
#                gg = T) +
#   geom_point(data = analysis_data, 
#              aes(x = crp, y = Cilp_in_MC/Cilp_tot), 
#              size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
#   annotate('text', -Inf, Inf,
#            label=paste("italic('p') == 0.0018"), parse=TRUE,
#            hjust=-0.3, vjust = 2, size = 3) +
#   labs(y = "% CILP-specific\nT cells in AC6", x = "CRP") + 
#   scale_y_continuous(labels = function(x) (x*100))
# 
# 
# ## VF in MC5 with rapid3
# analysis_data = mc_count_data %>%
#   dplyr::filter(metacluster == 5, !is.na(rapid3), Vimentin_Fibrinogen_tot >= 8)
# 
# antigenSpec_model = glm(Vimentin_Fibrinogen_in_MC/Vimentin_Fibrinogen_tot ~ 1 + rapid3 + age + sex,
#                         weights = Vimentin_Fibrinogen_tot,
#                         family = quasibinomial(),
#                         data = analysis_data)
# 
# fig13 = visreg(antigenSpec_model, 
#                xvar = "rapid3", 
#                ylim = range(analysis_data$Vimentin_Fibrinogen_in_MC/analysis_data$Vimentin_Fibrinogen_tot), 
#                rug = F, 
#                alpha = 1,
#                line=list(col="red"), 
#                scale = "response",
#                gg = T) +
#   geom_point(data = analysis_data, 
#              aes(x = rapid3, y = Vimentin_Fibrinogen_in_MC/Vimentin_Fibrinogen_tot), 
#              size = 1.5, shape = 1, position = position_jitter(width = 0.1)) + 
#   annotate('text', -Inf, Inf,
#            label=paste("italic('p') == 0.037"), parse=TRUE,
#            hjust=-0.3, vjust = 2, size = 3) +
#   labs(y = "% Vim+Fib-specific\nT cells in AC5", x = "Disease Activity Score\n(Rapid3)") + 
#   scale_y_continuous(labels = function(x) (x*100))
# 
# 
# ## VF in MC6 with Rapid3
# analysis_data = mc_count_data %>%
#   dplyr::filter(metacluster == 6, !is.na(rapid3), Vimentin_Fibrinogen_tot >= 8)
# 
# antigenSpec_model = glm(Vimentin_Fibrinogen_in_MC/Vimentin_Fibrinogen_tot ~ 1 + rapid3 + age + sex,
#                         weights = Vimentin_Fibrinogen_tot,
#                         family = quasibinomial(),
#                         data = analysis_data) 
# 
# fig14 = visreg(antigenSpec_model, 
#                xvar = "rapid3", 
#                ylim = range(analysis_data$Vimentin_Fibrinogen_in_MC/analysis_data$Vimentin_Fibrinogen_tot), 
#                rug = F, 
#                alpha = 1, 
#                line=list(col="red"), 
#                scale = "response",
#                gg = T) +
#   geom_point(data = analysis_data, 
#              aes(x = rapid3, y = Vimentin_Fibrinogen_in_MC/Vimentin_Fibrinogen_tot), 
#              size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
#   annotate('text', -Inf, Inf,
#            label=paste("italic('p') == 0.037"), parse=TRUE,
#            hjust=-0.3, vjust = 2, size = 3) +
#   labs(y = "% Vim+Fib-specific\nT cells in AC6", x = "Disease Activity Score\n(Rapid3)") + 
#   scale_y_continuous(labels = function(x) (x*100)) 
# 
# 
# ## VF in MC4 with DAS28_CRP3
# analysis_data = mc_count_data %>%
#   dplyr::filter(metacluster == 4, !is.na(das28_crp3_cat), Vimentin_Fibrinogen_tot >= 8)
# 
# antigenSpec_model = glm(Vimentin_Fibrinogen_in_MC/Vimentin_Fibrinogen_tot ~ 1 + das28_crp3_cat + age + sex,
#                         weights = Vimentin_Fibrinogen_tot,
#                         family = quasibinomial(),
#                         data = analysis_data) 
# 
# fig15 = visreg(antigenSpec_model, 
#                xvar = "das28_crp3_cat", 
#                ylim = range(analysis_data$Vimentin_Fibrinogen_in_MC/analysis_data$Vimentin_Fibrinogen_tot), 
#                rug = F, 
#                alpha = 1,
#                line=list(col="red"),  #, lty = 2
#                scale = "response", 
#                gg = T) +
#   geom_point(data = analysis_data, 
#              aes(x = ifelse(das28_crp3_cat == "< 2.6", 0.135,
#                             ifelse(das28_crp3_cat == "2.6 \u2264 3.2", 0.5, 0.855)), 
#                  y = Vimentin_Fibrinogen_in_MC/Vimentin_Fibrinogen_tot), 
#              size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
#   annotate('text', -Inf, Inf,
#            label=paste("italic('p') == '3.42e-6'"), parse=TRUE,
#            hjust=-0.3, vjust = 2, size = 3) +
#   labs(y = "% Vim+Fib-specific\nT cells in AC4", x = "Disease Activity Score\n(DAS28 - CRP3)") + 
#   scale_y_continuous(labels = function(x) (x*100))
# 
# 
# 
# 
# 
# quartz(file = "Figure_2_LRT.pdf", type = "pdf", width = 15, height = 7.5)
# cowplot::plot_grid(fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, 
#                    fig9, fig10, fig11, fig12, fig13, fig14, fig15,  
#                    align = "hv", ncol = 5)
# dev.off()
# 
# # # Wouldn't italicize the "p"
# # cairo_pdf("Figure_2.pdf", width = 12, height = 5)
# # cowplot::plot_grid(fig2A_pt1, fig2B_pt1, fig2C_pt1, fig2D, fig2A_pt2, fig2B_pt2, fig2C_pt2, align = "hv", ncol = 4)
# # dev.off()
# # 
# # # Wouldn't print unicode characters
# # cowplot::save_plot("Figure_2.pdf", test, ncol = 4, nrow = 2, base_height = 2.25, base_width = 3)
# 




#### RF and TNFi (left out of original) vs Tmr phenotypes ####
## CILP in MC8 with TNFi therapy
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 8, !is.na(TNFi), Cilp_tot >= 8)

antigenSpec_model = glm(Cilp_in_MC/Cilp_tot ~ 1 + TNFi + age + sex,
                        weights = Cilp_tot,
                        family = quasibinomial(),
                        data = analysis_data) 

fig2_pt0_A = visreg(antigenSpec_model, 
                   xvar = "TNFi", 
                   ylim = range(analysis_data$Cilp_in_MC/analysis_data$Cilp_tot), 
                   rug = F, 
                   line=list(col="red"),
                   scale = "response",
                   gg = T) +
  geom_point(data = analysis_data, 
             aes(x = ifelse(TNFi == "DMARD", 0.225, 0.775), 
                 y = Cilp_in_MC/Cilp_tot), 
             size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
  annotate('text', -Inf, Inf,
           label=paste("italic('p') == 0.0105"), parse=TRUE,
           hjust=-6, vjust = 2, size = 3) +
  labs(y = "% Cilp-specific\nT cells in AC8", x = "") + 
  scale_y_continuous(labels = function(x) (x*100), limits = c(-0.01, 1))


## Enolase in MC6 with VM_RF
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 6, !is.na(vm_rf_cont), Enolase_tot >= 8)

antigenSpec_model = glm(Enolase_in_MC/Enolase_tot ~ 1 + vm_rf_cont + age + sex,
                        weights = Enolase_tot,
                        family = quasibinomial(),
                        data = analysis_data)

fig2_pt0_B = visreg(antigenSpec_model, 
               xvar = "vm_rf_cont", 
               ylim = range(analysis_data$Enolase_in_MC/analysis_data$Enolase_tot), 
               rug = F, 
               line=list(col="red"),  #, lty = 2
               scale = "response", 
               gg = T) +
  geom_point(data = analysis_data, 
             aes(x = vm_rf_cont, y = Enolase_in_MC/Enolase_tot), 
             size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
  annotate('text', -Inf, Inf,
           label=paste("italic('p') == 0.0197"), parse=TRUE,
           hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% Enolase-specific\nT cells in AC6", x = "RF") + 
  scale_y_continuous(labels = function(x) (x*100), limits = c(-0.01, 1)) 


## VimFib in MC6 with VM_RF
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 4, !is.na(vm_rf_cont), Vimentin_Fibrinogen_tot >= 8)

antigenSpec_model = glm(Vimentin_Fibrinogen_in_MC/Vimentin_Fibrinogen_tot ~ 1 + vm_rf_cont + age + sex,
                        weights = Vimentin_Fibrinogen_tot,
                        family = quasibinomial(),
                        data = analysis_data)

fig2_pt0_C = visreg(antigenSpec_model, 
                    xvar = "vm_rf_cont", 
                    ylim = range(analysis_data$Vimentin_Fibrinogen_in_MC/analysis_data$Vimentin_Fibrinogen_tot), 
                    rug = F, 
                    line=list(col="red"),  #, lty = 2
                    scale = "response", 
                    gg = T) +
  geom_point(data = analysis_data, 
             aes(x = vm_rf_cont, y = Vimentin_Fibrinogen_in_MC/Vimentin_Fibrinogen_tot), 
             size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
  annotate('text', -Inf, Inf,
           label=paste("italic('p') == 0.0191"), parse=TRUE,
           hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% Vim-Fib-specific\nT cells in AC4", x = "RF") + 
  scale_y_continuous(labels = function(x) (x*100), limits = c(-0.01, 1)) 

quartz(file = file.path("Output", "191016_300dpi_plots", "dotplots", "Unshown_fdrSig_GLMs.pdf"), 
       type = "pdf", width = 9, height = 2.5, dpi = 300)
cowplot::plot_grid(fig2_pt0_A, fig2_pt0_B, fig2_pt0_C, align = "hv", ncol = 3)
dev.off()






#### Load data for specificity plots ####
# sp_data_path = file.path("Output", "190320", "heatmaps", 
#                          "191105specificity_counts_for_modeling.csv")
# 
# # Clean up data
# spec_count_data = read_csv(sp_data_path)
# 
# das28_relabel_vector = c(`Near remission` = "< 2.6",
#                          `Low` = "2.6 \u2264 3.2",
#                          `Moderate` = "3.2  \u2264 5.1",
#                          `High` = "> 5.1")
# 
# spec_count_data = spec_count_data %>%
#   mutate(disease_status = factor(disease_status, levels = c("HC", "RA")),
#          duration_cat = case_when(duration_cat == "Recent onset" ~ "\u2264 5",
#                                   duration_cat == "Long standing" ~ ">5") %>%
#            factor(levels = c("\u2264 5", ">5")),
#          rapid3_cat =  case_when(rapid3_cat == "Near remission" ~ "0 - 1",
#                                  rapid3_cat == "Low" ~ "1.3 - 2",
#                                  rapid3_cat == "Moderate" ~ "2.3 - 4",
#                                  rapid3_cat == "High" ~ "4.3 - 10",
#                                  is.na(rapid3_cat) ~ NA_character_) %>%
#            factor(levels = c("0 - 1", "1.3 - 2", "2.3 - 4", "4.3 - 10"), ordered = T),
#          das28_cat = factor(das28_cat, levels = c("Near remission", "Low", "Moderate", "High"), ordered = T) %>%
#            recode_factor(!!!das28_relabel_vector),
#          das28_esr4_cat = factor(das28_esr4_cat, levels = c("Near remission", "Low", "Moderate", "High"), ordered = T) %>%
#            recode_factor(!!!das28_relabel_vector),
#          das28_crp3_cat = factor(das28_crp3_cat, levels = c("Near remission", "Low", "Moderate", "High"), ordered = T) %>%
#            recode_factor(!!!das28_relabel_vector),
#          das28_crp4_cat = factor(das28_crp4_cat, levels = c("Near remission", "Low", "Moderate", "High"), ordered = T) %>%
#            recode_factor(!!!das28_relabel_vector),
#          crp_cat = factor(crp_cat, levels = c("Normal", "High")),
#          TNFi = factor(TNFi, levels = c("DMARD", "TNFi")),
#          biologic_cat = case_when(biologic_cat == T ~ "Biologic",
#                                   biologic_cat == F ~ "No biologic",
#                                   T ~ NA_character_) %>%
#            factor(levels = c("No biologic", "Biologic")),
#          age_cat = factor(age_cat, levels = c("Young", "Old")),
#          sex = factor(sex),
#          group = factor(group, levels = c("HC", "RO", "LS"))) 

spec_count_data = mc_count_data %>% 
  dplyr::select(-(metacluster:Cells_in_MC)) %>% 
  distinct() %>% 
  rowwise() %>%
  mutate(RA_tot = sum(Aggrecan_tot, Cilp_tot, Enolase_tot, Vimentin_Fibrinogen_tot))


#### Age & sex-controlled specificity models - WALD TEST ####
## Agg with Disease
analysis_data = spec_count_data %>%
  dplyr::filter(!is.na(disease_status), RA_tot >= 1)

specificity_model = glm(Aggrecan_tot/RA_tot ~ 1 + disease_status + age + sex,
                       weights = RA_tot,
                       family = quasibinomial(),
                       data = analysis_data) 

fig3A_pt1 = visreg(specificity_model, 
                   xvar = "disease_status", 
                   ylim = range(analysis_data$Aggrecan_tot/analysis_data$RA_tot), 
                   rug = F, 
                   line=list(col="red"),
                   scale = "response",
                   gg = T) +
  geom_point(data = analysis_data, 
             aes(x = ifelse(disease_status == "HC", 0.225, 0.775), 
                 y = Aggrecan_tot/RA_tot, 
                 shape = disease_status), 
             size = 1.5, position = position_jitter(width = 0.1)) +
  annotate('text', -Inf, Inf,
           label=paste("italic('p') == 0.001"), parse=TRUE,
           hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% RA Tmr+ cells\nspecific to Aggrecan", x = "") + 
  scale_shape_manual(values = setNames(c(19, 1), c("HC", "RA")), guide = F) +
  scale_y_continuous(labels = function(x) (x*100))

## CILP with Disease
analysis_data = spec_count_data %>%
  dplyr::filter(!is.na(disease_status), RA_tot >= 1)

specificity_model = glm(Cilp_tot/RA_tot ~ 1 + disease_status + age + sex,
                        weights = RA_tot,
                        family = quasibinomial(),
                        data = analysis_data) 

fig3A_pt2 = visreg(specificity_model, 
                   xvar = "disease_status", 
                   ylim = range(analysis_data$Cilp_tot/analysis_data$RA_tot), 
                   rug = F, 
                   line=list(col="red"),
                   scale = "response",
                   gg = T) +
  geom_point(data = analysis_data, 
             aes(x = ifelse(disease_status == "HC", 0.225, 0.775), 
                 y = Cilp_tot/RA_tot, 
                 shape = disease_status), 
             size = 1.5, position = position_jitter(width = 0.1)) +
  annotate('text', -Inf, Inf,
           label=paste("italic('p') == 0.027"), parse=TRUE,
           hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% RA Tmr+ cells\nspecific to CILP", x = "") + 
  scale_shape_manual(values = setNames(c(19, 1), c("HC", "RA")), guide = F) +
  scale_y_continuous(labels = function(x) (x*100))


## Eno with Disease
analysis_data = spec_count_data %>%
  dplyr::filter(!is.na(disease_status), RA_tot >= 1)

specificity_model = glm(Enolase_tot/RA_tot ~ 1 + disease_status + age + sex,
                        weights = RA_tot,
                        family = quasibinomial(),
                        data = analysis_data) 

fig3A_pt3 = visreg(specificity_model, 
                   xvar = "disease_status", 
                   ylim = range(analysis_data$Enolase_tot/analysis_data$RA_tot), 
                   rug = F, 
                   line=list(col="red"),
                   scale = "response",
                   gg = T) +
  geom_point(data = analysis_data, 
             aes(x = ifelse(disease_status == "HC", 0.225, 0.775), 
                 y = Enolase_tot/RA_tot, 
                 shape = disease_status), 
             size = 1.5, position = position_jitter(width = 0.1)) +
  annotate('text', -Inf, Inf,
           label=paste("italic('p') == 0.027"), parse=TRUE,
           hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% RA Tmr+ cells\nspecific to Enolase", x = "") + 
  scale_shape_manual(values = setNames(c(19, 1), c("HC", "RA")), guide = F) +
  scale_y_continuous(labels = function(x) (x*100))


## Agg with Age
analysis_data = spec_count_data %>%
  dplyr::filter(!is.na(disease_status), RA_tot >= 1)

specificity_model = glm(Aggrecan_tot/RA_tot ~ 1 + age,
                        weights = RA_tot,
                        family = quasibinomial(),
                        data = analysis_data) 

fig3_pt0 = visreg(specificity_model, 
                   xvar = "age", 
                   ylim = range(analysis_data$Aggrecan_tot/analysis_data$RA_tot), 
                   rug = F, 
                   line=list(col="red"),
                   scale = "response",
                   gg = T) +
  geom_point(data = analysis_data, 
             aes(x = age, 
                 y = Aggrecan_tot/RA_tot), 
             size = 1.5, shape = 1,
             position = position_jitter(width = 0.1)) +
  annotate('text', -Inf, Inf,
           label=paste("italic('p') == 0.042"), parse=TRUE,
           hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% RA Tmr+ cells\nspecific to Aggrecan", x = "Age") + 
  scale_y_continuous(labels = function(x) (x*100))


## Eno with Biologic
analysis_data = full_updated_phenos %>% 
  dplyr::select(subject = subj, biologic_cat = bio_cat, other_id) %>% 
  right_join(spec_count_data) %>%
  mutate(biologic_cat = ifelse(biologic_cat == F, "No biologic", "Biologic") %>%
           factor(levels = c("No biologic", "Biologic"))) %>%
  dplyr::filter(!is.na(biologic_cat), RA_tot >= 1)

specificity_model = glm(Enolase_tot/RA_tot ~ 1 + biologic_cat + age + sex,
                        weights = RA_tot,
                        family = quasibinomial(),
                        data = analysis_data) 

fig3B_pt1 = visreg(specificity_model, 
                   xvar = "biologic_cat", 
                   ylim = range(analysis_data$Enolase_tot/analysis_data$RA_tot), 
                   rug = F, 
                   line=list(col="red"),
                   scale = "response",
                   gg = T) +
  geom_point(data = analysis_data, 
             aes(x = ifelse(biologic_cat == "No biologic", 0.225, 0.775), 
                 y = Enolase_tot/RA_tot), 
             shape = 1, size = 1.5, 
             position = position_jitter(width = 0.1)) +
  annotate('text', -Inf, Inf,
           label=paste("italic('p') == 0.027"), parse=TRUE,
           hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% RA Tmr+ cells\nspecific to Enolase", x = "Treatment") + 
  scale_y_continuous(labels = function(x) (x*100))

## Explore for sanity check
# test_data = analysis_data %>%
#   dplyr::select(subject, biologic_cat, Aggrecan_tot, Cilp_tot,
#                 Enolase_tot, Vimentin_Fibrinogen_tot, RA_tot) %>%
#   mutate(pct_eno = Enolase_tot/RA_tot*100)
# subset = test_data %>%
#   dplyr::filter(subject %in% c("1b_028", "1a_029", "1b_054",
#                                "1a_028", "1a_033", "1a_049"))
# 
# subset %>%
#   dplyr::filter(subject %in% c("1a_028", "1a_049", "1a_029", "1b_054")) %>%
#   pivot_longer(c("Aggrecan_tot", "Enolase_tot", "Cilp_tot", "Vimentin_Fibrinogen_tot"), 
#                names_to = "specificity", values_to = "tmr_count") %>%
#   mutate(subject = factor(subject, levels = c("1a_029", "1b_054", "1a_028", "1a_049"))) %>%
#   ggplot(aes(x = subject, # biologic_cat
#              y = tmr_count, 
#              color = specificity)) +
#     geom_point(size = 2) +
#     # facet_wrap(~specificity) +
#     scale_color_manual("Specificity", values = c("#0D0887FF", "#9443A8FF", "#CC4678FF", "#F89441FF"))  + 
#     labs(y = "# Tmr+ cells", x = "Subject")


## Eno with Rapid3
analysis_data = spec_count_data %>%
  dplyr::filter(!is.na(rapid3), RA_tot >= 1)

specificity_model = glm(Enolase_tot/RA_tot ~ 1 + rapid3 + age + sex,
                        weights = RA_tot,
                        family = quasibinomial(),
                        data = analysis_data) 

fig3B_pt2 = visreg(specificity_model, 
                   xvar = "rapid3", 
                   ylim = range(analysis_data$Enolase_tot/analysis_data$RA_tot), 
                   rug = F, 
                   line=list(col="red"),
                   scale = "response",
                   gg = T) +
  geom_point(data = analysis_data, 
             aes(x = rapid3, 
                 y = Enolase_tot/RA_tot), 
             shape = 1, size = 1.5, 
             position = position_jitter(width = 0.1)) +
  annotate('text', -Inf, Inf,
           label=paste("italic('p') == 0.022"), parse=TRUE,
           hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% RA Tmr+ cells\nspecific to Enolase", x = "DISEASE  ACTIVITY SCORE\n(RAPID3)") + 
  scale_y_continuous(labels = function(x) (x*100))



quartz(file = file.path("Output", "191016_300dpi_plots", "dotplots", "Figure_3_with_pvals.pdf"),
       type = "pdf", width = 9, height = 5, dpi = 300)
cowplot::plot_grid(fig3A_pt1, fig3A_pt2, fig3A_pt3, fig3_pt0, fig3B_pt1, fig3B_pt2, align = "hv", ncol = 3)
dev.off()



#### Flu test 6/26/2020 ####

## Agg in MC5 with Duration
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 5, !is.na(duration_cat), Aggrecan_tot >= 8)

antigenSpec_model = glm(Aggrecan_in_MC/Aggrecan_tot ~ 1 + duration_cat + age + sex,
                        weights = Aggrecan_tot,
                        family = quasibinomial(),
                        data = analysis_data) 

visreg(antigenSpec_model, 
                   xvar = "duration_cat", 
                   ylim = range(analysis_data$Aggrecan_in_MC/analysis_data$Aggrecan_tot), 
                   rug = F, 
                   line=list(col="red"),
                   scale = "response",
                   gg = T) +
  geom_point(data = analysis_data, 
             aes(x = ifelse(duration_cat == "\u2264 5", 0.225, 0.775), 
                 y = Aggrecan_in_MC/Aggrecan_tot), 
             size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
  # annotate('text', -Inf, Inf,
  #          label=paste("italic('p') == 0.0026"), parse=TRUE,
  #          hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% Aggrecan-specific\nT cells in AC5", x = "") + 
  scale_y_continuous(labels = function(x) (x*100))

## Flu in MC5 with Duration
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 5, !is.na(duration_cat), Aggrecan_tot >= 8)

antigenSpec_model = glm(Influenza_in_MC/Influenza_tot ~ 1 + duration_cat + age + sex,
                        weights = Influenza_tot,
                        family = quasibinomial(),
                        data = analysis_data) 

visreg(antigenSpec_model, 
       xvar = "duration_cat", 
       ylim = range(analysis_data$Influenza_in_MC/analysis_data$Influenza_tot), 
       rug = F, 
       line=list(col="red"),
       scale = "response",
       gg = T) +
  geom_point(data = analysis_data, 
             aes(x = ifelse(duration_cat == "\u2264 5", 0.225, 0.775), 
                 y = Influenza_in_MC/Influenza_tot), 
             size = 1.5, shape = 1, position = position_jitter(width = 0.1)) +
  labs(y = "% Flu-specific\nT cells in AC5", x = "") + 
  scale_y_continuous(labels = function(x) (x*100))

## Flu in MC5 with Continuous Duration
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 5, !is.na(duration_cont), Aggrecan_tot >= 8)

antigenSpec_model = glm(Influenza_in_MC/Influenza_tot ~ 1 + duration_cont + age + sex,
                        weights = Influenza_tot,
                        family = quasibinomial(),
                        data = analysis_data) 

visreg(antigenSpec_model, 
                   xvar = "duration_cont", 
                   ylim = range(analysis_data$Influenza_in_MC/analysis_data$Influenza_tot), 
                   rug = F, 
                   line=list(col="red"), 
                   scale = "response",
                   gg = T) +
  geom_point(data = analysis_data, 
             aes(x = duration_cont, y = Influenza_in_MC/Influenza_tot), 
             size = 1.5, shape = 1, position = position_jitter(width = 0.2)) + 
  labs(y = "% Flu-specific\nT cells in AC5", x = "Disease Duration\n(Years)") + 
  scale_y_continuous(labels = function(x) (x*100)) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25))


## Cilp in MC1 with Duration
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 1, !is.na(duration_cont), Cilp_tot >= 8)

antigenSpec_model = glm(Cilp_in_MC/Cilp_tot ~ 1 + duration_cont + age + sex,
                        weights = Cilp_tot,
                        family = quasibinomial(),
                        data = analysis_data)

age_model = glm(Cilp_in_MC/Cilp_tot ~ 1 + age,
                weights = Cilp_tot,
                family = quasibinomial(),
                data = analysis_data)
  
visreg(antigenSpec_model, 
                   xvar = "duration_cont", 
                   ylim = range(analysis_data$Cilp_in_MC/analysis_data$Cilp_tot), 
                   rug = F, 
                   line=list(col="red"), 
                   scale = "response",
                   gg = T) +
  geom_point(data = analysis_data, 
             aes(x = duration_cont, y = Cilp_in_MC/Cilp_tot), 
             size = 1.5, shape = 1, position = position_jitter(width = 0.2)) + 
  labs(y = "% CILP-specific\nT cells in AC1", x = "Disease Duration\n(Years)") + 
  scale_y_continuous(labels = function(x) (x*100)) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25))

visreg(age_model, 
       xvar = "age", 
       ylim = range(analysis_data$Cilp_in_MC/analysis_data$Cilp_tot), 
       rug = F, 
       line=list(col="red"), 
       scale = "response",
       gg = T) +
  geom_point(data = analysis_data, 
             aes(x = age, y = Cilp_in_MC/Cilp_tot), 
             size = 1.5, shape = 1, position = position_jitter(width = 0.2)) + 
  labs(y = "% Flu-specific\nT cells in AC1", x = "Age\n(Years)") + 
  scale_y_continuous(labels = function(x) (x*100)) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25))


## Flu in MC1 with Duration
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 1, !is.na(duration_cont), Cilp_tot >= 8)

antigenSpec_model = glm(Influenza_in_MC/Influenza_tot ~ 1 + duration_cont + age + sex,
                        weights = Influenza_tot,
                        family = quasibinomial(),
                        data = analysis_data)

age_model = glm(Influenza_in_MC/Influenza_tot ~ 1 + age,
                weights = Influenza_tot,
                family = quasibinomial(),
                data = analysis_data)

visreg(antigenSpec_model, 
       xvar = "duration_cont", 
       ylim = range(analysis_data$Influenza_in_MC/analysis_data$Influenza_tot), 
       rug = F, 
       line=list(col="red"), 
       scale = "response",
       gg = T) +
  geom_point(data = analysis_data, 
             aes(x = duration_cont, y = Influenza_in_MC/Influenza_tot), 
             size = 1.5, shape = 1, position = position_jitter(width = 0.2)) + 
  labs(y = "% Influenza-specific\nT cells in AC1", x = "Disease Duration\n(Years)") + 
  scale_y_continuous(labels = function(x) (x*100)) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25))


## CD4s in MC1 with Duration
analysis_data = mc_count_data %>%
  dplyr::filter(metacluster == 1, !is.na(duration_cont), Cilp_tot >= 8)

antigenSpec_model = glm(NonSpec_in_MC/NonSpec_tot ~ 1 + duration_cont + age + sex,
                        weights = NonSpec_tot,
                        family = quasibinomial(),
                        data = analysis_data)

visreg(antigenSpec_model, 
       xvar = "duration_cont", 
       ylim = range(analysis_data$NonSpec_in_MC/analysis_data$NonSpec_tot), 
       rug = F, 
       line=list(col="red"), 
       scale = "response",
       gg = T) +
  geom_point(data = analysis_data, 
             aes(x = duration_cont, y = NonSpec_in_MC/NonSpec_tot), 
             size = 1.5, shape = 1, position = position_jitter(width = 0.2)) + 
  # annotate('text', -Inf, Inf,
  #          label=paste("italic('p') == 0.0342"), parse=TRUE,
  #          hjust=-0.3, vjust = 2, size = 3) +
  labs(y = "% Non-specific\nT cells in AC1", x = "Disease Duration\n(Years)") + 
  scale_y_continuous(labels = function(x) (x*100)) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25))

