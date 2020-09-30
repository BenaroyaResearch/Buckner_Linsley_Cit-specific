#### Set up #### 
setwd("~/Box/RA Cliff Metaclustering")
library(tidyverse)
library(openxlsx)
library(geepack)
library(broom)
library(ggcorrplot)
library(readxl)

options(warn = 1)

# Set up the ggplot default params
theme_set(theme_bw(12) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size=15, face="bold", margin = margin(10,0,10,0)),
                  axis.text.x = element_text(angle=45, hjust = 1)))

#### Load MC/pheno data #### 
mc_data_path = file.path("Output", "191211", "heatmaps", 
                         "191211metacluster_counts_for_modeling.csv")
mc_count_data = read_csv(mc_data_path)

load("20191210_combined_phenos.RData")
phenos <- abbrev_phenos

mc_count_data = mc_count_data %>%
  left_join(phenos, by = c("subject" = "subj")) %>%
  left_join(full_updated_phenos %>% 
              dplyr::select(subject = subj, biologic_cat = bio_cat, other_id)) %>% 
  rowwise() %>%
  mutate(biologic_cat = ifelse(biologic_cat == F, "No biologic", "Biologic") %>%
           factor(levels = c("No biologic", "Biologic")),
         ccp2_cont = case_when(disease_status == "HC" ~ NA_real_,
                               disease_status == "RA" ~ ccp2_cont),
         rf_cont = case_when(disease_status == "HC" ~ NA_real_,
                             disease_status == "RA" ~ rf_cont),
         RAspec_in_MC = Aggrecan_in_MC + Cilp_in_MC + 
           Enolase_in_MC + Vimentin_Fibrinogen_in_MC, 
         RAspec_tot = Aggrecan_tot + Cilp_tot + Enolase_tot + Vimentin_Fibrinogen_tot)
  
specificities = c("Aggrecan", "Cilp", "Enolase", "Vimentin_Fibrinogen", "Influenza",
                  "NonSpec") #, "RAspec")


#### Start with quasibinomial glm models testing cell phenos related to clin vars ####
logistic_prop_test = function(data_source, group_var, tmr_cutoff = 8, 
                              ag_specs, fdr_method = "holm", use_correction = T){
  mcs = unique(data_source$metacluster)
  results = list()
  options(warn=2) # This helps catch glm convergence errors below
  
  for(spec in ag_specs){
    tmr_in_mc = paste0(spec, "_in_MC")
    tmr_tot = paste0(spec, "_tot")
    
    for(mc in mcs){
      
      analysis_data = data_source %>%
        dplyr::filter(metacluster == mc, !is.na(get(group_var)), get(tmr_tot) >= tmr_cutoff)
      
      if(use_correction == T) {
        
        antigenSpec_model = try(glm(get(tmr_in_mc)/get(tmr_tot) ~ 1 + get(group_var) + age + sex,
                                    weights = get(tmr_tot),
                                    family = quasibinomial(),
                                    data = analysis_data), FALSE)
        
      } else {
        
        antigenSpec_model = try(glm(get(tmr_in_mc)/get(tmr_tot) ~ 1 + get(group_var), 
                                    weights = get(tmr_tot),
                                    family = quasibinomial(),
                                    data = analysis_data), FALSE)

      }
      
      failed <- inherits(antigenSpec_model,"try-error")
      
      if(failed || !antigenSpec_model$converged) {
        
        results[[spec]] = results[[spec]] %>%
          bind_rows(tibble(comparison = group_var,
                           metacluster = mc, 
                           odds_ratio = NA,
                           p_val = NA))
        
      } else {
        
        antigenSpec_res = tidy(antigenSpec_model)
        
        results[[spec]] = results[[spec]] %>%
          bind_rows(tibble(comparison = group_var,
                           metacluster = mc, 
                           odds_ratio = exp(coef(antigenSpec_model))[[2]],
                           p_val = antigenSpec_res$p.value[2]))
      }
      
    }
    
    results[[spec]] = results[[spec]] %>%
      mutate(fdr = p.adjust(p_val, fdr_method))
    
  }
  
  results = results %>%
    bind_rows(.id = "tmr") %>%
    arrange(fdr) %>%
    return()
}


#### Test the base models
sex_res = logistic_prop_test(mc_count_data, "sex", 8, specificities, use_correction = F) 
age_cont_res = logistic_prop_test(mc_count_data, "age", 8, specificities, use_correction = F)

# Plot the base models
uncorrected_models = bind_rows(sex_res, age_cont_res) %>%
  mutate(p_cat = case_when(p_val < 0.001 ~ "< 0.001",
                           p_val < 0.01 ~ "< 0.01",
                           p_val < 0.05 ~ "< 0.05",
                           p_val < 0.1 ~ "< 0.1",
                           T ~ "ns") %>%
           factor(levels = c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "ns")))

coolwarm_hcl <- colorspace::diverging_hcl(5,
                                          h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7))
# 10^seq(-.55, .55, length.out = 5)

ggplot(uncorrected_models, 
       aes(x = factor(metacluster, levels = c(1:8)), y = comparison, color = odds_ratio)) + 
  geom_point(aes(size = p_cat), show.legend = TRUE) +
  facet_wrap(~ factor(tmr, levels = c("NonSpec", "Aggrecan", "Cilp", "Enolase", "Vimentin_Fibrinogen")), ncol = 1) +
  scale_color_gradient2("OR", low = "navy", mid = coolwarm_hcl[3], 
                        high = "darkred", midpoint = 1, 
                        breaks = c(0.29, 0.54, 1, 1.87, 3.47),
                        limits = c(0.28, 3.5), oob = scales::squish) +
  scale_size_manual(values = setNames(c(6, 5, 4, 2.5, 1.25, 0), c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "ns", ""))) +
  labs(y = "", x = "Aligned cluster", color = "OR", size = "Uncorrected\np-value")
ggsave(file.path("Output", "191016_300dpi_plots", "corr_plots", 
                 paste0(format(Sys.Date(), "%y%m%d"), "_no_covariates_corr_plot.png")), 
       width = 4, height = 5, dpi = 300)

# Remove age/sex-uncorrected results to keep them separate from corrected models
rm(list = ls(pattern = "_res"))

#### Test the age & sex-adjusted models
disease_res = logistic_prop_test(mc_count_data, "disease_status", 8, specificities)  

duration_cat_res = logistic_prop_test(mc_count_data, "duration_cat", 8, specificities) 
duration_cont_res = logistic_prop_test(mc_count_data, "duration_cont", 8, specificities)  

r3_cat_res = logistic_prop_test(mc_count_data, "rapid3_cat", 8, specificities)  
r3_cont_res = logistic_prop_test(mc_count_data, "rapid3_cont", 8, specificities)  

d28_cat_res = logistic_prop_test(mc_count_data, "das28_cat", 8, specificities) 
d28_cont_res = logistic_prop_test(mc_count_data, "das28_cont", 8, specificities)  

biologic_cat_res = logistic_prop_test(mc_count_data, "biologic_cat", 8, specificities)

tnfi_org_cat_res = logistic_prop_test(mc_count_data, "tnfi_org_cat", 8, specificities) 

crp_cont_res = logistic_prop_test(mc_count_data, "crp_cont", 8, specificities)  
rf_cont_res = logistic_prop_test(mc_count_data, "rf_cont", 8, specificities)  

ccp2_cont_res = logistic_prop_test(mc_count_data, "ccp2_cont", 8, specificities)  


## Combine results
res_list = ls(pattern = "_res$")

corrected_models = sapply(res_list, function(res_df) get(res_df), simplify = F) %>%
  bind_rows() %>%
  mutate(comparison = str_replace_all(comparison, "_cat", "_categorical") %>%
           str_replace_all( "_cont", "_continuous") %>%
           factor(levels = c("tnfi_org_categorical", "biologic_categorical", "rf_continuous", "crp_continuous", 
                             "ccp2_continuous", "das28_categorical", "das28_continuous", "rapid3_categorical", 
                             "rapid3_continuous", "duration_categorical", "duration_continuous", "disease_status")),
         input_sep = case_when(comparison == "disease_status" ~ "top",
                               T ~ "bottom") %>%
           factor(levels = c("top", "bottom")),
         p_cat = case_when(p_val < 0.001 ~ "< 0.001",
                           p_val < 0.01 ~ "< 0.01",
                           p_val < 0.05 ~ "< 0.05",
                           p_val < 0.1 ~ "< 0.1",
                           is.na(p_val) ~ "",
                           T ~ "ns") %>%
           factor(levels = c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "ns", ""))) %>%
  arrange(tmr, metacluster)


ggplot(corrected_models %>% dplyr::filter(tmr == "NonSpec"), 
       aes(x = factor(metacluster, levels = c(1:8)), y = comparison, color = odds_ratio)) + 
  geom_point(aes(size = p_cat), show.legend = TRUE) +
  scale_color_gradient2("OR", low = "navy", mid = coolwarm_hcl[3], 
                        high = "darkred", midpoint = 1, 
                        breaks = c(0.29, 0.54, 1, 1.87, 3.47),
                        limits = c(0.28, 3.5), oob = scales::squish) +
  scale_size_manual(values = setNames(c(6, 5, 4, 2.5, 1.25, 0), c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "ns", ""))) +
  facet_grid(rows = vars(input_sep), scales = "free_y", space = "free_y") +
  labs(y = "", x = "Aligned cluster", title = "Nonspecific CD4+ Cells",
       color = "OR", size = "Uncorrected\np-value") +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank())
ggsave(file.path("Output", "191016_300dpi_plots", "corr_plots", 
                 paste0(format(Sys.Date(), "%y%m%d"), "_age_sex_cd4_corr_plot.svg")), 
       width = 6, height = 5.75, dpi = 300)

ggplot(corrected_models %>% dplyr::filter(tmr == "Aggrecan"), 
       aes(x = factor(metacluster, levels = c(1:8)), y = comparison, color = odds_ratio)) + 
  geom_point(aes(size = p_cat), show.legend = TRUE) +
  scale_color_gradient2("OR", low = "navy", mid = coolwarm_hcl[3], 
                        high = "darkred", midpoint = 1, 
                        breaks = c(0.29, 0.54, 1, 1.87, 3.47),
                        limits = c(0.28, 3.5), oob = scales::squish) +
  scale_size_manual(values = setNames(c(6, 5, 4, 2.5, 1.25, 0), c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "ns", ""))) +
  facet_grid(rows = vars(input_sep), scales = "free_y", space = "free_y") +
  labs(y = "", x = "Aligned cluster", title = "Aggrecan-specific CD4+ Cells",
       color = "OR", size = "Uncorrected\np-value") +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank())
ggsave(file.path("Output", "191016_300dpi_plots", "corr_plots", 
                 paste0(format(Sys.Date(), "%y%m%d"), "_age_sex_Agg_corr_plot.svg")), 
       width = 6, height = 5.75, dpi = 300)

ggplot(corrected_models %>% dplyr::filter(tmr == "Cilp"), 
       aes(x = factor(metacluster, levels = c(1:8)), y = comparison, color = odds_ratio)) + 
  geom_point(aes(size = p_cat), show.legend = TRUE) +
  scale_color_gradient2("OR", low = "navy", mid = coolwarm_hcl[3], 
                        high = "darkred", midpoint = 1, 
                        breaks = c(0.29, 0.54, 1, 1.87, 3.47),
                        limits = c(0.28, 3.5), oob = scales::squish) +
  scale_size_manual(values = setNames(c(6, 5, 4, 2.5, 1.25, 0), c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "ns", ""))) +
  facet_grid(rows = vars(input_sep), scales = "free_y", space = "free_y") +
  labs(y = "", x = "Aligned cluster", title = "CILP-specific CD4+ Cells",
       color = "OR", size = "Uncorrected\np-value") +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank())
ggsave(file.path("Output", "191016_300dpi_plots", "corr_plots", 
                 paste0(format(Sys.Date(), "%y%m%d"), "_age_sex_Cilp_corr_plot.svg")), 
       width = 6, height = 5.75, dpi = 300)

ggplot(corrected_models %>% dplyr::filter(tmr == "Enolase"), 
       aes(x = factor(metacluster, levels = c(1:8)), y = comparison, color = odds_ratio)) + 
  geom_point(aes(size = p_cat), show.legend = TRUE) +
  scale_color_gradient2("OR", low = "navy", mid = coolwarm_hcl[3], 
                        high = "darkred", midpoint = 1, 
                        breaks = c(0.29, 0.54, 1, 1.87, 3.47),
                        limits = c(0.28, 3.5), oob = scales::squish) +
  scale_size_manual(values = setNames(c(6, 5, 4, 2.5, 1.25, 0), c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "ns", ""))) +
  facet_grid(rows = vars(input_sep), scales = "free_y", space = "free_y") +
  labs(y = "", x = "Aligned cluster", title = "Enolase-specific CD4+ Cells",
       color = "OR", size = "Uncorrected\np-value") +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank())
ggsave(file.path("Output", "191016_300dpi_plots", "corr_plots", 
                 paste0(format(Sys.Date(), "%y%m%d"), "_age_sex_Eno_corr_plot.svg")), 
       width = 6, height = 5.75, dpi = 300)

ggplot(corrected_models %>% dplyr::filter(tmr == "Vimentin_Fibrinogen"), 
       aes(x = factor(metacluster, levels = c(1:8)), y = comparison, color = odds_ratio)) + 
  geom_point(aes(size = p_cat), show.legend = TRUE) +
  scale_color_gradient2("OR", low = "navy", mid = coolwarm_hcl[3], 
                        high = "darkred", midpoint = 1, 
                        breaks = c(0.29, 0.54, 1, 1.87, 3.47),
                        limits = c(0.28, 3.5), oob = scales::squish) +
  scale_size_manual(values = setNames(c(6, 5, 4, 2.5, 1.25, 0), c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "ns", ""))) +
  facet_grid(rows = vars(input_sep), scales = "free_y", space = "free_y") +
  labs(y = "", x = "Aligned cluster", title = "Vim/Fib-specific CD4+ Cells",
       color = "OR", size = "Uncorrected\np-value") +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank())
ggsave(file.path("Output", "191016_300dpi_plots", "corr_plots", 
                 paste0(format(Sys.Date(), "%y%m%d"), "_age_sex_VF_corr_plot.svg")), 
       width = 6, height = 5.75, dpi = 300)



#### Report specificity/clinical variable model results
res_df = bind_rows(uncorrected_models, 
                   corrected_models %>% droplevels()) %>%
  mutate(comparison = str_replace_all(comparison, "_cat$", "_categorical"),
         comparison = str_replace_all(comparison, "_cont$", "_continuous")) %>%
  arrange(tmr, comparison)

test_out_path = file.path("Output", "191211", "dotplots",
                          "quasibinomial_glm_phenotype_results.xlsx")

write.xlsx(list("sig_results" = res_df %>% dplyr::filter(fdr < 0.05), 
                "all_results" = res_df), test_out_path)



#### Restructure specificity data ####
spec_count_data = mc_count_data %>% 
  dplyr::select(-(metacluster:Cells_in_MC)) %>% 
  distinct() %>% 
  rowwise() %>%
  mutate(RA_tot = sum(Aggrecan_tot, Cilp_tot, Enolase_tot, Vimentin_Fibrinogen_tot))

ra_specificities = c("Aggrecan", "Cilp", "Enolase", "Vimentin_Fibrinogen")

# Remove previous model results to keep them separate from this analysis
rm(list = ls(pattern = "_res"))
              

logistic_spec_test = function(data_source, group_var, tmr_cutoff = 4, 
                              ag_specs, fdr_method = "holm", use_correction = T){
  results = list()
  options(warn=2) # This helps catch glm convergence errors below
  
  for(spec in ag_specs){
    tmr_tot = paste0(spec, "_tot")
    
    analysis_data = data_source %>%
      dplyr::filter(!is.na(get(group_var)), RA_tot >= tmr_cutoff)
    
    if(use_correction == T) {
      
      antigenSpec_model = try(glm(get(tmr_tot)/RA_tot ~ 1 + get(group_var) + age + sex,
                                  weights = RA_tot,
                                  family = quasibinomial(),
                                  data = analysis_data), FALSE)
      
    } else {
      
      antigenSpec_model = try(glm(get(tmr_tot)/RA_tot ~ 1 + get(group_var), 
                                  weights = RA_tot,
                                  family = quasibinomial(),
                                  data = analysis_data), FALSE)
      
    }
    
    failed <- inherits(antigenSpec_model,"try-error")
    
    if(failed || !antigenSpec_model$converged) {
      
      results[[spec]] = results[[spec]] %>%
        bind_rows(tibble(comparison = group_var,
                         odds_ratio = NA,
                         p_val = NA))
      
    } else {
      
      antigenSpec_res = tidy(antigenSpec_model)
      
      results[[spec]] = results[[spec]] %>%
        bind_rows(tibble(comparison = group_var,
                         odds_ratio = exp(coef(antigenSpec_model))[[2]],
                         p_val = antigenSpec_res$p.value[2]))
    }
    
  }

  results = results %>%
    bind_rows(.id = "tmr") %>%
    mutate(fdr = p.adjust(p_val, fdr_method)) %>%
    arrange(fdr) %>%
    return()
}


#### Test the base models
sex_res = logistic_spec_test(spec_count_data, "sex", 1, ra_specificities, use_correction = F) 
age_cont_res = logistic_spec_test(spec_count_data, "age", 1, ra_specificities, use_correction = F)

# Plot the base models
uncorrected_models = bind_rows(sex_res, age_cont_res) %>%
  mutate(p_cat = case_when(p_val < 0.001 ~ "< 0.001",
                           p_val < 0.01 ~ "< 0.01",
                           p_val < 0.05 ~ "< 0.05",
                           p_val < 0.1 ~ "< 0.1",
                           is.na(p_val) ~ "",
                           T ~ "ns") %>%
           factor(levels = c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "ns", "")))

ggplot(uncorrected_models, 
       aes(x = tmr, y = comparison, color = odds_ratio)) + 
  geom_point(aes(size = p_cat), show.legend = TRUE) +
  scale_color_gradient2("OR", low = "navy", mid = "#DDDDDD", 
                        high = "darkred", midpoint = 1, 
                        breaks = c(0.29, 0.54, 1, 1.87, 3.47),
                        limits = c(0.28, 3.5), oob = scales::squish) +
  scale_size_manual(values = setNames(c(6, 5, 4, 2.5, 1.25, 0), c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "ns", ""))) +
  labs(y = "", x = "Specificity", color = "OR", size = "Uncorrected\np-value")
ggsave(file.path("Output", "191016_300dpi_plots", "corr_plots", 
                 paste0(format(Sys.Date(), "%y%m%d"), "_no_covariates_specificity_corr_plot.png")), 
       width = 6, height = 5, dpi = 300)

# Remove age/sex-uncorrected results to keep them separate from corrected models
rm(list = ls(pattern = "_res"))

#### Test the age & sex-adjusted models
disease_res = logistic_spec_test(spec_count_data, "disease_status", 1, ra_specificities)  

duration_cat_res = logistic_spec_test(spec_count_data, "duration_cat", 1, ra_specificities) 
duration_cont_res = logistic_spec_test(spec_count_data, "duration_cont", 1, ra_specificities)  

r3_cat_res = logistic_spec_test(spec_count_data, "rapid3_cat", 1, ra_specificities)  
r3_cont_res = logistic_spec_test(spec_count_data, "rapid3_cont", 1, ra_specificities)  

d28_cat_res = logistic_spec_test(spec_count_data, "das28_cat", 1, ra_specificities) 
d28_cont_res = logistic_spec_test(spec_count_data, "das28_cont", 1, ra_specificities)  

biologic_cat_res = logistic_spec_test(spec_count_data, "biologic_cat", 1, ra_specificities)

tnfi_org_cat_res = logistic_spec_test(spec_count_data, "tnfi_org_cat", 1, ra_specificities) 

crp_cont_res = logistic_spec_test(spec_count_data, "crp_cont", 1, ra_specificities)  
rf_cont_res = logistic_spec_test(spec_count_data, "rf_cont", 1, ra_specificities)  

ccp2_cont_res = logistic_spec_test(spec_count_data, "ccp2_cont", 1, ra_specificities)  


## Combine results
res_list = ls(pattern = "_res$")

corrected_models = sapply(res_list, function(res_df) get(res_df), simplify = F) %>%
  bind_rows() %>%
  mutate(comparison = str_replace_all(comparison, "_cat$", "_categorical") %>%
           str_replace_all( "_cont$", "_continuous") %>%
           factor(levels = c("tnfi_org_categorical", "biologic_categorical", "rf_continuous", "crp_continuous", 
                             "ccp2_continuous", "das28_categorical", "das28_continuous", "rapid3_categorical", 
                             "rapid3_continuous", "duration_categorical", "duration_continuous", "disease_status")),
         input_sep = case_when(comparison == "disease_status" ~ "top",
                               T ~ "bottom") %>%
           factor(levels = c("top", "bottom")),
         p_cat = case_when(p_val < 0.001 ~ "< 0.001",
                           p_val < 0.01 ~ "< 0.01",
                           p_val < 0.05 ~ "< 0.05",
                           p_val < 0.1 ~ "< 0.1",
                           is.na(p_val) ~ "",
                           T ~ "ns") %>%
           factor(levels = c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "ns", ""))) 
 


ggplot(corrected_models, 
       aes(x = tmr, y = comparison, color = odds_ratio)) + 
  geom_point(aes(size = p_cat), show.legend = TRUE) +
  scale_color_gradient2("OR", low = "navy", mid = "#DDDDDD", 
                        high = "darkred", midpoint = 1, 
                        breaks = c(0.29, 0.54, 1, 1.87, 3.47),
                        limits = c(0.28, 3.5), oob = scales::squish) +
  scale_size_manual(values = setNames(c(6, 5, 4, 2.5, 1.25, 0), c("< 0.001", "< 0.01", "< 0.05", "< 0.1", "ns", ""))) +
  facet_grid(rows = vars(input_sep), scales = "free_y", space = "free_y") +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank()) +
  labs(y = "", x = "", title = "Proportion of RA Tmr+ cells\nwithin a given specificity",
       color = "OR", size = "Uncorrected\np-value")
ggsave(file.path("Output", "191016_300dpi_plots", "corr_plots", 
                 paste0(format(Sys.Date(), "%y%m%d"), "_age_sex_specificity_corr_plot.svg")), 
       width = 6, height = 5.75, dpi = 300)


#### Report specificity/clinical variable model results
res_df = bind_rows(uncorrected_models, corrected_models) %>%
  mutate(comparison = str_replace_all(comparison, "_cat", "_categorical")) %>%
  arrange(tmr, comparison)

test_out_path = file.path("Output", "190320", "dotplots",
                          "quasibinomial_glm_specificity_results.xlsx")

write.xlsx(list("sig_results" = res_df %>% dplyr::filter(fdr < 0.05), 
                "all_results" = res_df), test_out_path)


#### Report specificity dataframe for outside graphing use
sp_data_path = file.path("Output", "190320", "heatmaps", 
                         "191105specificity_counts_for_modeling.csv")

write_csv(spec_count_data, sp_data_path)



###############################################
#### Permutations for quasibinomial models ####
###############################################

hold_steady_annos = mc_count_data %>% 
  dplyr::select(subject:Cells_tot, RAspec_in_MC, RAspec_tot)
randomize_anno = mc_count_data %>% 
  dplyr::select(subject, other_id:biologic_cat) %>% 
  distinct()

permute_glms = function(var_tested, results_df, n_perm){

    original_results = results_df %>%
      mutate(perm_count = 0)

    for(i in 1:n_perm){
      permuted_data = randomize_anno %>%
        transform(subject = sample(subject)) %>%
        right_join(hold_steady_annos, by = "subject")
      
      permuted_output = logistic_prop_test(permuted_data, var_tested, 8, specificities)
      
      original_results = original_results %>%
        left_join(permuted_output, by = c("tmr", "comparison", "metacluster")) %>%
        mutate(perm_count = case_when(odds_ratio.x > 1 & 
                                        odds_ratio.y > odds_ratio.x ~ perm_count + 1, 
                                      odds_ratio.x < 1 & 
                                        odds_ratio.y < odds_ratio.x ~ perm_count + 1,
                                      T ~ perm_count + 0)) %>%
        dplyr::select(tmr, comparison, metacluster, odds_ratio = odds_ratio.x, 
                      p_val = p_val.x, fdr = fdr.x, perm_count)
    }
    
    original_results = original_results %>%
      mutate(perm_freq = perm_count/n_perm)
}

disease_perm = permute_glms("disease_status", disease_res, 10000)
duration_cat_perm = permute_glms("duration_cat", duration_cat_res, 10000) 
duration_cont_perm = permute_glms("duration_cont", duration_cont_res, 10000)  

r3_cat_perm = permute_glms("rapid3_cat", r3_cat_res, 10000)  
r3_cont_perm = permute_glms("rapid3_cont", r3_cont_res, 10000)  

crp_cont_perm = permute_glms("crp_cont", crp_cont_res, 10000)  






