library(tidyverse)
library(RColorBrewer)
library(readxl)

setwd("~/Box/RA Cliff Metaclustering")

dir.create("Output")
dir.create(file.path("Output", format(Sys.Date(), "%y%m%d")))
dir.create(file.path("Output", format(Sys.Date(), "%y%m%d"), "mosaic_plots"))
fname_prefix_mosaic <- file.path("Output", format(Sys.Date(), "%y%m%d"), "mosaic_plots", 
                                 format(Sys.Date(), "%y%m%d"))


# Reset ggplot defaults to match GLM figs
theme_set(theme_classic(10) + 
            theme(axis.title = element_text(size=10, face="bold"),
                  axis.text.x = element_text(size=10, face="bold", color = "black"),
                  axis.text.y = element_text(size=8, color = "black")))

validation_data = read_xlsx("DoD VM 200129.xlsx", sheet = 2)

count_data = validation_data %>% 
  dplyr::select(`Other id`, `Draw Sequence`, `Aggrecan cts`, 
                `CILP cts`, `Enolase cts`, `VF cts`) %>% 
  pivot_longer(ends_with(" cts"), names_to = "specificity", values_to = "Count") %>% 
  mutate(specificity = str_remove_all(specificity, " cts")) 

per_million_cd4_data = validation_data %>% 
  dplyr::select(`Other id`, `Draw Sequence`, `Aggrecan per mil CD4`, 
                `CILP per mil CD4`, `Enolase per mil CD4`, `VF per mil CD4`) %>% 
  pivot_longer(ends_with(" per mil CD4"), names_to = "specificity", values_to = "Count per million CD4") %>% 
  mutate(specificity = str_remove_all(specificity, " per mil CD4")) 

total_ra_per_million_cd4_data = validation_data %>% 
  dplyr::select(`Other id`, `Draw Sequence`, `Combined RA per mil CD4`) %>% 
  pivot_longer(ends_with(" per mil CD4"), names_to = "specificity", values_to = "Count per million CD4") %>% 
  mutate(specificity = str_remove_all(specificity, " per mil CD4")) 

## Option 1: Freq per million stacked, paired bar plots
ggplot(per_million_cd4_data, 
       aes(x = factor(`Draw Sequence`), 
           y = `Count per million CD4`, 
           fill = specificity)) +
  geom_bar(stat = "identity") +
  facet_grid(cols = vars(`Other id`), scales = "free_x", space = "free_x") +
  labs(x = "Draw #") +
  scale_fill_manual("Specificity", values = c("#0D0887FF", "#9443A8FF", "#CC4678FF", "#F89441FF")) 
ggsave("Fig3d_for_further_editing.pdf", width = 6, height = 4)

ggplot(count_data, 
       aes(x = `Draw Sequence`, 
           y = `Count`, 
           fill = specificity)) +
  geom_bar(stat = "identity") +
  facet_wrap(~`Other id`, nrow=1) +
  scale_fill_manual("Specificity", values = c("#0D0887FF", "#9443A8FF", "#CC4678FF", "#F89441FF")) 

## Option 2: Freq per total Tmr+ stacked, paired bar plots
ggplot(per_million_cd4_data, 
       aes(x = `Draw Sequence`, 
           y = `Count per million CD4`, 
           fill = specificity)) +
  geom_bar(position="fill", stat="identity") +
  facet_wrap(~`Other id`, nrow=1) +
  scale_fill_manual("Specificity", values = c("#0D0887FF", "#9443A8FF", "#CC4678FF", "#F89441FF")) 

ggplot(count_data, 
       aes(x = `Draw Sequence`, 
           y = `Count`, 
           fill = specificity)) +
  geom_bar(position="fill", stat="identity") +
  facet_wrap(~`Other id`, nrow=1) +
  scale_fill_manual("Specificity", values = c("#0D0887FF", "#9443A8FF", "#CC4678FF", "#F89441FF")) +
  labs(y = "% of Total RA Tmr+")

## Option 3:
ggplot(total_ra_per_million_cd4_data, 
       aes(x = `Draw Sequence`, 
           y = `Count per million CD4`, 
           fill = specificity)) +
  geom_bar(stat = "identity") +
  facet_wrap(~`Other id`, nrow=1) +
  scale_y_continuous(limits = c(0,17)) +
  scale_fill_manual("Specificity", values = c("darkgrey")) 

ggplot(per_million_cd4_data, 
       aes(x = `Draw Sequence`, 
           y = `Count per million CD4`, 
           color = specificity)) +
  geom_point(alpha = 0.8) +
  geom_line() +
  scale_x_continuous(breaks = c(1,2,3)) +
  scale_y_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 17)) +
  facet_wrap(~`Other id`, nrow=1) +
  scale_color_manual("Specificity", values = c("#0D0887FF", "#9443A8FF", "#CC4678FF", "#F89441FF")) 


