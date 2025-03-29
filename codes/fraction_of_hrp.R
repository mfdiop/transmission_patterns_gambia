

### **Step 1: Load Required Packages**

library(tidyverse)
library(ggpubr)
library(Hmisc)  # For Clopper-Pearson Exact Interval and Wilson CI

source("codes/functions.R")

fig.format <- c("pdf", "tiff")


### **Figure A: Distribution of Pairwise IBD (F Statistic)**
# Load Pairwise IBD data
ibd_data <- read_tsv("results/tables/IBD_GambiaOldandRecentSamples.tsv") %>% 
   rename(IBD = dcifer)

# Sample metadata
metadata <-read_tsv("data/metadata/metadata_to_use.tsv") %>%
   select(-c(2, 3, 7:9, 16:24)) %>%
   rename(indv = samples) %>% 
   mutate(Admin.level.1 = case_when(sites %in% c("Bessi", "Brikama", "Serrekunda") ~ 'West Coast', 
                                    sites %in% c("Chogen", "Yallal Ba") ~ 'North Bank',
                                    sites %in% c("Sinchu", "Dongoro Ba") ~ 'Lower River',
                                    sites %in% c("Farafenni", "Sare Seedy" ) ~ 'Central River',
                                    sites %in% c("Basse", "Njayel", "Madina Samako", 
                                                 "Sare Wuro", "Gunjur Kunta") ~ 'Upper River'),
          Year2 = if_else(Year2 == 1991, 1990, Year2),
          Admin.level.1 = factor(Admin.level.1, levels = c("West Coast", "North Bank", "Lower River",
                                                           "Central River", "Upper River")),
          sites = factor(sites, levels = c("Serrekunda", "Brikama", "Bessi", "Chogen", "Yallal Ba",
                                           "Sinchu", "Dongoro Ba", "Farafenni", "Sare Seedy",
                                           "Njayel", "Madina Samako", "Basse", "Sare Wuro", "Gunjur Kunta")))

counts(metadata, column = "Year2")

# Example: Metadata contains a column 'SampleID' with 720 sample names
samples <- metadata$indv

# Filter pairwise IBD dataset to include only pairs where both Sample1 and Sample2 are in metadata
filtered_ibd <- ibd_data %>%
   filter(p1 %in% samples & p2 %in% samples)

# Merge metadata with IBD data
ibd_metadata <- filtered_ibd %>%
   left_join(metadata, by = c("p1" = "indv")) %>%
   left_join(metadata, by = c("p2" = "indv"), suffix = c(".1", ".2"))

# location <- metadata %>% pull(sites) %>% unique() %>% sort()
rm(ibd_data)

# =========================================
# Considering only pairs with IBD >= 0.25
# =========================================

# # Perform ANOVA for differences in highly-related fractions across years
# anova_time <- aov(IBD ~ as.factor(Year2.1), data = data)
# summary(anova_time)


# ==========================================================
# Considering all pairs and compare less vs highly related
# ==========================================================
# Compute standard error (SE)
data_time <- calculate_relatedness_stats(ibd_metadata, threshold = 0.25, 
   filter_columns = c("Year2.1", "Year2.2"), group_column = "Year2.1") %>%
   rename(Year = Year2.1)  # Rename for clarity

# =============================================================
# Distribution of fraction of pairs above the chosen threshold
# =============================================================
frac_hrp_v1 <- ggplot(data_time, aes(x = as.factor(Year), y = HighlyRelatedFraction)) +
   geom_bar(stat = "identity", fill = "#0073C2FF", alpha = 0.8, color = "black", linewidth = 1) + # "#386CB0" "steelblue" "#73C2FF"
   geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.1, color = "black") +
   scale_y_continuous(expand = c(0, 0), limits = c(0, 0.035),
                      breaks = seq(0, 0.035, 0.005), #c(0, 0.01, 0.02, 0.03, 0.035),
                      labels = c("0", "0.005", "0.01", "0.015", "0.02", "0.025", "0.03", "0.035")) +
   ggnetwork::theme_blank() +
   labs(x = "", y = "Fraction of highly-related pairs") +
   theme(axis.text = element_text(colour = "black", size = 10, face = "bold"),
         axis.title = element_text(colour = "black", size = 16, face = "bold"),
         axis.line.y = element_line(linewidth = 1, lineend = "square"),
         axis.ticks = element_line(color = 'black', linewidth = 1),
         axis.ticks.length = unit(.25, "cm"))

print(frac_hrp_v1)

# ===============================================
# Mean IBD Distribution of related pairs by year
# ===============================================

frac_hrp_v0 <- ggplot(data_time, aes(x = as.factor(Year), y = MeanIBD)) +
   geom_bar(stat = "identity", fill = "#0073C2FF", alpha = 0.8, color = "black", linewidth = 1) +
   geom_errorbar(aes(ymin = L95CI, ymax = U95CI), width = 0.1, color = "black") +
   scale_y_continuous(expand = c(0, 0), limits = c(0, .04),
                      labels = c("0", "0.01", "0.02", "0.03", "0.04")) +
   # stat_compare_means(aes(group = as.factor(Year)), label = "p.signif") +
   ggnetwork::theme_blank() +
   labs(x = "", y = "Fraction of highly-related pairs") +
   theme(axis.text = element_text(colour = "black", size = 10, face = "bold"),
         axis.title = element_text(colour = "black", size = 16, face = "bold"),
         axis.line.y = element_line(linewidth = 1, lineend = "square"),
         axis.ticks = element_line(color = 'black', linewidth = 1),
         axis.ticks.length = unit(.25, "cm"))

print(frac_hrp_v0)

# =================================================================
# ---- (B) Partition by Study Site (City) with Distance Line ----

# Calculate the fraction of highly-related pairs per year
wth_sites <- calculate_relatedness_stats(ibd_metadata, threshold = 0.25, 
                                         filter_columns = c("sites.1", "sites.2"), 
                                         group_column = "sites.1") %>% 
   rename(sites = sites.1) %>% 
   add_column(type = "Within site")

btw_sites <- calculate_relatedness_stats(ibd_metadata, threshold = 0.25, 
                                         filter_columns = c("sites.1", "sites.2"), 
                                         group_column = "sites.1",
                                         filter_operator = "!=") %>% 
   mutate(sites = paste( sites.1, sites.2, sep = '-')) %>% 
   select(-sites.1, -sites.2) %>% 
   relocate(sites, .before = MeanIBD) %>% 
   add_column(type = "Between sites")

frac_sites <- rbind.data.frame(wth_sites, btw_sites) %>% 
   mutate(type = factor(type, levels = c("Within site", "Between sites"))) %>% 
   filter(HighlyRelatedFraction != 0 & !str_detect(sites, paste(c("Chogen", "Yallal Ba", "Sinchu", "Dongoro Ba", "Sare Seedy", "Njayel"), collapse = "|"))) %>%
   # filter(HighlyRelatedFraction != 0) %>% 
   ggplot(aes(x = sites, y = HighlyRelatedFraction)) +
   geom_bar(stat = "identity", aes(fill = type), alpha = 0.8, color = "black") +
   geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.1, color = "black") +
   scale_y_continuous(expand = c(0, 0), limits = c(0, 0.1),
                      breaks = seq(0, 0.1, 0.025),
                      labels = as.character(seq(0, 0.1, 0.025))) +
   scale_fill_manual(values = c("#0073C2FF", "#8F7700FF"), name = "") + # #EFC000FF
   ggnetwork::theme_blank() +
   labs(x = "", y = "Fraction of highly-related pairs") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 10, face = "bold"),
         axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
         axis.title = element_text(colour = "black", size = 16, face = "bold"),
         axis.line.y = element_line(linewidth = 1, lineend = "square"),
         axis.ticks.y = element_line(color = 'black', linewidth = 1),
         axis.ticks.length = unit(.25, "cm"),
         legend.text = element_text(size = 10, color = 'black', face = 'bold'),
         legend.position = "inside",
         legend.position.inside = c(.75, .75))

combined <- frac_hrp_v1 + frac_sites
print(combined)

ggsave("results/figures/fraction_of_hrp.tiff", width = 310, 
       height = 190, units = "mm", dpi = 600, plot = combined)


# # ===============================================
# # ---- (A) Partition by Time ----
# # Compute standard error (SE) and Calculate 
# # the fraction of highly-related pairs per year
# 
# data_time <- ibd_metadata %>%
#    filter(IBD >= .25) %>%
#    filter(Year2.1 == Year2.2) %>%
#    group_by(Year2.1) %>%
#    summarise(
#       fraction = mean(IBD),
#       sampleSize = n(),
#       SE = sqrt(fraction * (1 - fraction) / sampleSize)
#    ) %>% 
#    rename(Year = Year2.1)
# 
# # ========================================================
# # Extract highly related pairs within location
# wth_sites <- ibd_metadata %>%
#    filter(IBD >= .25) %>%
#    filter(sites.1 == sites.2) %>%
#    group_by(sites.1) %>%
#    summarise(
#       fraction = mean(IBD),
#       sampleSize = n(),
#       SE = sqrt(fraction * (1 - fraction) / sampleSize),
#       .groups = "drop") %>% 
#    rename(sites = sites.1) %>% 
#    add_column(type = "Within City")
# 
# # Extract highly related pairs between locations
# btw_sites <- ibd_metadata %>%
#    filter(IBD >= .25) %>%
#    filter(sites.1 != sites.2) %>%
#    group_by(sites.1, sites.2) %>%
#    summarise(
#       fraction = mean(IBD),
#       sampleSize = n(),
#       SE = sqrt(fraction * (1 - fraction) / sampleSize),
#       .groups = "drop") %>% 
#    mutate(sites = paste( sites.1, sites.2)) %>% 
#    select(-sites.1, -sites.2) %>% 
#    relocate(sites, .before = fraction) %>% 
#    add_column(type = "Between Cities")
# 
# # Plot the highly-related fraction by study site
# rbind.data.frame(wth_sites, btw_sites) %>% 
#    filter(sampleSize > 2) %>% 
#    ggplot(aes(x = sites, y = fraction)) +
#    geom_bar(stat = "identity", aes(fill = type), alpha = 0.8, color = 'black', linewidth = .5) + # fill = "darkorange"
#    geom_errorbar(aes(ymin = fraction - SE, ymax = fraction + SE), width = 0.2, color = "black") +
#    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 2, 0.25),
#                       labels = as.character(seq(0, 2, 0.25))) +
#    theme_blank() +
#    labs(x = "City (Ordered by Distance)",
#         y = "Fraction of Highly-Related Pairs",
#         caption = "Red dashed line represents inter-study site great-circle distance.") +
#    theme(axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
#          axis.text.y = element_text(colour = "black", size = 13, face = "bold"),
#          axis.title = element_text(colour = "black", size = 16, face = "bold"),
#          axis.line.y = element_line(linewidth = 1, lineend = "square"),
#          axis.ticks.y = element_line(color = 'black', linewidth = 1),
#          axis.ticks.length = unit(.25, "cm"))

# ===============================================================================
