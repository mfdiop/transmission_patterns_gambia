
library(tidyverse)
library(patchwork)

# Load summary data for DRC
rmcl <- read_tsv("results/tables/RMCOIL_gmb_samples.tsv") %>% select(1:3)
fws <- read_tsv("results/tables/Fws.tsv")
grc <- readxl::read_xlsx("results/tables/coi_from_GRC.xls") %>% type_convert()

# Load metadata file
metadata <- read_tsv("data/metadata/metadata_to_use.tsv") %>%
   select(-c(2, 3, 7:9, 16:24)) %>%
   rename(indv = samples) %>% 
   mutate(Admin.level.1 = case_when(sites %in% c("Bessi", "Brikama", "Serrekunda") ~ 'West Coast', 
                                    sites %in% c("Chogen", "Yallal Ba") ~ 'North Bank',
                                    sites %in% c("Sinchu", "Dongoro Ba") ~ 'Lower River',
                                    sites %in% c("Farafenni", "Sare Seedy" ) ~ 'Central River',
                                    sites %in% c("Basse", "Njayel", "Madina Samako", 
                                                 "Sare Wuro", "Gunjur Kunta") ~ 'Upper River'),
           Admin.level.1 = factor(Admin.level.1, levels = c("West Coast", "North Bank", "Lower River",
                                                           "Central River", "Upper River")),
          sites = factor(sites, levels = c("Serrekunda", "Brikama", "Bessi", "Chogen", "Yallal Ba",
                                           "Sinchu", "Dongoro Ba", "Farafenni", "Sare Seedy",
                                           "Njayel", "Madina Samako", "Basse", "Sare Wuro", "Gunjur Kunta")))

# Combined COI data and metadata
data <- rmcl %>% 
   left_join(., fws, by = c("indv" = "Indiv")) %>%
   left_join(., grc, by = c("indv" = "SampleId")) %>%
   left_join(., metadata) %>%
   select(-c(5:7)) %>% rename(lat = Admin.level1.latitude,  long = Admin.level1.longitude) %>% 
   # Create a new column to classify COI as Clonal (COI=1) or Polygenomic (COI>=2)
   mutate(McCOIL = ifelse(is.na(McCOIL), COI, McCOIL),
          coi_category = ifelse(COI == 1, "COI=1", "COI>=2"),
          Year2 = ifelse(Year2 == 1991, 1990, Year2),
          Year2 = ifelse(Year2 == 2013, 2014, Year2))

# writexl::write_xlsx(data, "data/metadata/complexity_of_infections.xlsx")

# 3. Plotting the scatter plot
plot_coi <- ggplot(data, aes(x = Admin.level.1, y = as.character(Year2), color = coi_category)) +
   geom_jitter(width = 0.25, height = 0.3, size = 2, alpha = 0.6) +  # Jitter to avoid overplotting
   scale_color_manual(values = c("COI=1" = "blue", "COI>=2" = "red")) +
   # scale_y_continuous(breaks = c(1984, 1990, 2001, 2008, 2014, 2015),
   #                    labels = as.character(c(1984, 1990, 2001, 2008, 2014, 2015))) +
   labs(x = "", y = "", color = "") +
   theme_minimal() +
   theme(axis.text.x = element_text(color = "black", face = "bold", size = 10, angle = 45, hjust = 1),
         axis.text.y = element_text(color = "black", face = "bold", size = 10),
         axis.line = element_line(color = "black"),
         axis.ticks = element_line(color = 'black', linewidth = .5),
         axis.ticks.length = unit(.2, "cm"),
         # axis.title = element_text(color = "black", face = "bold", size = 11),
         # plot.title = element_text(hjust = 0.5, face = "bold"),
         panel.grid = element_blank(),
         legend.text = element_text(color = "black", face = "bold", size = 12), #legend.position= c(0.15, 0.8) 
         )

print(plot_coi)

ggsave("results/figures/DistributionSamples_complexityOfInfection.tiff", 
       width = 150, height = 90, units = "mm", dpi = 600, plot = plot_coi)

# ================================
# Study sites levels
# ================================
coi.sites <- ggplot(data, aes(x = sites, y = McCOIL)) +
   geom_boxplot() +
   geom_jitter(width = 0.2, height = 0.2, alpha = 0.5, size = 2.5) +  # Jitter to avoid overplotting 
   labs(x = "", y = "Complexity of infections (COI)") +
   theme_minimal() +
   theme(axis.text.x = element_text(color = "black", face = "bold", size = 8, angle = 45, hjust = 1),
         axis.text.y = element_text(color = "black", face = "bold", size = 8),
         axis.line = element_line(color = "black", linewidth = 1, lineend = "square"),
         axis.title = element_text(color = "black", face = "bold", size = 12),
         axis.ticks = element_line(color = 'black', linewidth = 1),
         axis.ticks.length = unit(.2, "cm"))

ggsave("results/figures/coi_studySites.tiff", plot = coi.sites,
       width = 150, height = 90, units = "mm", dpi = 600)

# ================================
# Administrative levels
# ================================
sign_test <- compare_means(IBD ~ site, data = test_sign)

comparisons <- list(
   c("West Coast", "North Bank"),
   c("West Coast", "Lower River"),
   c("West Coast", "Central River"),
   c("West Coast", "Upper River"),
   c("North Bank", "Lower River"),
   c("North Bank", "Central River"),
   c("North Bank", "Upper River"),
   c("Lower River", "Central River"),
   c("Lower River", "Upper River"),
   c("Central River", "Upper River")
)

coi.adm <- data %>% 
   ggplot(aes(x = Admin.level.1, y = McCOIL)) +
   geom_boxplot() +
   geom_jitter(width = 0.2, height = 0.2, size = 2.5, alpha = 0.5) +  # Jitter to avoid overplotting 
   labs(x = "", y = "Complexity of infections (COI)")  +
   # stat_compare_means(method = "kruskal.test", label.y = max(data$McCOIL) + 0.5) +  # Kruskal-Wallis test
   stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons, 
                      vjust = 0.2, step_increase = 0, textsize = 4, tip.length = .02) +
   # geom_signif(comparisons = my_comparisons, annotations= sign_test$p.signif, step_increase = 0.1, textsize = 5) +
   scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
   theme_classic() +
   theme(axis.text.x = element_text(color = "black", face = "bold", size = 10, angle = 45, hjust = 1),
         axis.text.y = element_text(color = "black", face = "bold", size = 10),
         axis.line = element_line(color = "black", linewidth = 1, lineend = "square"),
         axis.title = element_text(color = "black", face = "bold", size = 11),
         axis.ticks = element_line(color = 'black', linewidth = 1),
         axis.ticks.length = unit(.25, "cm"),
         text = element_text(size = 14, face = "bold"))

print(coi.adm)

ggsave("results/figures/coi_admDiv.tiff", plot = coi.adm,
       width = 150, height = 90, units = "mm", dpi = 600)

# Statistical test
ks.test(data$McCOIL, "pnorm", mean(data$McCOIL), sd(data$McCOIL))

# Identify groups with variation
valid_groups <- data %>%
   group_by(Admin.level.1) %>%
   filter(sd(McCOIL) > 0)  # Keep only groups with nonzero standard deviation

# Run Shapiro-Wilk only on valid groups
normality_results <- valid_groups %>%
   group_by(Admin.level.1) %>%
   summarise(p_value = shapiro.test(McCOIL)$p.value)

print(normality_results)

# ================================
# ================================
data %>% 
   ggplot(aes(x = as.character(Year2), y = McCOIL, group = 1)) +
   # geom_boxplot() +
   geom_jitter(width = 0.3, height = 0.3, size = 2, alpha = 0.7) +  # Jitter to avoid overplotting 
   # scale_x_continuous(breaks = c(1984, 1990, 2001, 2008, 2014, 2015), expand = c(0.02,0),
   #                    labels = as.character(c(1984, 1990, 2001, 2008, 2014, 2015))) +
   labs(x = "", y = "Complexity of infections (COI)") +
   theme_minimal() +
   theme(axis.text.x = element_text(color = "black", face = "bold", size = 10), # , angle = 45, hjust = 1
         axis.text.y = element_text(color = "black", face = "bold", size = 10),
         axis.line = element_line(color = "black", linewidth = 1, lineend = "square"),
         axis.title = element_text(color = "black", face = "bold", size = 13),
         axis.ticks = element_line(color = 'black', linewidth = 1),
         axis.ticks.length = unit(.25, "cm"))

ggsave("results/figures/coi_byyear.tiff", dpi = 600, 
       width = 150, height = 90, units = "mm")

# ================================
# ================================
ggplot(data, aes(x = sites, y = 1-fws)) +
   geom_boxplot() +
   geom_jitter(width = 0.25, size = 1, alpha = 0.8) +  # height = 0.3, 
   scale_y_continuous(expand = c(0.05,0), limits = c(0,1)) + # , breaks = seq(0, 1, .1), labels = as.character(seq(0, 1, .1))
   labs(x = "", y = expression("1 - Fws")) +
   theme_classic() +
   theme(axis.text.x = element_text(color = "black", face = "bold", size = 8, angle = 45, hjust = 1),
         axis.text.y = element_text(color = "black", face = "bold", size = 8),
         axis.line = element_line(color = "black", linewidth = 1, lineend = "square"),
         axis.title = element_text(color = "black", face = "bold", size = 13),
         axis.ticks = element_line(color = 'black', linewidth = 1),
         axis.ticks.length = unit(.2, "cm"))

ggsave("results/figures/fws_studySites.tiff", dpi = 600, 
       width = 150, height = 90, units = "mm")

# ================================
# ================================
ggplot(data, aes(x = Admin.level.1, y = 1-fws)) +
   geom_boxplot() +
   geom_jitter(width = 0.25, size = 2, alpha = 0.8) +  # height = 0.3, 
   scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
   stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = comparisons, 
                      vjust = 0.2, step_increase = 0, textsize = 4, tip.length = .02) +
   labs(x = "", y = expression("1 - Fws")) +
   theme_classic() +
   theme(axis.text.x = element_text(color = "black", face = "bold", size = 8, angle = 45, hjust = 1),
         axis.text.y = element_text(color = "black", face = "bold", size = 8),
         axis.line = element_line(color = "black", linewidth = 1, lineend = "square"),
         axis.title = element_text(color = "black", face = "bold", size = 13),
         axis.ticks = element_line(color = 'black', linewidth = 1),
         axis.ticks.length = unit(.2, "cm"))

ggsave("results/figures/fws_admdiv.tiff", dpi = 600, 
       width = 150, height = 90, units = "mm")

# ================================
# ================================
ggplot(data, aes(x = as.character(Year2), y = 1-fws, group = 1)) +
   # geom_boxplot() +
   geom_jitter(width = 0.25, size = 1, alpha = 0.8) +  # height = 0.3, 
   scale_y_continuous(expand = c(0.02,0), limits = c(0,1)) + # , breaks = seq(0, 1, .1), labels = as.character(seq(0, 1, .1))
   labs(x = "", y = expression("1 - Fws")) +
   theme_classic() +
   theme(axis.text.x = element_text(color = "black", face = "bold", size = 10), # , angle = 45, hjust = 1
         axis.text.y = element_text(color = "black", face = "bold", size = 10),
         axis.line = element_line(color = "black", linewidth = 1, lineend = "square"),
         axis.title = element_text(color = "black", face = "bold", size = 13),
         axis.ticks = element_line(color = 'black', linewidth = 1),
         axis.ticks.length = unit(.2, "cm"))

ggsave("results/figures/fws_year.tiff", dpi = 600, 
       width = 150, height = 90, units = "mm")

# ================================
# ================================
ggplot(data, aes(x = Admin.level.1, y = fws, color = Admin.level.1)) +
   geom_boxplot(outlier.shape = NA, fill = NA, color = "black", width = 0.3, linewidth = .7,
                notch = TRUE, staplewidth = 0.5, varwidth = TRUE) +
   geom_jitter(width = 0.25, size = 2, alpha = 0.8) +  # height = 0.3, 
   # # Add error bars (mean ± standard error)
   # stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, linewidth = 1.2) +
   # 
   # # Add mean points
   # stat_summary(fun = mean, geom = "point", size = 3, shape = 18) +
   
   # Add boxplot with median and whiskers (IQR)
   # geom_boxplot(aes(fill = Admin.level.1), alpha = 0.3, outlier.shape = NA, width = 0.5) +
   # Add error bars (median line and IQR whiskers only)
   
   labs(x = "", y = expression("1 - Fws"), color = "") +
   theme_classic() +
   theme(axis.text.x = element_text(color = "black", face = "bold", size = 10, angle = 45, hjust = 1),
         axis.text.y = element_text(color = "black", face = "bold", size = 10),
         axis.line = element_line(color = "black", linewidth = 1, lineend = "square"),
         axis.title = element_text(color = "black", face = "bold", size = 13),
         axis.ticks = element_line(color = 'black', linewidth = 1),
         axis.ticks.length = unit(.2, "cm"),
         legend.position = "none")

ggsave("results/figures/fws_admdiv_v1.tiff", dpi = 600, 
       width = 150, height = 90, units = "mm")

# ================================
data %>% 
   ggplot(aes(x=1-fws)) +
   geom_histogram(fill="gray20", binwidth = .05) +
   scale_x_continuous(expand = c(0,0)) +
   scale_y_continuous(expand = c(0,0)) +
   labs(x = "Within-host genetic diversity (Fws)", y = "Number of samples") +
   theme_classic() +
   theme(legend.position = "none",
         axis.line = element_line(color = 'black', linewidth = 1),
         axis.text = element_text(size = 8, color = 'black', face = "bold"),
         axis.title = element_text(size = 10, color = 'black', face = "bold"),
         axis.ticks = element_line(color = 'black', linewidth = .9),
         axis.ticks.length = unit(.20, "cm")) 

ggsave("results/figures/fws_distribution.tiff", dpi = 600, 
       width = 150, height = 90, units = "mm")

data %>% 
   ggplot(aes(x=McCOIL)) +
   geom_histogram(fill="gray20", binwidth = 1) +
   scale_x_continuous(breaks = seq(1,3,1), labels = as.character(seq(1,3,1)), expand = c(0,0)) +
   scale_y_continuous(expand = c(0,0)) +
   labs(x = "Complexity of infection (COI)", y = "Number of samples") +
   theme_classic() +
   theme(legend.position = "none",
         axis.line = element_line(color = 'black', linewidth = 1),
         axis.text = element_text(size = 10, color = 'black', face = "bold"),
         axis.title = element_text(size = 12, color = 'black', face = "bold"),
         axis.ticks = element_line(color = 'black', linewidth = .9),
         axis.ticks.length = unit(.20, "cm")) 

ggsave("results/figures/coi_distribution.tiff", dpi = 600, 
       width = 150, height = 90, units = "mm")

df <- data %>% group_by(COI) %>% count()

# Basic piechart
ggplot(df, aes(x="", y=n, fill=COI)) +
   geom_bar(stat="identity", width=1, color="white") +
   coord_polar("y", start=0) +
   theme_void() + 
   theme(legend.position="none") +
   geom_text(aes(y = -1, label = n), color = "white", size=6) +
   scale_fill_continuous()




data <- readxl::read_xlsx("data/metadata/complexity_of_infections.xlsx")

data %>% 
   dplyr::mutate(coi = as.character(coi),
                 COI = as.character(COI),
                 McCOIL = as.character(McCOIL),
                 McCOIL_modif = as.character(McCOIL_modif)
   ) %>% 
   ggplot() + 
   geom_boxplot(aes(x = McCOIL_modif, y = fws, group = McCOIL_modif)) + 
   labs(x = "COI", y = "Fws") +
   theme_linedraw() +
   theme(legend.position = "none",
         axis.line = element_line(color = 'black'),
         axis.text = element_text(size = 8, color = 'black'),
         title = element_text(size = 10, color = 'black', face = "bold"))
