
# https://github.com/laurabenestan/DAPC
# https://rpubs.com/cfb0001/777801

library(hierfstat)
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(tidyverse)
library(adegenet)
library(reshape2)
library(ggpubr)

fig.format <- c("pdf", "tiff")

source("codes/functions.R")

# Load metadata file
metadata <- read_tsv("data/metadata/metadata_to_use.tsv") %>%
   select(-c(2, 3, 7:9, 16:24)) %>%
   mutate(Admin.level.1 = case_when(sites %in% c("Bessi", "Brikama", "Serrekunda") ~ 'West Coast', 
                                    sites %in% c("Chogen", "Yallal Ba") ~ 'North Bank',
                                    sites %in% c("Sinchu", "Dongoro Ba") ~ 'Lower River',
                                    sites %in% c("Farafenni", "Sare Seedy" ) ~ 'Central River',
                                    sites %in% c("Basse", "Njayel", "Madina Samako", 
                                                 "Sare Wuro", "Gunjur Kunta") ~ 'Upper River'),
          Admin.level.1 = factor(Admin.level.1, levels = c("West Coast", "North Bank", "Lower River",
                                                           "Central River", "Upper River")))

# input <- "data/reference/gambia.vcf"

# Create output file name
output <- "dapc"

# filter_vcf(input, output)

data.genind <- vcf_genind(paste0(output, ".vcf"), metadata)

# #==========================
# # Compute pairwise FST 
# # (Weir & Cockerham 1984)
# #==========================
# fst <- pairwise.fst(data.genind)
# 
# fst.plot <- plot.fst(fst)
# 
# # Save plot
# for (i in fig.format) {
#    ggsave(sprintf("results/figures/pop-differentiation.%s", i), plot = fst.plot,
#           units = "mm", width = 170, height = 130, dpi = 600)
# }

#================
# Perform DAPC
#================
set.seed(271088)

# Replace missing data with the mean allele frequencies
df <- tab(data.genind, NA.method = "mean")

#==============
# Perform PCA
#==============

pca <- dudi.pca(df, scannf = FALSE, scale = TRUE, nf = 10)

# Analyse how much percent of genetic variance is explained by each PC
percent <- pca$eig/sum(pca$eig)*100

barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)",
        ylim = c(0,2), names.arg = round(percent, 1))

# Create a data.frame containing individual coordinates
pca_coords <- as.data.frame(pca$li)

#=======================
# Visualize PCA results
#=======================

# Add a column containing individuals
pca_coords$Ind <- rownames(pca_coords)

# Add a column with the site IDs
pca_coords$Site <- factor(data.genind$pop, levels = c("Serrekunda", "Brikama", "Bessi", "Chogen", "Yallal Ba",
                                                      "Sinchu", "Dongoro Ba", "Farafenni", "Sare Seedy",
                                                      "Njayel", "Madina Samako", "Basse", "Sare Wuro", "Gunjur Kunta"))

# Define colour palette
cols <- c(brewer.pal(nPop(data.genind), "Set3"), brewer.pal(nPop(data.genind), "Accent"))

# colors <- sample(cols, nPop(data.genind), replace = FALSE)

colors <- c("#FFFFB3", "#CCEBC5", "#BC80BD", "#386CB0", "#80B1D3", "#FCCDE5", "#FDC086", "#FB8072",
            "#FDB462", "#F0027F", "#BF5B17", "#8DD3C7", "#7FC97F", "#D9D9D9")

axes <- colnames(pca_coords)[1:10]

# Custom x and y labels
xlab <- paste(axes[1], " (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab <- paste(axes[2], " (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Custom theme for ggplot2
ggtheme <- theme(axis.text = element_text(colour = "black", size = 12, face = "bold"),
                 axis.title = element_text(colour = "black", size = 14, face = "bold"),
                 axis.line = element_line(linewidth = 1),
                 panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
                 panel.background = element_blank(),
                 axis.ticks = element_line(color = 'black', linewidth = 1),
                 axis.ticks.length = unit(.25, "cm"),
                 legend.text = element_text(colour = "black", size = 12, face = "bold"),
                 legend.title = element_text(colour = "black", size = 14, face = "bold"))

# Scatter plot PC 1 vs. 2
plot1 <- ggplot(data = pca_coords, aes(x = Axis1, y = Axis2)) +
   # points
   geom_point(aes(fill = Site), shape = 21, size = 5, stroke = 1) +
   geom_hline(yintercept = 0) +
   geom_vline(xintercept = 0) +
   # colouring
   scale_fill_manual(values = colors, name = "") +
   # scale_colour_manual(values = cols) +
   # custom labels
   labs(x = xlab, y = ylab) + ggtheme

plot2 <- ggplot(data = pca_coords, aes(x = Axis1, y = Axis3)) +
   # points
   geom_point(aes(fill = Site), shape = 21, size = 5, color = "black", stroke = 1) +
   geom_hline(yintercept = 0) +
   geom_vline(xintercept = 0) +
   # colouring
   scale_fill_manual(values = colors, name = "") +
   # scale_colour_manual(values = cols) +
   # custom labels
   labs(x = xlab, y = ylab) + ggtheme

plot3 <- ggplot(data = pca_coords, aes(x = Axis2, y = Axis4)) +
   # points
   geom_point(aes(fill = Site), shape = 21, size = 5, color = "black", stroke = 1) +
   geom_hline(yintercept = 0) +
   geom_vline(xintercept = 0) +
   # colouring
   scale_fill_manual(values = colors, name = "") +
   # scale_colour_manual(values = cols) +
   # custom labels
   labs(x = xlab, y = ylab) + ggtheme

print(plot3)


# Step 1: Detecting the optimal number of clusters using K-means
grp_find <- find.clusters(data.genind, max.n.clust = 40)  # Test up to 10 clusters

# Visualize BIC values to determine the best number of clusters
plot(grp_find$Kstat, type="b", xlab="Number of clusters (K)", ylab="BIC", main="Optimal K Selection")

# Step 2: Run DAPC on the dataset using the optimal cluster number
opt_k <- grp_find$K  # Choose K with lowest BIC
dapc_results <- dapc(data.genind, pop = grp_find$grp, n.pca = 50, n.da = opt_k - 1)  # Retain 50 PCs initially

saveRDS(dapc_results, "results/tables/dapc_results_k3.rds")

scatter(dapc_results)


# Perform cross validation to find the optimal number of PCs to retain in DAPC
crossval <- xvalDapc(df, data.genind$pop, result = "groupMean", xval.plot = TRUE)

# Number of PCs with best stats (lower score = better)
numPCs <- as.numeric(crossval$`Number of PCs Achieving Lowest MSE`)

# Run a DAPC using population IDs as priors
dapc1 <- dapc(data.genind, data.genind$pop, n.pca = numPCs, n.da = NULL)
scatter(dapc1)

# Perform DAPC analysis
dapc_res <- dapc(data.genind, grp_find$grp)  # Perform DAPC
scatter(dapc_res)

# # Extract DAPC results
# dapc_df <- as.data.frame(dapc_res$ind.coord)  # DAPC coordinates
# dapc_df$group <- grp_find$grp  # Add group information
# dapc_df$sites <- factor(data.genind$pop, levels = sort(unique(metadata$sites)))  # Add species information
# 
# # Add a column containing individuals
# dapc_df$samples <- rownames(dapc_df)
# 
# # Analyse how much percent of genetic variance is explained by each axis
# percent = dapc_res$eig/sum(dapc_res$eig)*100
# 
# # make axis labels
# x_lab <- sprintf("Axis 1 (%s%%)", round(percent[1], digits = 1))
# y_lab <- sprintf("Axis 2 (%s%%)", round(percent[2], digits = 1))
# z_lab <- sprintf("Axis 3 (%s%%)", round(percent[3], digits = 1))
# 
# 
# # produce 2D scatterplot
# plot1 <- ggplot(data = dapc_df, aes(x = LD1, y = LD2)) + 
#    theme_bw(base_size = 8) +
#    # stat_ellipse(type = "t", show.legend = FALSE, aes(color = sites)) + # level = .75, 
#    geom_point(aes(fill = sites), shape = 21, size = 5, color = "black", stroke = 1) + 
#    xlab(x_lab) + ylab(y_lab) +
#    scale_fill_manual(values = colors, name = "Study sites") + # , guide = "none"
#    theme(legend.title = element_text(size = 14, face = "bold"), # , hjust = 0.5
#          legend.text = element_text(size = 12, face = "bold", colour = "#000000"),
#          axis.text = element_text(size = 12, face = "bold", color = "#000000"),
#          axis.title = element_text(size = 14, face = "bold", color = "#000000"))
# 
# print(plot1)
# 
# plot2 <- ggplot(data = dapc_df, aes(x = LD1, y = LD3)) +
#    geom_point(aes(fill = sites), shape = 21, size = 5, color = "black", stroke = 1) + 
#    stat_ellipse(type = "t", show.legend = FALSE) + # level = .75, 
#    xlab(x_lab) + ylab(z_lab) + 
#    theme_bw(base_size = 8) +
#    scale_fill_manual(values = colors, name = "Study sites") + 
#    theme(legend.title = element_text(size = 14, face = "bold"), 
#          legend.text = element_text(size = 12, face = "bold", colour = "#000000"),
#          axis.text = element_text(size = 12, face = "bold", color = "#000000"),
#          axis.title = element_text(size = 14, face = "bold", color = "#000000"))
# 
# plot3 <- ggplot(data = dapc_df, aes(x = LD2, y = LD3)) +
#    geom_point(aes(fill = sites), shape = 21, size = 5, color = "black", stroke = 1) + 
#    stat_ellipse(type = "t", show.legend = FALSE) +
#    xlab(y_lab) + ylab(z_lab) +
#    theme_bw(base_size = 8) +
#    scale_fill_manual(values = colors, name = "Study sites") + 
#    theme(legend.title = element_text(size = 14, face = "bold"), 
#          legend.text = element_text(size = 12, face = "bold", colour = "#000000"),
#          axis.text = element_text(size = 12, face = "bold", color = "#000000"),
#          axis.title = element_text(size = 14, face = "bold", color = "#000000"))
# 
# 
# # produce 3D scatterplot
# plot4 <- bobfunctions2::gg3d_scatterplot(dapc_df$LD1, dapc_df$LD2, 
#                                          dapc_df$LD3, colour = dapc_df$sites, 
#                                          size = 5, theta = 125, d = 2.0, phi = 5, axis_on = TRUE,
#                                          x_lab = "", y_lab = "", z_lab = "", tick_length = 0.5,
#                                          axis_lab_dist = 5, axis_lab_size = 3,
#                                          grid_size = 0.15, zero_line_size = 0.6) +
#    scale_color_manual(values = colors, name = "Study sites") +
#    theme(legend.title = element_text(size = 14, face = "bold", hjust = 0.5),
#          legend.text = element_text(size = 10, face = "bold", colour = "#000000"),
#          axis.text = element_text(size = 14, face = "bold", color = "#000000"),
#          axis.title = element_text(size = 12, face = "bold", color = "#000000"))
# 
# print(plot4)
# 
# ggsave(paste0("results/new_data/", genename, "_3D.pdf"),
#        plot = plot4, width = 130, height = 130, units = "mm")


