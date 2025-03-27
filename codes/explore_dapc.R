
# Load required packages
library(adegenet)
library(factoextra) # for visualization
library(tidyverse)

source("codes/functions.R")

set.seed(271088) # for reproducibility

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

# Replace missing data with the mean allele frequencies
df <- tab(data.genind, NA.method = "mean")

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

# Extract IDs from genind object & Remove duplicates
indiv <- adegenet::indNames(data.genind)

# Extract Species
pop <- metadata[which(metadata$samples %in% indiv),] %>% 
   pull(sites)

data_kmeans <- data.frame(sites = pop, grp_find$grp)
names(data_kmeans)
colnames(data_kmeans) = c("sites","K")
write.table(data_kmeans, "results/tables/Individuals_clusters_BIC_ind.txt", quote=F, row.names = T, col.names = T)

# Match geographic coordinates of each individual to its inferred genetic group in order to 
# see if there is any clustering relatively to the sampling location.
geo <-  metadata[which(metadata$samples %in% indiv),]
# kmean_geo <- merge(data_kmeans, geo, by = c("rownames(data_kmeans)" = "samples"))
kmean_geo <- data_kmeans %>% 
   rownames_to_column(., var = "samples") %>% 
   left_join(., geo, by = "samples")

write.table(kmean_geo, "results/tables/Individuals_clusters_locations.txt", 
            quote=F, row.names = F, col.names = T)

# Step 2: Run DAPC on the dataset using the optimal cluster number
opt_k <- grp_find$K  # Choose K with lowest BIC

dapc_results <- dapc(data.genind, pop = grp_find$grp, n.pca = 50, n.da = opt_k - 1)  # Retain 50 PCs initially
scatter(dapc_results)
saveRDS(dapc_results, "results/tables/dapc_results.rds")

# Step 3: Cross-validation to find the optimal number of PCs to retain
set.seed(270188)

cross_val <- xvalDapc(tab(data.genind, NA.method = "mean"), grp_find$grp, n.pca.max = 50, training.set = 0.9, n.rep = 1000)

optimal_pcs <- as.numeric(cross_val$`Number of PCs Achieving Lowest MSE`)  # Optimal PC number based on lowest RMSE

# Step 4: Perform final DAPC with the optimal number of PCs
dapc_final <- dapc(data.genind, pop = grp_find$grp, n.pca = optimal_pcs, n.da = opt_k - 1)

# Step 5: Visualization of DAPC results
scatter(dapc_final, posi.da = "bottomright", posi.pca = "topleft", cex = 1.5, cstar = 0.75)

# Alternative ggplot visualization
dapc_df <- as.data.frame(dapc_final$ind.coord)
dapc_df$Cluster <- as.factor(grp_find$grp)

cluster <- ggplot(dapc_df, aes(x = LD1, y = LD2, fill = Cluster)) +
   geom_hline(yintercept = 0, linewidth = 1) +
   geom_vline(xintercept = 0, linewidth = 1) +
   geom_point(size = 3, shape = 21, color = 'black', stroke = 1) +
   scale_fill_manual(values = c("#386CB0", "#FB8072", "#7FC97F"), name = "") +
   theme_bw() +
   labs(x = "Linear Discriminant 1", y = "Linear Discriminant 2") +
   theme(legend.position = "right",
         legend.text = element_text(colour = "black", size = 12, face = "bold"),
         axis.text = element_text(color = 'black', size = 10, face = 'bold'),
         axis.title = element_text(color = 'black', size = 12, face = 'bold'),
         axis.line = element_line(color = 'black', linewidth = 1, lineend = "square"),
         axis.ticks = element_line(color = 'black', linewidth = 1),
         axis.ticks.length = unit(.2, "cm"))

ggsave("results/figures/dapc_gambia_v1.tiff", plot = cluster,
       width = 110, height = 90, units = "mm", dpi = 600)

