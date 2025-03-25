
library(tidyverse)
library(patchwork)

# Function
calc_biallelic_heterozygosity <- function(p_alt) {
   return (2 * p_alt * (1 - p_alt))
}

fig.format <- c("pdf", "tiff")


# Format data to meet dcifer requirements
input <- "data/reference/senegam_pf.vcf"

# Load metadata file
metadata <- read_tsv("data/metadata/metadata_to_use.tsv")

# Save samples from the Gambia
samples <-  "data/metadata/gambian.isolates.txt"

metadata %>% 
   pull(samples) %>% 
   write.table(., samples, col.names = FALSE,
               row.names = FALSE, quote = FALSE)

# Extract samples from the VCF
output <- "data/reference/gambia.vcf"

system(paste0("bcftools view -S ", samples, " -o ", output, " ", input))

# Read Output VCF
vcf <- vcfR::read.vcfR(output)

# extract information
loci <- vcf@fix[,1:2] %>%
   tibble::as_tibble() %>%
   dplyr::mutate(POS = as.numeric(POS), Key = 1:dplyr::n()) %>%
   dplyr::select(c("CHROM", "POS", "Key"))

# tidy up the DRC data to long format
geno <- vcfR::extract_gt_tidy(vcf) %>% 
   dplyr::select(-c(gt_GQ, gt_MIN_DP, gt_PGT, gt_PID, gt_PL, gt_RGQ, gt_SB))


# Format data into one long data frame
combined_long <- geno %>%
   # now lets merge the loci information with the individual level information
   dplyr::full_join(x = loci, y = ., by  = "Key") %>%
   # don't need Key anymore
   dplyr::select(-c("Key")) %>% 
   mutate(gt_GT = case_when(gt_GT == "0/0" ~ 0,
                            gt_GT == "0/1" ~ 0.5,
                            gt_GT == "1/1" ~ 1)) %>% 
   drop_na(gt_GT) %>% 
   arrange(Indiv)

# Estimate number of heterozygous sites
n_hets <- combined_long %>% 
   dplyr::group_by(Indiv) %>% 
   dplyr::summarise(
      # multiple by 2 to assume diploid genome
      n_sites = length(gt_GT),
      n_hets = sum(gt_GT == 0.5, na.rm=T))

# Compute population allele frequencies
plaf <- combined_long %>% 
   dplyr::group_by(CHROM, POS) %>% 
   dplyr::summarise(
      # multiple by 2 to assume diploid genome
      PLAF = sum(gt_GT * 2, na.rm = T) / (2* sum(!is.na(gt_GT))),
      .groups = "drop")

# calculate het for pop
plaf <- plaf %>% 
   dplyr::mutate(
      het_plaf = purrr::map_dbl(PLAF, calc_biallelic_heterozygosity))


# extract GT information and get WSAF for every sample at each loci
combined_long <- combined_long %>% 
   # lets make some new variables
   dplyr::mutate(
      rad = purrr::map_dbl(gt_AD, function(x){vcfR::masplit(
         as.matrix(x), record = 1, sort = FALSE, decreasing = FALSE)}),
      # get alternate allele depth 
      aad = purrr::map_dbl(gt_AD, function(x){vcfR::masplit(
         as.matrix(x), record = 2, sort = FALSE, decreasing = FALSE)}),
      # calculate within-sample reference allele freq
      wsaf = aad/(rad + aad),
      wsaf = ifelse(is.nan(wsaf), NA, wsaf)) %>% # occurs when 0/0
    
   # now let's select the variables that we want
   dplyr::select(c("CHROM", "POS", "Indiv", "gt_GT", "wsaf")) 


combined_long <- combined_long %>% 
   # lets make some new variables
   dplyr::mutate(
      het_wsaf = calc_biallelic_heterozygosity(wsaf))

# now calculate fws
Fws <- dplyr::full_join(combined_long, plaf, by = c("CHROM", "POS")) %>% 
   dplyr::mutate(Fwsloci = het_wsaf / het_plaf) %>% 
   dplyr::group_by(Indiv) %>% 
   dplyr::summarise(
      fws = mean( 1 - het_wsaf / het_plaf, na.rm = T),
      .groups = "drop")


Fws <- dplyr::left_join(Fws, n_hets) %>% 
   relocate(n_sites, .before = fws)

# Save Fws results
write_tsv(Fws, "results/tables/Fws.tsv", col_names = TRUE)


# Plot
plot1 <- Fws %>% 
   ggplot() +
   geom_histogram(aes(x = fws, y = (after_stat(count)/sum(after_stat(count)))*100),
                  color = "#000000", fill = "#d9d9d9", linewidth = .8) +
   labs(x = "Within-host genetic diversity (Fws)", y = "Frequency (%)") +
   scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
   theme_minimal() +
   theme(legend.position = "none",
         axis.line = element_line(color = 'black', linewidth = 1),
         axis.text = element_text(size = 11, color = 'black'),
         axis.title = element_text(size = 14, color = 'black', face = "bold"),
         axis.ticks = element_line(color = 'black', linewidth = .9),
         axis.ticks.length = unit(.20, "cm")) 

# Save plot
for (i in fig.format) {
   ggsave(sprintf("results/figures/figure1_fws.%s", i), width = 110, 
          height = 90, units = "mm", dpi = 600, plot = plot1)
}


