
filter_vcf <- function(input, data.out){
   
   # Annotate VCF by modifying ID column
   system(paste0("bcftools annotate --set-id +'%CHROM:%POS' ", input, " > gambia.vcf"))
   
   #==================
   # Quality Control
   #==================
   geno = 0.05; maf = 0.01; mind = 0.05; 
   LD.window = 50; LD.step = 10; LD.correlation = 0.2
   
   system(paste0("plink --vcf  ", input,
                 " --geno ", geno,
                 " --mind ", mind,
                 " --maf ", maf,
                 " --make-bed --allow-extra-chr",
                 " --out ", data.out, "_qc"))
   
   system(paste0("plink --bfile ", data.out, 
                 "_qc -indep-pairwise ", LD.window, " ",
                 LD.step, " ", LD.correlation,
                 " --allow-extra-chr",
                 " --out ", data.out))
   
   system(paste0("plink --bfile ", data.out, 
                 "_qc --extract ", data.out, ".prune.in",
                 " --recode vcf-iid --allow-extra-chr", 
                 " --out ", data.out))
   
   system(paste0("rm -rf dapc_* dapc.prune* dapc.nosex"))
}

vcf_genind <- function(input_vcf, metadata){
   
   #===================
   # Read filtered VCF
   #===================
   data <- vcfR::read.vcfR(input_vcf, verbose = FALSE)
   
   #===================
   # Convert VCF 
   # to a genind object
   #===================
   data.genind <- vcfR::vcfR2genind(data, return.alleles = TRUE, NA.char = "./.") 
   
   # Extract IDs from genind object & Remove duplicates
   indiv <- adegenet::indNames(data.genind)
   
   # Extract Species
   pop <- metadata[which(metadata$samples %in% indiv),] %>% 
      pull(sites)
   
   adegenet::pop(data.genind) <- as.factor(pop)
   
   return(data.genind)
}


pairwise.fst <- function(x){
   
   gene_fst <- hierfstat::genet.dist(x, method = "WC84")
   
   # Convert dist object to data.frame
   fst.matrix <- as.matrix(gene_fst)
   
   ## Sort column names
   fst.matrix <- fst.matrix[order(rownames(fst.matrix)), order(colnames(fst.matrix))]
   ind <- which( upper.tri(fst.matrix), arr.ind = TRUE)
   fst.df <- data.frame(Site1 = dimnames(fst.matrix)[[2]][ind[,2]],
                        Site2 = dimnames(fst.matrix)[[1]][ind[,1]],
                        Fst = fst.matrix[ ind ] %>% round(digits = 3))
   
   # Convert minus values to zero
   fst.df$Fst[fst.df$Fst < 0] <- 0
   
   return(fst.df)
}

plot.fst <- function(data){
   
   # Fst italic label
   fst.label <- expression(italic("F")[ST])
   
   # Extract middle Fst value for gradient argument
   mid <- max(data$Fst) / 2
   
   # Plot heatmap
   data %>%
      ggplot(aes(x = Site1, y = Site2, fill = Fst)) +
      geom_tile(colour = "black") +
      geom_text(aes(label = Fst), color="black", size = 4) +
      scale_fill_gradient2(low = "blue", mid = "pink", high = "red",
                           midpoint = mid, name = fst.label, 
                           limits = c(0, max(data$Fst))) + 
      # scale_x_discrete(labels = c(glue("*An. coluzzii*"), glue("*An. gambiae*"), 
      #                             glue("*Hybrid C-G*")), expand = c(0,0)) +
      scale_y_discrete( expand = c(0,0), position = "right") +
      theme(legend.position = c(0.15, 0.7),
            legend.title = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 10, face = "bold", colour = "#000000"),
            axis.text.x = element_text(size = 10, vjust = 0.5, face = "bold", color = "#000000"),
            axis.text.y = element_text(size = 10, face = "bold", color = "#000000"),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_blank())
}


# Summarise IBD data

summary_data <- function(data, thres = 0.25) {
   # Ensure the Ibs threshold is between 0 and 1
   if (thres < 0 || thres > 1) {
      stop("thresreshold must be a value between 0 and 1.")
   }
   
   # Dynamically create column names based on the Ibs threshold
   threshold_col <- paste0("PairCount.", as.character(thres))
   threshold_group_per_col <- paste0("TotalPairCount.", as.character(thres))
   threshold_location_per_col <- paste0("LocationPairCount.", as.character(thres))
   
   # Summarize the data using the melted IBD data frame
   # pairs_data <- data %>%
   #    dplyr::group_by(sites.1, sites.2) %>%
   #    dplyr::summarise(
   #       mean_IBD = round(mean(IBD), 2),
   #       median_IBD = round(stats::median(IBD), 2),
   #       TotalPairCount = dplyr::n(),
   #       !!threshold_col := sum(IBD >= thres),
   #       .groups = "drop") %>%
   #    dplyr::mutate(sites.1 = as.character(sites.1), sites.2 = as.character(sites.2),
   #                  pair_location = ifelse(sites.1 < sites.2,
   #                                         paste(sites.1, sites.2, sep = " - "),
   #                                         paste(sites.2, sites.1, sep = " - "))) %>% 
   #    dplyr::ungroup() %>% 
   
   pairs_data <- data %>%
      select(c(1:3, 5, 14)) %>% 
      dplyr::mutate(sites.1 = as.character(sites.1), 
                    sites.2 = as.character(sites.2),
                    pair_location = ifelse(sites.1 < sites.2, 
                                           paste(sites.1, sites.2, sep = " - "), 
                                           paste(sites.2, sites.1, sep = " - "))) %>% 
      left_join(., coordinates, by = c("sites.1" = "sites")) %>% 
      left_join(coordinates, by = c("sites.2" = "sites"), suffix = c(".1", ".2")) %>%
      dplyr::group_by(pair_location, lat.1, long.1, lat.2, long.2) %>%
      dplyr::summarise(
         mean_IBD = round(mean(IBD), 2),
         median_IBD = round(stats::median(IBD), 2),
         TotalPairCount = dplyr::n(),
         !!threshold_col := sum(IBD >= thres),  # Count how many pairs exceed threshold
         .groups = "drop") %>% 
      ungroup()

   
   total_threshold_pairs <- sum(pairs_data[[threshold_col]])
   
   pairs_data <- pairs_data %>%
      dplyr::mutate(
         !!threshold_group_per_col := round(.data[[threshold_col]] / total_threshold_pairs * 100, 2), # Group percentage
         !!threshold_location_per_col := round(.data[[threshold_col]] / TotalPairCount * 100, 2) # Location-specific percentage
      )
   
   return(pairs_data)
}



# # ==================================================================
# # Also, you can use K-means clustering to find the K with minimum BIC
# 
# maxK <- 10
# myMat <- matrix(nrow=10, ncol=maxK)
# colnames(myMat) <- 1:ncol(myMat)
# for(i in 1:nrow(myMat)){
#    grp <- find.clusters(popdata, n.pca = 40, choose.n.clust = FALSE,  max.n.clust = maxK)
#    myMat[i,] <- grp$Kstat
# }
# 
# library(reshape2)
# my_df <- melt(myMat)
# colnames(my_df)[1:3] <- c("Group", "K", "BIC")
# my_df$K <- as.factor(my_df$K)
# head(my_df)
# 
# # The K vs. BIC plot
# library(ggplot2)
# p1 <- ggplot(my_df, aes(x = K, y = BIC))
# p1 <- p1 + geom_boxplot()
# p1 <- p1 + theme_bw()
# p1 <- p1 + xlab("Number of groups (K)")
# p1
# 
# 
# ### conduct DAPC
# vcf_dapc <- dapc(popdata, n.pca = 40, n.da = 5)
# 
# ### IMPORTANT: find the optimal # of PCs to retain using the 'a-score' criterion
# temp <- optim.a.score(vcf_dapc)
# 
# #* n = 6 PCs = optimal!
# 
# # Cross validation method to find # of PCs
# xval = xvalDapc(x, popdata$pop, n.pca.max=10, training.set=0.9,
#                 result="groupMean", center=TRUE, scale=FALSE,
#                 n.rep=100, n.pca=NULL, parallel="snow", ncpus=4)
# 
# ## Number of PCs with best stats
# xval$`Number of PCs Achieving Highest Mean Success`
# 
# 
# ### DAPC1 vs 2 scatter
# scatter(vcf_dapc, col = mycol, xax = 1, yax = 2, cex = 2, scree.da=FALSE, legend = FALSE, grp = pop(popdata))
# 
# ### DAPC1 vs 3 scatter
# scatter(vcf_dapc, col = mycol, xax = 1, yax = 3, cex = 2, scree.da=FALSE, legend = FALSE, grp = pop(popdata))
# 
# 
# ### Check the assignment plot
# assignplot(vcf_dapc)
# 
# ### plot membership probabilities
# compoplot(vcf_dapc, lab="",ncol=1,col=mycol)

counts <- function(data, column = "sites"){
   data %>% 
      group_by(across(all_of(column))) %>% 
      count()
}

# Define the function for calculating highly related pairs and their confidence intervals
calculate_relatedness_stats <- function(data, threshold = 0.25, filter_columns, group_column, filter_operator = "==") {
   library(Hmisc)
   
   # # Create a filter expression dynamically
   # filter_expr <- paste0("get(filter_columns[1]) ", filter_operator, " get(filter_columns[2])")
   # 
   # # Apply the filter
   # filtered_data <- data %>%
   #    filter(eval(parse(text = filter_expr))) %>%
   #    mutate(HighlyRelated = ifelse(IBD > threshold, 1, 0))
   # 
   # # Determine the grouping variable dynamically
   # if (filter_operator == "==") {
   #    grouped_data <- filtered_data %>%
   #       group_by(across(all_of(group_column)))
   # } else {
   #    grouped_data <- filtered_data %>%
   #       group_by(across(all_of(filter_columns)))  # Groups by both columns if operator is not "=="
   # }
   
   # Summarization step (same for both cases)
   grouped_data <- data %>%
      dplyr::mutate(HighlyRelated = ifelse(IBD > threshold, 1, 0)) %>% 
      dplyr::group_by(across(all_of(filter_columns))) %>%
      summarise(
         sampleSize = n(),
         MeanIBD = mean(IBD),
         sdIBD = sd(IBD, na.rm = TRUE),
         seIBD = sdIBD / sqrt(n()),
         U95CI = MeanIBD + 1.96 * seIBD, 
         L95CI = MeanIBD - 1.96 * seIBD,
         nbre_hrp = sum(HighlyRelated),
         fraction = (mean(IBD > threshold)),
         SE = sqrt(fraction * (1 - fraction) / sampleSize), 
         ymin = fraction - SE, 
         ymax = fraction + SE,
         
         # Wilson Score Interval
         CI_Lower = binconf(fraction * sampleSize, sampleSize, method = "wilson")[2],  
         CI_Upper = binconf(fraction * sampleSize, sampleSize, method = "wilson")[3],   
         
         # Clopper-Pearson Exact Interval
         lower_ci = binconf(fraction * sampleSize, sampleSize, method = "exact")[2],  
         upper_ci = binconf(fraction * sampleSize, sampleSize, method = "exact")[3],   
         
         # Bayesian Shrinkage (Empirical Bayes)
         # corectedFraction = EbayesThresh::ebayesthresh(HighlyRelatedRaw / SampleSize),
         correctedFraction = ifelse(fraction != 0, EbayesThresh::ebayesthresh(fraction, 
                                                                              prior = "laplace",
                                                                              bayesfac = FALSE), 
                                    NA),
         
         .groups = "drop")
   
   total_hrp <- sum(grouped_data$nbre_hrp)
   
   grouped_data <- grouped_data %>%
      dplyr::mutate(
         fraction_within_hrp = round((nbre_hrp / total_hrp)*100, 2), # Group percentage
         fraction_hrp = round((nbre_hrp / sampleSize)*100, 2) # Location-specific percentage
      )
   
   return(grouped_data)
}


# ========================================================
recap_v1 <- function(ibd.data){
   
   household <- ibd.data %>%
      filter(Household == "Within household") %>%
      mutate(site = "Household") %>%
      select(p1, p2, IBD, site)
   
   compound <- ibd.data %>%
      filter(Compound == "Within Compound") %>%
      mutate(site = "Compound") %>%
      select(p1, p2, IBD, site)
   
   village <- ibd.data %>%
      filter(Village == "Within Village") %>%
      mutate(site = "Village") %>%
      select(p1, p2, IBD, site)
   
   region <- ibd.data %>%
      filter(Region == "Within Region") %>%
      mutate(site = "Region") %>%
      select(p1, p2, IBD, site)
   
   results <- rbind.data.frame(household, compound, village, region)
   return(results)
}











