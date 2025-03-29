


gmb <- sf::st_read( "../../01_data/meta/shapeFiles/geoBoundaries-GMB-ADM3_simplified.shp", quiet = TRUE) 

mtdtclst <- metadata %>% 
   distinct(sites, .keep_all = TRUE) %>% 
   select(-c(indv, locations, year, Admin.level.1, Country.latitude, Country.longitude, Year2)) %>% 
   rename(lat = Admin.level1.latitude, long = Admin.level1.longitude)

hrp.btw.sites <- calculate_relatedness_stats(ibd_metadata, threshold = 0.8, 
                                         filter_columns = c("sites.1", "sites.2"), 
                                         group_column = "sites.1",
                                         filter_operator = "!=") %>% 
   
   left_join(mtdtclst, by = c("sites.1" = "sites")) %>%
   left_join(mtdtclst, by = c("sites.2" = "sites"), suffix = c(".1", ".2")) %>% 
   filter(HighlyRelatedFraction != 0)

# =========
# Make plot
# =========
plot.map <- ggplot() +
   ggplot2::geom_sf(data = gmb, fill = "gray90", col = NA) +
   ggplot2::geom_sf(data = gmb, fill = NA, col = "gray20", linewidth = .5) +
   ggplot2::geom_sf(data = gmb, fill = NA, col = "black", linewidth = .8) +
   ggplot2::labs(x = "Longitude", y = "Latitude") +
   theme_bw()

print(plot.map)

mapplot <- plot.map +
      # Connectivity between study sites
   geom_curve(data = hrp.btw.sites, lineend = "round", color = '#330099', curvature = .2, # '#8F7700FF'
              aes(x = long.1, y = lat.1, xend = long.2, yend = lat.2, linewidth = HighlyRelatedPairs)) + # show.legend = TRUE, 
   scale_linewidth_continuous(name = "Fraction") +
      # Add non-overlapping text labels
   ggrepel::geom_label_repel(data = mtdtclst, aes(x = long, y = lat, label = sites),
                             size = 4, fontface = "bold", box.padding = unit(1.5, "lines"),
                             angle = 45, #box.padding = 0.5,  # Increases space around labels
                             point.padding = 0.3,  # Adjusts label-to-point distance
                             max.overlaps = Inf,  # Ensures all labels are placed
                             force = 5,  # Adjusts the repelling strength
                             direction = "both") + # Allows movement in all directions
   
   geom_point(data = mtdtclst, aes(x = long, y = lat), fill = "steelblue1", size = 4,
              color = "black", show.legend = FALSE, shape = 21, stroke = 1) +
   
   # scale_color_viridis_c("Pairwise IBD", option="C", direction = 1, values = c(0,1)) +
   theme(axis.line = element_line(color = 'black', linewidth = .7, lineend = 'square'),
         axis.text = element_text(size = 10, color = 'black', face = "bold"),
         axis.title = element_text(size = 14, color = 'black', face = "bold"),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         axis.ticks = element_line(color = 'black', linewidth = 1),
         axis.ticks.length = unit(.2, "cm"),
         panel.background = element_rect(fill = "gray60"), 
         legend.text = element_text(size = 10, color = 'black', face = "bold"),
         legend.title = element_text(size = 12, color = 'black', face = "bold"))

ggsave("results/figures/map_of_fraction_hrp.tiff", width = 310, 
       height = 190, units = "mm", dpi = 600, plot = mapplot)
