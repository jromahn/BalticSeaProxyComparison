figure_theme <-theme_light() +theme(axis.title = element_text(hjust=1, size= 20), 
                                    axis.text = element_text(size = 18), 
                                    axis.text.x = element_text(angle = 45,hjust=1), 
                                    axis.line = element_line(size = 0.2),
                                    plot.title = element_text(face = "bold", size = 22 ),
                                    plot.subtitle = element_text(size = 20),
                                    strip.text = element_text(size = 18), # facet_wrap
                                    panel.grid.major.x = element_line(colour = "grey28", size = 0.08),
                                    panel.grid.major.y = element_line(colour = "grey28", size = 0.1), 
                                    #panel.grid.minor = element_blank(),
                                    legend.position="bottom",
                                    legend.text = element_text(size=16), 
                                    legend.title = element_text(size=18))

## calculate stats per replicate
# data.frame to calc stats from and how to name it in df
# col - samples, rows - asvs
community_stats <- function(input_df, title){
  stats <- data.frame(type= title, sample_no=length(colnames(input_df)) ,
                      total.reads = sum(rowSums(input_df)), total.asvs= sum(1*(rowSums(input_df)>0)), replicates = sum(1*(colSums(input_df)>0)),
                      mean.reads = mean(rowSums(input_df)), median.reads= median(rowSums(input_df)), 
                      sd.reads =sd(rowSums(input_df)), max.reads = max(rowSums(input_df)), min.reads= min(rowSums(input_df)),
                      mean.asvs = mean(rowSums(1*(input_df>0))), median.asvs = median(rowSums(1*(input_df>0))), 
                      sd.asvs =sd(rowSums(1*(input_df>0))), max.asvs= max(rowSums(1*(input_df>0))), min.asvs = min(rowSums(1*(input_df>0))))
  
  return(stats)
}

# ordination plot of communtiy data (col <- asvs, rows <- samples) 
#input: nmds data, original data with more information, 2 columns within data, merging column (replicate)--> rownames (before )!!!!
# column defines after which column the color is chosen
# column defines after which column the shape is chosen
# before doing metaDMs rownames have to be the same to column after which it should be merged


plot_simpleNDMS_color <- function(NMDS_result, data, column_c,merge_by){
  replicate_coordinates <- data.frame(NMDS_result["points"])
  species_coords <- data.frame(NMDS_result["species"])
  
  replicate_coordinates$replicate <- rownames(replicate_coordinates)
  replicate_coordinates <- merge(replicate_coordinates, data, by =merge_by, all.x = TRUE)
  #criteria color
  colnames(replicate_coordinates)[colnames(replicate_coordinates) == column_c] <- 'crit_color'
  replicate_coordinates$crit_color <- as.factor(replicate_coordinates$crit_color )
  
  # Now we can build the plot! 
  
  p <- ggplot() +
    geom_point(data = replicate_coordinates, aes(x = points.MDS1, y = points.MDS2, 
                                                 color = crit_color), size=2) + 
    annotate(geom = "label", 
             x = min(replicate_coordinates$points.MDS1)+1, 
             y =min(replicate_coordinates$points.MDS2), 
             size =4,
             label = paste("Stress: ", round(NMDS_result$stress, digits = 3)))+ #+
    labs(color = column_c) + theme( axis.text = element_text(size = 12)) + theme_light()
  
  list_data <- list(replicate_coordinates, p)
  return(list_data)
}

plot_simpleNDMS_colored_symbol <- function(NMDS_result, data, column_c,column_s, merge_by){
  replicate_coordinates <- data.frame(NMDS_result["points"])
  species_coords <- data.frame(NMDS_result["species"])
  
  replicate_coordinates$replicate <- rownames(replicate_coordinates)
  replicate_coordinates <- merge(replicate_coordinates, data, by =merge_by, all.x = TRUE)
  #criteria color
  colnames(replicate_coordinates)[colnames(replicate_coordinates) == column_c] <- 'crit_color'
  colnames(replicate_coordinates)[colnames(replicate_coordinates) == column_s] <- 'crit_shape'
  replicate_coordinates$crit_color <- as.factor(replicate_coordinates$crit_color )
  
  # Now we can build the plot! 
  
  p <- ggplot() +
    geom_point(data = replicate_coordinates, aes(x = points.MDS1, y = points.MDS2, 
                                                 color = crit_color, shape =crit_shape), size=2) + 
    annotate(geom = "label", 
             x = min(replicate_coordinates$points.MDS1)+1, 
             y =min(replicate_coordinates$points.MDS2), 
             size =4,
             label = paste("Stress: ", round(NMDS_result$stress, digits = 3)))+ #+
    labs(color = column_c) + theme( axis.text = element_text(size = 12)) + theme_light()
  #xlim(min(replicate_coordinates$points.MDS1)*0.9, max(replicate_coordinates$points.MDS1)*1.1) +
  #ylim(min(replicate_coordinates$points.MDS2)*0.9, max(replicate_coordinates$points.MDS2)*1.1) 
  
  #print(p)
  list_data <- list(replicate_coordinates, p)
  return(list_data)
}