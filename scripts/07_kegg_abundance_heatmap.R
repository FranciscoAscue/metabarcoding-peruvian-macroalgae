kegg_abundance_heatmap <- function (abundance, daa_results_df, Group, ko_to_kegg = FALSE, 
    p_values_threshold = 0.05, order = "group", select = NULL, 
    p_value_bar = TRUE, colors = NULL, x_lab = NULL) 
{
  missing_pathways <- daa_results_df[is.na(daa_results_df$pathway_name), 
                                     "feature"]
  if (length(missing_pathways) > 0) {
    message("The following pathways are missing annotations and have been excluded: ", 
            paste(missing_pathways, collapse = ", "))
    message("You can use the 'pathway_annotation' function to add annotations for these pathways.")
  }
  daa_results_df <- daa_results_df[!is.na(daa_results_df[,x_lab]), ]
  if (is.null(x_lab)) {
    if (ko_to_kegg == TRUE) {
      x_lab <- "pathway_name"
    }
    else {
      x_lab <- "description"
    }
    if (is.null(daa_results_df$pathway_name) & is.null(daa_results_df$description)) {
      message("Please use pathway_annotation to annotate the daa_results_df")
    }
  }
  if (!(x_lab %in% colnames(daa_results_df))) {
    message("There is no x_lab you defined in daa_results_df")
  }
  if (nlevels(factor(daa_results_df$method)) != 1) {
    message("There are more than one method in daa_results_df$method, please filter it.")
  }
  if (is.null(colors)) {
    colors <- c("#d93c3e", "#3685bc", "#6faa3e", "#e8a825", 
                         "#c973e6", "#ee6b3d", "#2db0a7", "#f25292")[1:nlevels(as.factor(Group))]
  }
  errorbar_abundance_mat <- as.matrix(abundance)
  daa_results_filtered_df <- daa_results_df[daa_results_df$p_adjust < p_values_threshold, ]
  if (!is.null(select)) {
    daa_results_filtered_sub_df <- daa_results_filtered_df[daa_results_filtered_df$feature %in% select, ]
  }
  else {
    daa_results_filtered_sub_df <- daa_results_filtered_df
  }  
  if (nrow(daa_results_filtered_sub_df) == 0) {
    stop("The feature with statistically significance is zero, pathway_errorbar can't do the visualization.")
  }
  relative_abundance_mat <- apply(t(errorbar_abundance_mat),1, function(x) x/sum(x))
  sub_relative_abundance_mat <- relative_abundance_mat[rownames(relative_abundance_mat) %in% daa_results_filtered_sub_df$feature, ]
  error_bar_matrix <- cbind(sample = colnames(sub_relative_abundance_mat), 
                            group = Group, t(sub_relative_abundance_mat))
  error_bar_df <- as.data.frame(error_bar_matrix)
  error_bar_df$group <- factor(Group, levels = levels(as.factor(Group)))
  error_bar_pivot_longer_df <- tidyr::pivot_longer(error_bar_df, -c(sample, group))
  error_bar_pivot_longer_tibble <- mutate(error_bar_pivot_longer_df, group = as.factor(group))
  error_bar_pivot_longer_tibble$sample <- factor(error_bar_pivot_longer_tibble$sample)
  error_bar_pivot_longer_tibble$name <- factor(error_bar_pivot_longer_tibble$name)
  error_bar_pivot_longer_tibble$value <- as.numeric(error_bar_pivot_longer_tibble$value)
  error_bar_pivot_longer_tibble_summarised <- error_bar_pivot_longer_tibble %>% 
    group_by(name, group) %>% summarise(mean = mean(value), sd = stats::sd(value))
  error_bar_pivot_longer_tibble_summarised <- error_bar_pivot_longer_tibble_summarised %>% mutate(group2 = "nonsense")
  switch(order, 
    p_values = {order <- order(daa_results_filtered_sub_df$p_adjust)}, 
    name = {order <- order(daa_results_filtered_sub_df$feature)}, 
    group = {daa_results_filtered_sub_df$pro <- 1
    for (i in levels(error_bar_pivot_longer_tibble_summarised$name)) {
      error_bar_pivot_longer_tibble_summarised_sub <- error_bar_pivot_longer_tibble_summarised[error_bar_pivot_longer_tibble_summarised$name == i, ]
      pro_group <- error_bar_pivot_longer_tibble_summarised_sub[error_bar_pivot_longer_tibble_summarised_sub$mean == max(error_bar_pivot_longer_tibble_summarised_sub$mean),]$group
      pro_group <- as.vector(pro_group)
      daa_results_filtered_sub_df[daa_results_filtered_sub_df$feature == i, ]$pro <- pro_group
    }
    order <- order(daa_results_filtered_sub_df$pro, daa_results_filtered_sub_df$p_adjust)
  }, 
  pathway_class = { 
    if (!"pathway_class" %in% colnames(daa_results_filtered_sub_df)) {
      stop("Please use pathway_annotation function to annotate the pathway_daa results")
    }
    order <- order(daa_results_filtered_sub_df$pathway_class,daa_results_filtered_sub_df$p_adjust)
  }, { order <- order })
  daa_results_filtered_sub_df <- daa_results_filtered_sub_df[order,]
  error_bar_pivot_longer_tibble_summarised_ordered <- data.frame(name = NULL,group = NULL, mean = NULL, sd = NULL)
  for (i in daa_results_filtered_sub_df$feature) {
    error_bar_pivot_longer_tibble_summarised_ordered <- rbind(error_bar_pivot_longer_tibble_summarised_ordered, 
                                                              error_bar_pivot_longer_tibble_summarised[error_bar_pivot_longer_tibble_summarised$name == i, ])
  }
  if (ko_to_kegg == FALSE) {
    error_bar_pivot_longer_tibble_summarised_ordered[, x_lab] <- rep(daa_results_filtered_sub_df[,x_lab], each = length(levels(factor(error_bar_pivot_longer_tibble_summarised_ordered$group))))
  }
  if (ko_to_kegg == TRUE) {
    error_bar_pivot_longer_tibble_summarised_ordered$pathway_class <- rep(daa_results_filtered_sub_df$pathway_class, 
                                                                          each = length(levels(factor(error_bar_pivot_longer_tibble_summarised_ordered$group))))
  }
  error_bar_pivot_longer_tibble_summarised_ordered$name <- factor(error_bar_pivot_longer_tibble_summarised_ordered$name, 
                                                                  levels = rev(daa_results_filtered_sub_df$feature))
  
  daa_results_filtered_sub_df <- cbind(daa_results_filtered_sub_df, 
                                       negative_log10_p = -log10(daa_results_filtered_sub_df$p_adjust), 
                                       group_nonsense = "nonsense", log_2_fold_change = NA)
  
  for (i in daa_results_filtered_sub_df$feature) {
    mean <- error_bar_pivot_longer_tibble_summarised_ordered[error_bar_pivot_longer_tibble_summarised_ordered$name %in% i, ]$mean
    daa_results_filtered_sub_df[daa_results_filtered_sub_df$feature == i, ]$log_2_fold_change <- log2(mean[1]/mean[2])
  }
  daa_results_filtered_sub_df$feature <- factor(daa_results_filtered_sub_df$feature, 
                                                levels = rev(daa_results_filtered_sub_df$feature))
  error_bar_pivot_longer_tibble_summarised_ordered$Mean_SD <- as.numeric(error_bar_pivot_longer_tibble_summarised_ordered$mean) + 
    as.numeric(error_bar_pivot_longer_tibble_summarised_ordered$sd)
  datos_reshaped <- dcast(error_bar_pivot_longer_tibble_summarised_ordered, pathway_name ~ group, value.var = "Mean_SD", fun.aggregate = sum)
  rownames(datos_reshaped) <- datos_reshaped$pathway_name
  group_nonsense = "nonsense"
  heatm <- merge(datos_reshaped, daa_results_filtered_sub_df[,c("pathway_name","log_2_fold_change")], by= "pathway_name")
  rownames(heatm) <- heatm$pathway_name
  
  columnas_numericas <- sapply(heatm, is.numeric)
  tmp_numericas <- heatm[, columnas_numericas]
  tmp <- daa_results_filtered_sub_df
  tmp$p_adjust <- as.character(tmp$p_adjust)
  tmp$unique <- nrow(tmp) - seq_len(nrow(tmp)) + 1
  tmp$p_adjust <- substr(tmp$p_adjust, 1, 5)
  datos_largos <- pivot_longer(datos_reshaped, 
                               cols = -pathway_name, 
                               names_to = "group", 
                               values_to = "value")
    ggplot(datos_largos, aes(x = group, y = pathway_name, fill = value)) + geom_tile() +
    scale_fill_gradient2(low = "#2471A3", mid = "#FFFFFF", high = "#E74C3C", midpoint = median(datos_largos$value)) +
    theme_minimal() + labs(fill = "Relative\nAbundance", x = "Genus of Macroalgae", y = "Pathway Name") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    daa_results_filtered_sub_df %>% ggplot2::ggplot(ggplot2::aes(feature, 
    log_2_fold_change, fill = group_nonsense)) + 
    ggplot2::geom_bar(stat = "identity",position = ggplot2::position_dodge(width = 0.8), width = 0.8) + 
    ggplot2::labs(y = "log2 fold change", x = NULL) + GGally::geom_stripped_cols() + 
    ggplot2::scale_fill_manual(values = "#76D7C4") + ggplot2::scale_color_manual(values = "#76D7C4") + 
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = "dashed",color = "black") + ggprism::theme_prism() + 
    ggplot2::scale_y_continuous(expand = c(0,0), guide = "prism_offset_minor") + 
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(), 
                   axis.line.y = ggplot2::element_blank(), axis.line.x = ggplot2::element_line(size = 0.5), 
                   axis.ticks.x = ggplot2::element_line(size = 0.5), panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(), axis.text = 
                     ggplot2::element_text(size = 10,color = "black"), axis.text.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 10, color = "black",margin = ggplot2::margin(b = 6)), 
                   axis.title.x = ggplot2::element_text(size = 11, color = "black", hjust = 0.5), legend.position = "non") + 
    ggplot2::coord_flip() +
    tmp %>% ggplot2::ggplot(ggplot2::aes(group_nonsense, p_adjust)) + ggplot2::geom_text(ggplot2::aes(group_nonsense, 
                                                                                                      unique, label = p_adjust), size = 3.5, color = "black", 
                                                                                         fontface = "bold", family = "sans") + ggplot2::labs(y = "p-value (adjusted)") + 
    ggplot2::scale_y_discrete(position = "right") + ggprism::theme_prism() + 
    ggplot2::theme(axis.ticks = ggplot2::element_blank(), 
                   axis.line = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank(), 
                   panel.grid.major.x = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                   axis.text = ggplot2::element_blank(), plot.margin = ggplot2::unit(c(0,0.2, 0, 0), "cm"), axis.title.y = 
                   ggplot2::element_text(size = 11, color = "black", vjust = 0), axis.title.x = ggplot2::element_blank(), 
                   legend.position = "non")
}