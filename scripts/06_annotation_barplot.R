calculate_stats <- function(row) {
  mean_val <- mean(row)
  sd_val <- sd(row)
  c(mean = mean_val, sd = sd_val)
}

color_palette <- c("#050f2c","#003666","#00aeff","#3369e7","#8e43e7","#b84592","#ff4f81",
                  "#ff6c5f","#ffc168","#2dde98","#1cc7d0","#d20962","#f47721","#7ac143",
                  "#00a78e","#ce181e","#fdb813","#88aca1","#788cb6","#057855","#cd595a",
                  "#f8dfc2","#0af167")

barplot_annotation <- function(listKo,kegg_tmp,color_palette){
    listEntry <- data.frame(KO = listKo, NAME = rep(NA,length(listKo)), 
                        CLASS = rep(NA,length(listKo)))  
    for (i in 1:length(listKo)) {
      tryCatch({ kegg_annot <- keggGet(listEntry$KO[i])[[1]][c("NAME","CLASS")]
        listEntry$NAME[i] <- kegg_annot$NAME
        listEntry$CLASS[i] <- kegg_annot$CLASS
      }, error = function(e) {cat("Error en la iteración ", listEntry$KO[i],"\n")})}
    errorbar_abundance_mat <- as.matrix(kegg_tmp)
    relative_abundance_mat <- apply(t(errorbar_abundance_mat),1, function(x) x/sum(x))
    stats_mat <- t(apply(relative_abundance_mat, 1, calculate_stats))
    final_mat <- cbind(relative_abundance_mat, stats_mat) 
    final_df <- as.data.frame(final_mat)
    final_df$KO <- rownames(final_df)
    via1_2 <- listEntry
    via1_2$VIA1 <- rep(NA,251)
    via1_2$VIA2 <- rep(NA,251)
    annot_string <- strsplit(via1_2$CLASS, split = ";")
    for( i in 1:251){
      via1_2$VIA1[i] <- annot_string[[i]][1]
      via1_2$VIA2[i] <- annot_string[[i]][2]
    }
    final_df_anotado <- merge(final_df, via1_2, by ="KO")
    final_df_anotado_da <- drop_na(final_df_anotado)
    final_df_anotado_filter <- final_df_anotado_da[final_df_anotado_da$VIA1 != "Human Diseases",]

    ############### Barplot of kegg anotation
    final_df_anotado_filter %>% group_by(VIA1) %>% summarise(n = n()/115)

    tmp_Cp <- final_df_anotado_filter %>% filter(VIA1 == "Cellular Processes (CP)")
    tmp_M <- final_df_anotado_filter %>% filter(VIA1 == "Metabolism (M)")
    tmp_eip <- final_df_anotado_filter %>% filter(VIA1 == "Environmental Information Processing (EIP)")
    tmp_gip <- final_df_anotado_filter %>% filter(VIA1 == "Genetic Information Processing (GIP)")
    tmp_os <- final_df_anotado_filter %>% filter(VIA1 == "Organismal Systems (OS)")

    ggplot(data = final_df_anotado_filter, aes(x = VIA1, y = mean + sd, fill = VIA2)) +
      geom_bar(stat = "identity", position = "fill") +
      scale_fill_manual(values = color_palette) +
      theme_minimal() +
      labs(
        x = " ",
        y = "Relative abundance",
        fill = "Predicted pathways" 
      )  + theme(legend.position = "bottom") -> bar_annot_plot
      return(bar_annot_plot) 
}