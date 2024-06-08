
pcoa_plot <- function(pcoa_result, metadata_matrix, my_pal, color_label = "Host Macroalgae", 
  shape_label = "Location", criteria = "Genero" ){

    pcoa_df <- data.frame(Sample = rownames(pcoa_result$points),
                          PCoA1 = pcoa_result$points[, 1],
                          PCoA2 = pcoa_result$points[, 2],
                          Genero = metadata_matrix$Genero_update ,
                          Grupo = metadata_matrix$Grupo,
                          Procedencia= metadata_matrix$Procedencia)
    if( criteria = "Genero"){
          plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Genero, shape = Procedencia)) 
    }else{
          plot <-  ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Grupo, shape = Procedencia)) 
    }

      plot + geom_point(size = 4) + theme_base() + labs(
        x = paste("PCoA1 (", round(pcoa_result$eig[1] / sum(pcoa_result$eig) * 100, 2), "%)", sep = ""),
        y = paste("PCoA2 (", round(pcoa_result$eig[2] / sum(pcoa_result$eig) * 100, 2), "%)", sep = ""),
        color = color_label, shape = shape_label) 
       + scale_color_manual(values = c(my_pal))
       + scale_fill_manual(values = c(paste(my_pal, sep = ""))) -> plot_save

       return(plot_save)
    }
}