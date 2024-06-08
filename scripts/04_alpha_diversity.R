library(vegan)
library(phyloseq)
library(ggplot2)
library(ggsignif)



efective_numbers <- function(abundance_matrix, list, target = "sample",
                             diversity = "Simpson"){
  
  if(target == "total"){
    tmp <- t(abundance_matrix[,rownames(list2)])
    total_abundance <- rowSums(tmp)
    
    simpson_index <- diversity(total_abundance, index = "simpson")
    simpson_effective_numbers <- 1 / (1 - simpson_index)
    
    shannon_index <- diversity(total_abundance, index = "shannon")
    shannon_effective_numbers <- exp(shannon_index)
    
    
    text <- paste0("shannon_effective_numbers : ", shannon_effective_numbers,
        " simpson_effective_numbers : ", simpson_effective_numbers)
    
    return(text)
  }
  
  if(target == "asv"){
    total_abundance <- abundance_matrix[,rownames(list)]
    #total_abundance <- rowSums(tmp)
  }
  
  if(target == "sample"){
    total_abundance <- t(abundance_matrix[,rownames(list)])
    #total_abundance <- rowSums(tmp)
  }
 
  if(diversity == "Simpson"){
    simpson_index <- diversity(total_abundance, index = "simpson")
    simpson_effective_numbers <- 1 / (1 - simpson_index)
    #cat("Simpson Effective numbers: ", simpson_effective_numbers)
    return(simpson_effective_numbers)
  }
  
  if(diversity == "Shannon"){
    shannon_index <- diversity(total_abundance, index = "shannon")
    shannon_effective_numbers <- exp(shannon_index)
   # cat("Shannon Effective numbers: ",shannon_effective_numbers)
    return(shannon_effective_numbers)
  }
}

barplot_diversity <- function(efn1, efn2, pal_color, y_lab = "Simpson effective numbers"){
  datos <- data.frame(
    Samples = factor(c(names(efn1), names(efn2))),
    Effective_Numbers = c(as.vector(efn1), as.vector(efn2)),
    Location = factor(c(rep("La Punta", length(efn1)), rep("Pucusana", length(efn2))))
  )
  print(pairwise.t.test(datos$Effective_Numbers, datos$Location, p.adjust.method = "bonferroni"))
  
  ggplot(datos, aes(x = Location, y = Effective_Numbers, fill = Location)) +
    geom_boxplot() +
    geom_jitter(shape = 18, position = position_jitter(0.2), size = 2, alpha = 0.7) +  # Añadir puntos individuales
    geom_signif(comparisons = list(c("La Punta", "Pucusana")), 
                map_signif_level = TRUE, 
                test = "t.test") +
    scale_fill_manual(values = Location) +  
    labs(title = "" ,
         x = "",
         y = y_lab) +
    theme_minimal() 
}