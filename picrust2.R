#####################################################################
################### Picrust2 results processing #####################
#####################################################################
 
############### Install and load libraries 
source("scripts/01_install.R", local = TRUE)

############### Load results
kegg_abund <- ko2kegg_abundance(file = "results/picrust2/results/KO_metagenome_out/pred_metagenome_unstrat.tsv")
#head(kegg_abund)

############### Counting Normalization
cpm_norm <- cpm(kegg_abund) ## cpm Normalization insight samples
thresh <- cpm_norm > 1000 ## Threshold 1000 reads minimum 
table(rowSums(thresh))
keep <- rowSums(thresh) >= 18 ## almost in 18 samples appear counts
kegg_tmp <- kegg_abund[names(which(keep)),]

############### Load metadata
metadata_f <- read_delim("results/picrust2/metadata.tsv",delim = "\t", escape_double = FALSE,trim_ws = TRUE)
as.tibble(metadata_f)
group <- "Grupo"

############### Kegg anotation and barplot
source("scripts/06_annotation_barplot.R", local = TRUE)
listKo <- rownames(kegg_tmp)
barplot_annotation(listKo,kegg_tmp,color_palette)

############### Differential Analysis using ALDEx2 method for Host Macroalgae
source("scripts/07_kegg_abundance_heatmap.R", local = TRUE)
daa_results_df_grupo <-
  pathway_daa(abundance = kegg_tmp,
    metadata = metadata_f,
    group = "Grupo",
    daa_method = "ALDEx2",
    select = NULL,
    reference = NULL
  )
daa_results_df_grupo_select  <-  daa_results_df_grupo[daa_results_df_grupo$method == "ALDEx2_Kruskal-Wallace test", ]
daa_results_df_grupo_heatmap <-  pathway_annotation(pathway = "KO",daa_results_df = daa_results_df_grupo_select,ko_to_kegg = TRUE)
kegg_abundance_heatmap( abundance = kegg_abund,
                    daa_results_df= daa_results_df_grupo_heatmap,
                    Group = metadata_f$Grupo,
                    p_values_threshold = 0.05,
                    order = "p_values",
                    select = NULL,
                    ko_to_kegg = FALSE,
                    p_value_bar = TRUE,
                    x_lab = "pathway_name")

############### Differential Analysis using ALDEx2 method for Gender of Macroalgae
daa_results_df_gender <-
  pathway_daa( abundance = kegg_tmp,
    metadata = metadata_f,
    group = "Genero_update",
    daa_method = "ALDEx2",
    select = NULL,
    reference = NULL
  )

daa_results_df_gender_select  <- daa_results_df_gender[daa_results_df_gender$method == "ALDEx2_Kruskal-Wallace test", ] 
daa_results_df_gender_heatmap <-  pathway_annotation(pathway = "KO", daa_results_df = daa_results_df_gender_select[c(1:45),],ko_to_kegg = TRUE)

kegg_abundance_heatmap( abundance = kegg_abund,
                    daa_results_df= daa_results_df_gender_heatmap,
                    Group = metadata_f$Genero_update,
                    p_values_threshold = 0.05,
                    order = "p_values",
                    select = NULL,
                    ko_to_kegg = FALSE,
                    p_value_bar = TRUE,
                    x_lab = "pathway_name")