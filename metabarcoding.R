#########################################################################
#### Bacteria associated with macroalgae from the Peruvian coast:    ####
#### A functional prediction study based on 16S rRNA metagenome data ####
#########################################################################

############### Install and load libraries 
source("scripts/01_install.R", local = TRUE)

############### Quality check 
source("scripts/02_listfastq.R", local = TRUE)
lecturas <- list_fastq(pattern = c("1.fastq.gz","2.fastq.gz"))
plotQualityProfile(c(lecturas$lf[3],lecturas$lr[3]))

############### FFilter low-quality reads
source("scripts/03_filtereads.R", local = TRUE)
log_filter <- filter_reads(name = lecturas$name, lf = lecturas$lf, 
             lr = lecturas$lr, trunc = 300)

##############################################################
################## "Metabarcoding processing #################
##############################################################

############### Error modeling 
filtF <- file.path("data/processed_data/filtered_F", paste0(lecturas$name, "_filt_1.fastq.gz"))
filtR <- file.path("data/processed_data/filtered_R", paste0(lecturas$name, "_filt_2.fastq.gz"))
errR <-learnErrors(filtR, multithread = T)
errF <-learnErrors(filtF, multithread = T)
#plotErrors(errF,nominalQ = T)
#plotErrors(errR,nominalQ = T)
save(errF, errR, file = "results/error_model.RDATA")

############### ASV inference 
dadaF <- dada(filtF, err = errF, multithread = T)
dadaR <- dada(filtR, err = errR, multithread = T)
save(dadaR, dadaF, file = "results/dadasFR_270_FV.RDATA")

############### Denoising data 
pareadas <- mergePairs(dadaF, filtF, dadaR, filtR, verbose = T)
#head(pareadas[[1]])
seqTab <- makeSequenceTable(pareadas)
table(nchar(getSequences(seqTab)))
seqtab_nochim <- removeBimeraDenovo(seqTab, method = "consensus"
                                    , multithread = T,
                                    verbose = T)
#head(seqtab_nochim)
#sum(seqtab_nochim/sum(seqTab))

############### Taxonomy assignment  
ruta_clasificador <- "data/reference/silva_nr99_v138_train_set.fa.gz"
taxa <- assignTaxonomy(seqtab_nochim, ruta_clasificador , multithread = TRUE)
taxa_print <- taxa
rownames(taxa_print) <- NULL
  
############### Metadata preparation
met <- read.table("data/metadata.tsv", sep = "\t", header = TRUE, row.names = "Samples")
met <- met[lecturas$name,]
rownames(met) <- paste0(lecturas$name, "_filt_1.fastq.gz")

############### Create phyloseq object 
#rownames(otu_table(seqtab_nochim, taxa_are_rows = F))
phyloseq_ob <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = F),
                        sample_data(met),
                        tax_table(taxa))
 
############### Add DNA sequences
dna <- Biostrings::DNAStringSet(taxa_names(phyloseq_ob))
names(dna) <- taxa_names(phyloseq_ob)
ps <- merge_phyloseq(phyloseq_ob, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

############### Export ASV sequences in FASTA format
names(dna) <- taxa_names(ps)
writeXStringSet(dna, "results/picrust2/asv-seq.fna")

############### Remove contaminants and Export BIOM format
contamdf.freq <- isContaminant(ps, method="frequency", conc="quant_reading")
ps.noncontam <- prune_taxa(!contamdf.freq$contaminant, ps)
otu12 <- as(otu_table(ps.noncontam), "matrix")
md22 <- make_biom(otu12, matrix_element_type = "int")
write_biom(md22, "results/picrust2/otu_norarefy.biom")

############### Export ASV names
write.table(rownames(ps.noncontam@otu_table), "results/picrust2/asv.table", row.names = FALSE, col.names = FALSE)

############### Export taxonomy table
TAX <- as(tax_table(ps), "matrix")
write.csv(TAX,"results/taxa-algae16.csv")

############### Export ASV table
ASV <- as(otu_table(ps),"matrix")
ASV <- t(ASV)
colnames(ASV) <- gsub("_filt_1.fastq.gz","",colnames(ASV), fixed = TRUE)
write.csv(ASV, "results/asv-table.csv")

############### Save phyloseq object in Rdata
saveRDS(ps, "results/objeto_phyloseq_270_FV.RDS")

############### Statistics
summary(sample_sums(ps))
sort(sample_sums(ps))

############### Save rarefied phyloseq object in RData
set.seed(1)
ps_rar <- rarefy_even_depth(ps, min(sample_sums(ps)))

############### Save phyloseq object in RData
saveRDS(ps_rar,"results/phyloseq_rar_270_FV.RDS")

############### Remove chloroplast taxa and NA Genus
ps_rar <-subset_taxa(ps_rar,Order!="Chloroplast")
ps_rar <-subset_taxa(ps_rar,!is.na(Genus))

############### Export taxa counts by taxonomy category
taxonomy<-ps_rar@tax_table
count_genus <-  taxonomy %>% group_by(Genus) %>% 
summarise( n = n()) %>% mutate(Percent = 100*(n/sum(n))) %>% arrange(desc(Percent))
count_class <-  taxonomy %>% group_by(Class) %>% 
summarise( n = n()) %>% mutate(Percent = 100*(n/sum(n))) %>% arrange(desc(Percent))
count_phylum <-  taxonomy %>% group_by(Phylum) %>% 
summarise( n = n()) %>% mutate(Percent = 100*(n/sum(n))) %>% arrange(desc(Percent))
write.csv(count_genus, "results/genus_abundance_200224.csv")
write.csv(count_class, "results/class_abundance_200224.csv")
write.csv(count_phylum, "results/phylum_abundance_200224.csv")
Phylum<-taxonomy%>% count (Phylum)
write.csv(Phylum,"results/tabla_Phylum_ASV.csv")
Class<-taxonomy%>% count (Class)
write.csv(Class,"results/tabla_Class_ASV.csv")
Genus<-taxonomy%>% count (Genus)
write.csv(Genus,"results/tabla_genus_ASV.csv")
tmp <- as.data.frame(ps_rar@tax_table) %>% dplyr::filter(Family %in% c("Flavobacteriaceae","Saprospiraceae" ,
                                                                       "Rhodobacteraceae","Hyphomonadaceae"))
tmp <- tmp %>% dplyr::group_by(Family,Genus) %>% dplyr::summarise(n = n()) %>% dplyr::filter(n > 20)
write_csv(tmp, "results/tabla_family_plot_count.csv")


#############################################################
################### Metabarcoding Plots #####################
#############################################################

############### Taxonomy composition by Location
ps_rar_dedupe <- ps_dedupe(ps_rar, vars = c("Genero","Procedencia"),
  method = "readcount",   verbose = TRUE,n = 1,
  .keep_group_var = FALSE,
  .keep_readcount = FALSE,
  .message_IDs = FALSE,
  .label_only = FALSE,
  .keep_all_taxa = FALSE )
ps_rar_dedupe <- phyloseq_validate(ps_rar_dedupe)
ps_rar_dedupe <- tax_fix(ps_rar_dedupe)
ps_rar_dedupe <- ps_rar_dedupe %>% tax_fix(unknowns = c("endosymbionts"))
ps_rar_dedupe <- ps_rar_dedupe %>% tax_fix(unknowns = c("Unknown Family"))

for( taxa_level in c("Phylum", "Class", "Genus")){
  ps_rar_dedupe %>% comp_barplot(taxa_level, n_taxa = 19, merge_other = TRUE) +
  facet_wrap(vars(Procedencia), scales = "free") + coord_flip() + 
  ggtitle( "Relative abundance of ASV") + theme(axis.ticks.y = element_blank(), 
    strip.text = element_text(face = "bold"))
}

############### Alpha diversity (Shannon and Simpson)
source("scripts/04_alpha_diversity.R", local = TRUE)
threshold_abundance <- 0.01
ps_rar_filtrado <- prune_taxa(taxa_abundance >= threshold_abundance, ps_rar)
abundance_matrix <- as.matrix(otu_table(ps_rar_filtrado))
species = c("Rhodomelaceae","Gymnogongrus durvillei","Chondracantus chamissoi","Ulva sp")
list <- metadata %>% filter(Genero %in% species, Procedencia == "La Punta") %>%  select(Samples)
lpunta_shannon <- efective_numbers(abundance_matrix, list, diversity = "Shannon")
lpunta_simpson <- efective_numbers(abundance_matrix, list, diversity = "Simpson")
list2 <- metadata %>% filter(Genero %in% species, Procedencia == "Pucusana") %>% select(Samples)
pucusana_shannon <- efective_numbers(abundance_matrix, list2, diversity = "Shannon")
pucusana_simpson <- efective_numbers(abundance_matrix, list2, diversity = "Simpson")
Location = c("Pucusana" = "#5C80BC", "La Punta" = "#E8C547")
barplot_diversity(efn1 = lpunta_shannon, efn2 = pucusana_shannon, pal_color = Location,y_lab ="Shannon effective numbers" )
barplot_diversity(efn1 = lpunta_simpson, efn2 = pucusana_simpson, pal_color = Location)

############### pallete colours 
my_colour = list( "Group of Macroalgae" = c(Rhodophyta = "#DE3163", Chlorophyta = "#2ECC71",Phaeophyta ="#DFFF00"),
  "Genus of Macroalgae" = c("Asterfilopsis furcelatus" = "#5e4fa2", "Chondracantus chamissoi" = "#9e0142",
                            "Rhodomelaceae" = "#d53e4f", "Gymnogongrus durvillei" = "#f46d43" , 
                            "Porphyra sp" = "#fdae61", "Phyllymenia acletoi" = "#fee08b",
                            "Ulva sp" = "#3288bd", "Caulerpa filiformis" = "#66c2a5",
                            "Mazaella canaliculata" = "#ffffbf", "Macrocystis pyrifera" = "#e6f598",
                            "Codium fragile" = "#abdda4"),
  Location = c(Pucusana = "#5C80BC", "La Punta" = "#E8C547")
)

############### Beta diversity GUniFrac Heatmap
random_tree = rtree(ntaxa(ps_rar_filtrado), rooted=TRUE, tip.label=taxa_names(ps_rar_filtrado))
ps_rar_filtrado@phy_tree <- random_tree
result <- GUniFrac(t(ps_rar_filtrado@otu_table), ps_rar_filtrado@phy_tree)
dist_matrix <- result$unifracs[, , "d_1"]
metadata_matrix <- as.data.frame(ps_rar_filtrado@sam_data)
rownames(metadata_matrix) <- metadata_f$Samples
#columns <- c("Genus","Location","Group of Macroalgae","Genus of Macroalgae")
macroalgae_group <- metadata_matrix[,c("Genus of Macroalgae","Group of Macroalgae")]
procedencia_column = data.frame(Location = metadata_matrix$Location)
rownames(procedencia_column) <- rownames(metadata_matrix)
pheatmap(dist_matrix, cutree_rows = 6, main = "Gunifrac pairwise distances Heatmap",
         clustering_method = "average", border_color = NA, fontsize_row = 7,
         cluster_cols = F, cluster_rows = F,
         annotation_row = macroalgae_group,cellwidth = NA,fontsize_col = 7,
         annotation_col = procedencia_column, treeheight_row = 45, angle_col = 315,
         annotation_colors = my_colour,treeheight_col = 25, 
         color = viridis(60))

############### Beta diversity GUniFrac PCoA
source("scripts/05_beta_pcoa_diversity.R", local = TRUE)
pcoa_result <- cmdscale(dist_matrix, eig = TRUE)
pcoa_plot(pcoa_result, metadata_matrix, my_colour$`Genus of Macroalgae`, 
  color_label = "Host Macroalgae", shape_label = "Location")
pcoa_plot(pcoa_result, metadata_matrix, my_colour$`Group of Macroalgae`, 
  color_label = "Group of Macroalgae", shape_label = "Location", criteria = "Group")

############### DAPC - Adegenet 
my_pal <- c( "#5e4fa2",  "#9e0142", "#d53e4f",  "#f46d43" , "#fdae61","#fee08b","#3288bd","#66c2a5", "#ffffbf", "#e6f598", "#abdda4")
bray <- vegdist(otu_table(ps_rar_filtrado),"bray")
DAPC_METALAGA <- as.data.frame(as.matrix(bray))
grp <- find.clusters(DAPC_METALAGA)
DAPC <- dapc(DAPC_METALAGA, grp$grp)
my_metalga <- as.data.frame(DAPC$ind.coord)
my_metalga$Group <- DAPC$grp
my_metalga$id <- row.names(DAPC$tab)
metadata_matrix$id <- row.names(metadata_matrix)
GGPLOT_DAPC <- merge(my_metalga, metadata_matrix, by = "id")
colnames(GGPLOT_DAPC) <- c("id","LD1","LD2","Group","Genus","Location",
                           "Group of Macroalgae","Genus2")
ggplot(GGPLOT_DAPC,aes(x = LD1,y = LD2,color = Genus2,shape = Location)) + geom_point(size = 4) + 
  scale_color_manual(values = c(my_colour$`Genus of Macroalgae`)) + 
  scale_fill_manual(values = c(paste(my_pal, sep = ""))) +
  labs(color = "Genus of Macroalgae")
