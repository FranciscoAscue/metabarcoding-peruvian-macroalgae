library(ggvenn)
library(phyloseq)

filter_uniq_ASV <- function(otu_table_t, sample, min_count = 0){
  index_pos <- match(sample,colnames(otu_table_t))
  tmp <- otu_table_t[otu_table_t[,index_pos]> min_count ,index_pos]
  return(tmp)
}

inter_specie_asv <- function(otu_table, vector){
  tmp <- as.data.frame(ps_rar@otu_table)
  data_venn<- tmp[,vector]
  data_venn <- data_venn[rowSums(data_venn[])> 0,]
  for( i in 1:length(vector)){
    indx <- which(data_venn[,i] > 0)
    ind0 <- which(data_venn[,i] == 0)
    data_venn[indx,i] <- rownames(data_venn)[indx]
    data_venn[ind0,i] <- rep(NA,length(data_venn[ind0,i]))
  }
  plot <- ggVennDiagram(data_venn, label_alpha =0.4, label = "both")
  return(list(DATA = data_venn, PLOT = plot))
}

  

sum_asv_gender <- function(otu_table, vector){
  tmp <- as.data.frame(ps_rar@otu_table)
  data<- tmp[,vector]
  data <- data[rowSums(data[])> 0,]
  return(rowSums(data))
}

Rhodophyta <- sum_asv_gender(ps_rar@otu_table, c(vector_chondracantus_pucusanad, vector_chondracantus_lapunta))
Phaeophyceae <- sum_asv_gender(ps_rar@otu_table, vector_macrocystis)
Chlorophyta <- sum_asv_gender(ps_rar@otu_table, c(vector_ulva_pucusana, vector_ulva_lapunta,"PBMD006D"))
names(Rhodophyta)

sum(Rhodophyta,Phaeophyceae,Chlorophyta)

lista <- list(Rhodophyta = names(Rhodophyta) , Phaeophyceae = names(Phaeophyceae),
              Chlorophyta = names(Chlorophyta))

ggvenn(
  lista,
  fill_color = c("#DE3163", "#DFFF00", "#2ECC71"),
  stroke_size = 0.5, set_name_size = 4
)

#### Intesection of Rhodophyta, Phaeophyceae and Chlorophyta
uniq1 <- intersect(names(Chlorophyta),names(Phaeophyceae))
uniq2 <- intersect(names(Rhodophyta),names(Phaeophyceae))
uniq3 <- intersect(names(Rhodophyta),names(Chlorophyta))
uniq <- intersect(uniq2,uniq3)

tax_intersec <- ps_rar@tax_table[uniq1,]
table(as.data.frame(tax_intersec)$Phylum)


tax_intersec <- ps_rar@tax_table[uniq2,]
table(as.data.frame(tax_intersec)$Phylum)

tax_intersec <- ps_rar@tax_table[uniq3,]
table(as.data.frame(tax_intersec)$Phylum)


tax_intersec <- ps_rar@tax_table[uniq,]
table(as.data.frame(tax_intersec)$Genus)

write.csv(as.data.frame(tax_intersec), "190_intersec.csv")
getwd()

otu_ps <- t(ps@otu_table)
filter_uniq_ASV(otu_ps,sample = "16S-01-ARM-A")