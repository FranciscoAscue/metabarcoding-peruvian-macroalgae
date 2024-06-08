list_fastq <- function(ruta = getwd(), pattern = c("_1.fastq", "_2.fastq")){
  #Set working directory
  setwd(ruta)
  fastq= paste0(ruta,"/data/rawdata")
  list.files(fastq)
  lecturasf <- sort(list.files(fastq, pattern = pattern[1], full.names = T))
  lecturasr <- sort(list.files(fastq, pattern = pattern[2], full.names = T))
  nombres_muestras <- sapply(strsplit(basename(lecturasf),"_"), `[`,1)
  return(list(lf = lecturasf, lr = lecturasr, name = nombres_muestras))
  
}
