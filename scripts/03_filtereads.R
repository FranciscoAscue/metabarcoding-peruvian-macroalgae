
#----filtrado y recorte de lecturas----

filter_reads <- function(name, lf, lr, trunc){
  
  filtF <- file.path(paste0(getwd(),"/data/processed_data"), "filtered_F", 
                     paste0(name, "_filt_1.fastq.gz"))
  filtR <- file.path(paste0(getwd(),"/data/processed_data"), "filtered_R", 
                     paste0(name, "_filt_2.fastq.gz"))
  names(filtF) <- name
  names(filtR) <- name

  out <- filterAndTrim(lf,filtF,lr,filtR, 
                       truncLen = c(trunc,270), 
                       maxN = 0,
                       maxEE = c(2,2),
                       truncQ = 2,
                       rm.phix = T,
                       compress = T,
                       minQ = 2,
                       multithread = T)
  
  recovery <- round(out[2]/out[1]*100,2)
  log <- paste0("Secuencias recuperadas : ",recovery, "%")
  cat(log)
  return(list(LOG=log,OUT=out))
}


