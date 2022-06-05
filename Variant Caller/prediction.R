getSubset = function(folder){
  if (!requireNamespace("RWeka", quietly = TRUE)) {
    install.packages('RWeka')
  }
  library(RWeka)
  
  #retrieve ARFF file 
  arff=read.arff(Sys.glob(file.path(getwd(), "*prediction.arff")))
  
  #remove empty TRUTH table
  arff = subset(arff, select = -c(TRUTH))
  
  #keep only records that are deemed to be true 
  arff2 = subset(arff, `predicted TRUTH` == TRUE)
  
  return(arff2)
}


writeToBed = function(file){
  bed = subset(file, select = c('Chr', 'START_POS_REF', 'END_POS_REF'))
  bed = bed[order(bed$START_POS_REF),]
  bed = bed[order(bed$Chr),]
  write.table(bed, paste0(getwd(), "/", "newprediction.bed"), row.names = F, quote = F, sep = '\t')
}

setwd("C:/Users/mrlim/OneDrive/Desktop/NUS/Y2S2/CS4220/Data2/real2_part2")
arff2 = getSubset()
writeToBed(arff2)

