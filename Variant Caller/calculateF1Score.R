# Calculate F1 stats

get.F1.values = function() {
  bed = data.frame(read.table(Sys.glob(file.path(getwd(), "*truth*.bed")) ,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
  bed$TRUTH = T
  prediction = data.frame(read.table(Sys.glob(file.path(getwd(), "prediction.bed")) ,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
  prediction$predict = T
  final = merge(x=bed,y=prediction,all.x=TRUE, all.y=TRUE)
  final$TRUTH[is.na(final$TRUTH)]= FALSE
  final$predict[is.na(final$predict)]= FALSE
  df = final
  
  table = data.frame(matrix(nrow = 1,ncol = 7))
  colnames(table) = c('Sample', 'TP','FP','FN','Precision','Recall','F1')
  
  table$Sample = 'test'
  
  table$TP = nrow(df[df$predict==TRUE & df$TRUTH==TRUE,])
  table$FP = nrow(df[df$predict==TRUE & df$TRUTH==FALSE,])
  table$FN = nrow(df[df$predict==FALSE & df$TRUTH==TRUE,])
  
  table$Precision = table$TP/(table$TP + table$FP)
  table$Recall = table$TP/(table$TP + table$FN)
  table$F1 = (2*table$Precision*table$Recall)/(table$Precision + table$Recall)
  
  return(table)
  
}

setwd("C:/Users/mrlim/OneDrive/Desktop/NUS/Y2S2/CS4220/Data2/real2_part1")
table = get.F1.values()
table
