helper = function(folder, algo){
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install()
  }
  if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
    BiocManager::install(c("VariantAnnotation"))
  }
  
  if (!requireNamespace("RWeka", quietly = TRUE)) {
    install.packages('RWeka')
  }
  library(tidyr)
  library(RWeka)
  library(VariantAnnotation)
  library(farff)
  
  
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  mutect2.features = c('NLOD', 'MQ', 'DP', 'MQRankSum', 'HCNT' ,'ClippingRankSum', 'ECNT', 'FS', 'GC', 'HRun','HaplotypeScore','MAX_ED','MIN_ED','MQ0','PON','QD','ReadPosRankSum','TLOD')
  freebayes.features = c('AB','ABP','AC','AF','AN','AO','CIGAR','DECOMPOSED','DP','DPB','DPRA','END','EPP','EPPR','GTI','LEN','MEANALT','MIN_DP','MQM','MQMR','NS','NUMALT','ODDS','OLD_VARIANT','PAIRED','PAIREDR','PAO','PQA','PQR','PRO','QA','QR','RO','RPL','RPP','RPPR','RPR','RUN','SAF','SAP','SAR','SOMATIC','SRF','SRR','TYPE','technology.illumina')
  varscan.features = c('DP','GPV','SOMATIC','SPV','SS','SSC')
  vardict.features = c('AF','DP','END','LSEQ','MSI','MSILEN','RSEQ','SAMPLE','SHIFT3','SOMATIC','SOR','SSF','STATUS','TYPE','VD')
  
  if (algo == 'mutect2'){
    all.features = c(mutect2.features)
    vcfFilePath <- Sys.glob(file.path(getwd(), folder, "*mutect2.vcf.gz"))
    x<-list(vcfFilePath)
    vcf.tbi <- Sys.glob(file.path(getwd(), folder, "*mutect2*gz.tbi"))
    tbi<-list(vcf.tbi)
  }
  
  if (algo == 'freebayes'){
    all.features = c(freebayes.features)
    vcfFilePath <- Sys.glob(file.path(getwd(), folder, "*freebayes.vcf.gz"))
    x<-list(vcfFilePath)
    tbi <- Sys.glob(file.path(getwd(), folder, "*freebayes*gz.tbi"))
    tbi<-list(tbi)
  }
  
  if (algo == 'varscan'){
    all.features = c(varscan.features)
    vcfFilePath <- Sys.glob(file.path(getwd(), folder, "*varscan.vcf.gz"))
    x<-list(vcfFilePath)
    tbi <- Sys.glob(file.path(getwd(), folder,"*varscan*gz.tbi"))
    tbi<-list(tbi)
  }
  
  if (algo == 'vardict'){
    all.features = c(vardict.features)
    vcfFilePath <- Sys.glob(file.path(getwd(), folder, "*vardict.vcf.gz"))
    x<-list(vcfFilePath)
    tbi <- Sys.glob(file.path(getwd(), folder,"*vardict*gz.tbi"))
    tbi<-list(tbi)
  }
  

  print("Reading Bed File")
  #reading bed file
  bed = data.frame(read.table(Sys.glob(file.path(getwd(), '/', folder, "*truth.bed")) ,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
  #adding truth label to all the bed file elements
  bed$TRUTH = TRUE
  names(bed) <- c("Chr","START_POS_REF","END_POS_REF", "TRUTH")
  
  
  svp<-ScanVcfParam(info=all.features, samples=suppressWarnings(scanVcfHeader(x[[1]])@samples))

  
  print("Reading VCF Files")
  vcf<- suppressWarnings(readVcf(tbi[[1]], genome=seqinfo(scanVcfHeader(x[[1]])), svp))

  
  # remove duplicates from vcfs
  print("Removing duplicates")
  vcf=vcf[!duplicated(rowRanges(vcf))]
  
  H=header(vcf)

  # sample name
  sampleid.t = H@header@listData$PEDIGREE$Derived
  sampleid.n = H@header@listData$PEDIGREE$Original
  
  # get rows that are SNVs
  snv<-isSNV(vcf, singleAltOnly=FALSE)

  # extract passed calls from each caller
  if (algo == 'mutect2'){
    pass<- vcf@fixed$FILTER=="PASS" | vcf@fixed$FILTER=="MinAF"
  } else {
    pass<- vcf@fixed$FILTER=="PASS" 
  }
  

  # get passed snv calls from all callers
  snv_pass=c(rowRanges(vcf[pass&snv]),
             ignore.mcols=T)
  snv_pass= unique(snv_pass)
  
  snv.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf), type="equal"))
  
  # tabulate caller features
  print('Tabulating caller features')
  if (length(snv.index)!=0) {
    na=F
    
    # convert vcf to VRanges then to Granges, keep metadata columns
    vr<- as(vcf[snv.index], "VRanges")
    vr=GenomicRanges::split(vr, vr@sampleNames)
    gr=GRanges(vr[[sampleid.t]])
    mcols(gr)=cbind(mcols(gr), data.frame(REF=ref(vr[[sampleid.t]]), ALT=alt(vr[[sampleid.t]]), T_totalDepth=totalDepth(vr[[sampleid.t]]), 
                                                T_refDepth=refDepth(vr[[sampleid.t]]), T_altDepth=altDepth(vr[[sampleid.t]]),N_totalDepth=totalDepth(vr[[sampleid.n]]), 
                                                N_refDepth=refDepth(vr[[sampleid.n]]), N_altDepth=altDepth(vr[[sampleid.n]]), stringsAsFactors=F))
    
    gr <- unique(gr)
    gr$FILTER=vcf[snv.index]@fixed$FILTER
    
  } else {
    print("Warning: There are no passed SNVs.")
    na=T
  }
  
  ## merge 4 vcfs and meta data
  print("Merging data")
  cols = c(names(vcf@info@listData),'REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
  
  names= paste(cols, sep="")
  
  meta_data=data.frame(snv_pass)[,1:3]
  meta_data[c(names)]=NA
  
  if (na==F) {
    meta_data[Biostrings::match(gr, snv_pass), names]=data.frame(mcols(gr)[,cols])
  }
  
  # extract reference allele
  ref=meta_data["REF"] 
  ref.ind=which(!is.na(ref), arr.ind=T) # find non-NA alleles
  ref.ind=ref.ind[order(ref.ind[,1]),] #order array indices by row number
  ref.ind=ref.ind[!duplicated(ref.ind[,1]),] #get first non-NA value for each row
  meta_data$REF=ref[ref.ind]
  
  # extract alternate allele
  alt=meta_data[c("ALT")]
  suppressWarnings(alt$ALT[(!is.na(alt$ALT) & nchar(alt$ALT)!=1)] <- substrRight(alt$ALT[(!is.na(alt$ALT) & nchar(alt$ALT)!=1)],1))
  alt.ind=which(!is.na(alt), arr.ind=T) #find all non-NA alleles
  alt.ind=alt.ind[order(alt.ind[,1]),]  #order array indices by row number
  alt.ind=alt.ind[!duplicated(alt.ind[,1]),]# take the first non-NA allele of each row
  meta_data$ALT=alt[alt.ind]  
  
  # check for list objects in columns
  for (i in 1:ncol(meta_data)) {
    if(class(meta_data[,i])=='list'){
      meta_data[,i] = unlist(meta_data[,i])
    }
  }
  
  # make filters logical
  if (algo == "mutect2"){
    meta_data$FILTER[meta_data$FILTER != "PASS" & meta_data$FILTER != "MinAF"] <- FALSE
  } else{
    meta_data$FILTER[meta_data$FILTER != "PASS"] <- FALSE
  }
  
  meta_data$FILTER[is.na(meta_data$FILTER)] <- FALSE
  meta_data$FILTER[meta_data$FILTER == "PASS"] <- TRUE
  meta_data$FILTER[meta_data$FILTER == "MinAF"] <- TRUE
  meta_data$FILTER <- as.logical(meta_data$FILTER)
  
  # get all ref/alternates
  meta_data$REF_MFVdVs<-paste(meta_data$REF, sep ="/")
  meta_data$ALT_MFVdVs<-paste(meta_data$ALT, sep ="/")
  
  # sample name
  meta_data$Sample_Name <- sampleid.t
  
  #chromosome tag
  meta_data$seqnames = gsub('chr','',meta_data$seqnames)
  
  meta_data[,all.features] <- lapply (meta_data[,all.features], as.numeric)
  
  # keep important columns and rename columns
  feature.cols = c(all.features)
  
  parse_snv <- meta_data[,c("seqnames","start","end","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                            "FILTER",
                            all.features)]
  colnames(parse_snv) <- c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                           "FILTER2",
                           feature.cols)
  
  #merge table with truth
  print("Merging table with truth")
  parse_snv_with_truth = merge(x=parse_snv,y=bed,all.x=TRUE)
  
  #non bed data replaced with FALSE
  parse_snv_with_truth$TRUTH[is.na(parse_snv_with_truth$TRUTH)]= FALSE
  
  print("Sorting data")
  # sort table
  parse_snv_with_truth = parse_snv_with_truth[order(parse_snv_with_truth$START_POS_REF),]
  parse_snv_with_truth = parse_snv_with_truth[order(parse_snv_with_truth$Chr),]
  return(parse_snv_with_truth)
}

writeArff = function (algo){
  folders = c("real1", "syn1", "syn2", "syn3", "syn4", "syn5")
  
  df = data.frame()
  
  for (folder in folders){
    temp.df = helper(folder, algo)
    df = rbind(df, temp.df)
  }
  write.arff(parse_snv_with_truth, file = paste0(getwd(), "/", algo, ".arff"), eol = "\n")
  return(df)
}

setwd("C:/Users/mrlim/OneDrive/Desktop/NUS/Y2S2/CS4220/Data2")
parse_snv_with_truth = writeArff("vardict")
