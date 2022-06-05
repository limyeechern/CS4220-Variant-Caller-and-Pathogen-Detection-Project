writeArff = function(folder){
  
  
  setwd(paste0(getwd(), "/", folder))
  # Check for packages
  
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
  
  
  # Features to include
  mutect2.features = c('NLOD', 'MQ', 'DP', 'MQRankSum', 'HCNT')
  freebayes.features = c('DPRA','MQMR','MQM','ODDS','DP','DPB')
  varscan.features = c('SSC', 'SPV')
  vardict.features = c('AF','SSF', 'MSI')
  all.features = c(mutect2.features, freebayes.features, varscan.features, vardict.features)
  
  names.features = c(paste0(mutect2.features,'_Mutect2'),
                     paste0(freebayes.features,'_Freebayes'),
                     paste0(varscan.features,'_Varscan'),
                     paste0(vardict.features,'_Vardict'))
  
  print("Reading Bed File")
  #reading bed file
  bed = data.frame(read.table(Sys.glob(file.path(getwd(), "*truth.bed")) ,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
  #adding truth label to all the bed file elements
  bed$TRUTH = TRUE
  names(bed) <- c("Chr","START_POS_REF","END_POS_REF", "TRUTH")
  
  # Scan for vcf files
  mutect2 <- Sys.glob(file.path(getwd(), "*mutect2.vcf.gz"))
  freebayes <- Sys.glob(file.path(getwd(), "*freebayes.vcf.gz"))
  varscan <- Sys.glob(file.path(getwd(), "*varscan.vcf.gz"))
  vardict <- Sys.glob(file.path(getwd(), "*vardict.vcf.gz"))
  x<-list(mutect2,freebayes,varscan,vardict)
  
  mutect2.tbi <- Sys.glob(file.path(getwd(), "*mutect2*gz.tbi"))
  freebayes.tbi <- Sys.glob(file.path(getwd(), "*freebayes*.gz.tbi"))
  varscan.tbi <- Sys.glob(file.path(getwd(), "*varscan*gz.tbi"))
  vardict.tbi <- Sys.glob(file.path(getwd(), "*vardict*gz.tbi"))
  tbi<-list(mutect2.tbi,freebayes.tbi,varscan.tbi,vardict.tbi)
  
  # ScanVCF with required parameters
  svp_m<-ScanVcfParam(info=mutect2.features, samples=suppressWarnings(scanVcfHeader(x[[1]])@samples))
  svp_f<-ScanVcfParam(info=freebayes.features,samples=suppressWarnings(scanVcfHeader(x[[2]])@samples))
  svp_vs<-ScanVcfParam(info=varscan.features, samples=suppressWarnings(scanVcfHeader(x[[3]])@samples))
  svp_vd<-ScanVcfParam(info=vardict.features, samples=suppressWarnings(scanVcfHeader(x[[4]])@samples))
  
  
  # Read in VCF files and features
  print("Reading VCF Files")
  vcf_m2<- suppressWarnings(readVcf(tbi[[1]], genome=seqinfo(scanVcfHeader(x[[1]])), svp_m))
  vcf_f<- suppressWarnings(readVcf(tbi[[2]], genome=seqinfo(scanVcfHeader(x[[2]])), svp_f))
  vcf_vs<- suppressWarnings(readVcf(tbi[[3]], genome=seqinfo(scanVcfHeader(x[[3]])), svp_vs))
  vcf_vd<- suppressWarnings(readVcf(tbi[[4]], genome=seqinfo(scanVcfHeader(x[[4]])), svp_vd))
  
  # remove duplicates from vcfs
  print("Removing duplicates")
  vcf_m2=vcf_m2[!duplicated(rowRanges(vcf_m2))]
  vcf_f=vcf_f[!duplicated(rowRanges(vcf_f))]
  vcf_vs=vcf_vs[!duplicated(rowRanges(vcf_vs))]
  vcf_vd=vcf_vd[!duplicated(rowRanges(vcf_vd))]
  
  H_m2=header(vcf_m2)
  H_f=header(vcf_f)  
  H_vs=header(vcf_vs) 
  H_vd=header(vcf_vd)
  
  # sample name
  sampleid.t = H_m2@header@listData$PEDIGREE$Derived
  sampleid.n = H_m2@header@listData$PEDIGREE$Original
  
  # get rows that are SNVs
  snv_m2<-isSNV(vcf_m2, singleAltOnly=FALSE)
  snv_f<-isSNV(vcf_f, singleAltOnly=FALSE)
  snv_vs<-isSNV(vcf_vs, singleAltOnly=FALSE)
  snv_vd<-isSNV(vcf_vd, singleAltOnly=FALSE)  
  
  # extract passed calls from each caller
  pass_m2<- vcf_m2@fixed$FILTER=="PASS" | vcf_m2@fixed$FILTER=="MinAF"
  pass_f<- vcf_f@fixed$FILTER=="PASS" 
  pass_vs<- vcf_vs@fixed$FILTER=="PASS"
  pass_vd<- vcf_vd@fixed$FILTER=="PASS"
  
  # get passed snv calls from all callers
  snv_pass=c(rowRanges(vcf_m2[pass_m2&snv_m2]), 
             rowRanges(vcf_f[pass_f&snv_f]),  
             rowRanges(vcf_vs[pass_vs&snv_vs]), 
             rowRanges(vcf_vd[pass_vd&snv_vd]),
             ignore.mcols=T)
  snv_pass= unique(snv_pass)
  
  snv.m2.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_m2), type="equal"))
  snv.f.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_f), type="equal"))
  snv.vs.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_vs), type="equal"))
  snv.vd.index <- subjectHits(findOverlaps(snv_pass, rowRanges(vcf_vd), type="equal"))
  
  # tabulate caller features
  print('Tabulating caller features')
  if (length(snv.m2.index)!=0) {
    m2.na=F
    
    # convert vcf to VRanges then to Granges, keep metadata columns
    vr_m2<- as(vcf_m2[snv.m2.index], "VRanges")
    vr_m2=GenomicRanges::split(vr_m2, vr_m2@sampleNames)
    gr_m2=GRanges(vr_m2[[sampleid.t]])
    mcols(gr_m2)=cbind(mcols(gr_m2), data.frame(REF=ref(vr_m2[[sampleid.t]]), ALT=alt(vr_m2[[sampleid.t]]), T_totalDepth=totalDepth(vr_m2[[sampleid.t]]), 
                                                T_refDepth=refDepth(vr_m2[[sampleid.t]]), T_altDepth=altDepth(vr_m2[[sampleid.t]]),N_totalDepth=totalDepth(vr_m2[[sampleid.n]]), 
                                                N_refDepth=refDepth(vr_m2[[sampleid.n]]), N_altDepth=altDepth(vr_m2[[sampleid.n]]), stringsAsFactors=F))
    
    gr_m2 <- unique(gr_m2)
    gr_m2$FILTER=vcf_m2[snv.m2.index]@fixed$FILTER
    
  } else {
    print("Warning: There are no passed SNVs from MuTect2.")
    m2.na=T
  }
  
  if (length(snv.f.index)!=0) {
    f.na=F  
    vr_f<- suppressWarnings(as(vcf_f[snv.f.index], "VRanges"))
    vr_f=GenomicRanges::split(vr_f, vr_f@sampleNames)
    gr_f=GRanges(vr_f[[sampleid.t]])
    mcols(gr_f)=cbind(mcols(gr_f), data.frame(REF=ref(vr_f[[sampleid.t]]), ALT=alt(vr_f[[sampleid.t]]), T_totalDepth=totalDepth(vr_f[[sampleid.t]]), 
                                              T_refDepth=refDepth(vr_f[[sampleid.t]]), T_altDepth=altDepth(vr_f[[sampleid.t]]),N_totalDepth=totalDepth(vr_f[[sampleid.n]]), 
                                              N_refDepth=refDepth(vr_f[[sampleid.n]]), N_altDepth=altDepth(vr_f[[sampleid.n]]), stringsAsFactors=F))
    
    gr_f <- unique(gr_f)
    gr_f$FILTER=vcf_f[snv.f.index]@fixed$FILTER
    
  } else {
    print("Warning: There are no passed SNVs from FreeBayes.")
    f.na=T
  }
  
  if (length(snv.vs.index)!=0) {
    vs.na=F  
    vr_vs<- suppressWarnings(as(vcf_vs[snv.vs.index], "VRanges"))
    vr_vs=GenomicRanges::split(vr_vs, vr_vs@sampleNames)
    gr_vs=GRanges(vr_vs[[sampleid.t]])
    mcols(gr_vs)=cbind(mcols(gr_vs), data.frame(REF=ref(vr_vs[[sampleid.t]]), ALT=alt(vr_vs[[sampleid.t]]), T_totalDepth=totalDepth(vr_vs[[sampleid.t]]), 
                                                T_refDepth=refDepth(vr_vs[[sampleid.t]]), T_altDepth=altDepth(vr_vs[[sampleid.t]]),N_totalDepth=totalDepth(vr_vs[[sampleid.n]]), 
                                                N_refDepth=refDepth(vr_vs[[sampleid.n]]), N_altDepth=altDepth(vr_vs[[sampleid.n]]), stringsAsFactors=F))
    
    gr_vs <- unique(gr_vs)
    gr_vs$FILTER=vcf_vs[snv.vs.index]@fixed$FILTER
    
  } else {
    print("Warning: There are no passed SNVs from VarScan.")
    vs.na=T
  }
  
  if (length(snv.vd.index)!=0) {
    vd.na=F  
    vr_vd<- suppressWarnings(as(vcf_vd[snv.vd.index], "VRanges"))
    vr_vd=GenomicRanges::split(vr_vd, vr_vd@sampleNames)
    gr_vd=GRanges(vr_vd[[sampleid.t]])
    mcols(gr_vd)=cbind(mcols(gr_vd), data.frame(REF=ref(vr_vd[[sampleid.t]]), ALT=alt(vr_vd[[sampleid.t]]), T_totalDepth=totalDepth(vr_vd[[sampleid.t]]), 
                                                T_refDepth=refDepth(vr_vd[[sampleid.t]]), T_altDepth=altDepth(vr_vd[[sampleid.t]]),N_totalDepth=totalDepth(vr_vd[[sampleid.n]]), 
                                                N_refDepth=refDepth(vr_vd[[sampleid.n]]), N_altDepth=altDepth(vr_vd[[sampleid.n]]), stringsAsFactors=F))
    
    gr_vd <- unique(gr_vd)
    gr_vd$FILTER=vcf_vd[snv.vd.index]@fixed$FILTER
    
  } else {
    print("Warning: There are no passed SNVs from VarDict.")
    vd.na=T
  }
  
  ## merge 4 vcfs and meta data
  print("Merging data")
  m2.cols = c(names(vcf_m2@info@listData),'REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
  f.cols = c(names(vcf_f@info@listData),'REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
  vs.cols = c(names(vcf_vs@info@listData),'REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
  vd.cols = c(names(vcf_vd@info@listData),'REF','ALT','T_totalDepth','T_refDepth','T_altDepth','N_totalDepth','N_refDepth','N_altDepth','FILTER')
  
  names_m2= paste(m2.cols, "_Mutect2", sep="")
  names_f= paste(f.cols, "_Freebayes", sep="")
  names_vs= paste(vs.cols, "_Varscan", sep="")
  names_vd= paste(vd.cols, "_Vardict", sep="")
  
  meta_data=data.frame(snv_pass)[,1:3]
  meta_data[c(names_m2,names_f,names_vs,names_vd)]=NA
  
  #do not merge when index=0, caller.na=T
  if (m2.na==F) {
    meta_data[Biostrings::match(gr_m2, snv_pass), names_m2]=data.frame(mcols(gr_m2)[,m2.cols])
  }
  if (f.na==F) {
    meta_data[Biostrings::match(gr_f, snv_pass), names_f]=data.frame(mcols(gr_f)[,f.cols])
  }
  if (vs.na==F) {
    meta_data[Biostrings::match(gr_vs, snv_pass), names_vs]=data.frame(mcols(gr_vs)[,vs.cols])
  }
  if (vd.na==F) {
    meta_data[Biostrings::match(gr_vd, snv_pass), names_vd]=data.frame(mcols(gr_vd)[,vd.cols])
  }
  
  # extract reference allele
  ref=meta_data[,c("REF_Mutect2", "REF_Freebayes", "REF_Varscan", "REF_Vardict")] 
  ref.ind=which(!is.na(ref), arr.ind=T) # find non-NA alleles
  ref.ind=ref.ind[order(ref.ind[,1]),] #order array indices by row number
  ref.ind=ref.ind[!duplicated(ref.ind[,1]),] #get first non-NA value for each row
  meta_data$REF=ref[ref.ind]
  
  # extract alternate allele
  alt=meta_data[,c("ALT_Mutect2", "ALT_Freebayes", "ALT_Varscan", "ALT_Vardict")]
  suppressWarnings(alt$ALT_Mutect2[(!is.na(alt$ALT_Mutect2) & nchar(alt$ALT_Mutect2)!=1)] <- substrRight(alt$ALT_Mutect2[(!is.na(alt$ALT_Mutect2) & nchar(alt$ALT_Mutect2)!=1)],1))
  suppressWarnings(alt$ALT_Freebayes[(!is.na(alt$ALT_Freebayes) & nchar(alt$ALT_Freebayes)!=1)] <- substrRight(alt$ALT_Freebayes[(!is.na(alt$ALT_Freebayes) & nchar(alt$ALT_Freebayes)!=1)],1))
  suppressWarnings(alt$ALT_Varscan[(!is.na(alt$ALT_Varscan) & nchar(alt$ALT_Varscan)!=1)] <- substrRight(alt$ALT_Varscan[(!is.na(alt$ALT_Varscan) & nchar(alt$ALT_Varscan)!=1)],1))
  suppressWarnings(alt$ALT_Vardict[(!is.na(alt$ALT_Vardict) & nchar(alt$ALT_Vardict)!=1)] <- substrRight(alt$ALT_Vardict[(!is.na(alt$ALT_Vardict) & nchar(alt$ALT_Vardict)!=1)],1))
  suppressWarnings(alt$ALT_Strelka2[(!is.na(alt$ALT_Strelka2) & nchar(alt$ALT_Strelka2)!=1)] <- substrRight(alt$ALT_Strelka2[(!is.na(alt$ALT_Strelka2) & nchar(alt$ALT_Strelka2)!=1)],1))
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
  meta_data$FILTER_Mutect2[meta_data$FILTER_Mutect2 != "PASS" & meta_data$FILTER_Mutect2 != "MinAF"] <- FALSE
  meta_data$FILTER_Mutect2[is.na(meta_data$FILTER_Mutect2)] <- FALSE
  meta_data$FILTER_Mutect2[meta_data$FILTER_Mutect2 == "PASS"] <- TRUE
  meta_data$FILTER_Mutect2[meta_data$FILTER_Mutect2 == "MinAF"] <- TRUE
  meta_data$FILTER_Mutect2 <- as.logical(meta_data$FILTER_Mutect2)
  
  meta_data$FILTER_Freebayes[meta_data$FILTER_Freebayes != "PASS"] <- FALSE
  meta_data$FILTER_Freebayes[is.na(meta_data$FILTER_Freebayes)] <- FALSE
  meta_data$FILTER_Freebayes[meta_data$FILTER_Freebayes == "PASS"] <- TRUE
  meta_data$FILTER_Freebayes <- as.logical(meta_data$FILTER_Freebayes)
  
  meta_data$FILTER_Vardict[meta_data$FILTER_Vardict != "PASS"] <- FALSE
  meta_data$FILTER_Vardict[is.na(meta_data$FILTER_Vardict)] <- FALSE
  meta_data$FILTER_Vardict[meta_data$FILTER_Vardict == "PASS"] <- TRUE
  meta_data$FILTER_Vardict <- as.logical(meta_data$FILTER_Vardict)
  
  meta_data$FILTER_Varscan[meta_data$FILTER_Varscan != "PASS"] <- FALSE
  meta_data$FILTER_Varscan[is.na(meta_data$FILTER_Varscan)] <- FALSE
  meta_data$FILTER_Varscan[meta_data$FILTER_Varscan == "PASS"] <- TRUE
  meta_data$FILTER_Varscan <- as.logical(meta_data$FILTER_Varscan)
  
  # get all ref/alternates
  meta_data$REF_MFVdVs<-paste(meta_data$REF_Mutect2,meta_data$REF_Freebayes,meta_data$REF_Vardict,meta_data$REF_Varscan,meta_data$REF_Strelka2, sep ="/")
  meta_data$ALT_MFVdVs<-paste(meta_data$ALT_Mutect2,meta_data$ALT_Freebayes,meta_data$ALT_Vardict,meta_data$ALT_Varscan,meta_data$ALT_Strelka2, sep ="/")
  
  # sample name
  meta_data$Sample_Name <- sampleid.t
  
  #chromosome tag
  meta_data$seqnames = gsub('chr','',meta_data$seqnames)
  
  meta_data[,names.features] <- lapply (meta_data[,names.features], as.numeric)
  
  # keep important columns and rename columns
  feature.cols = c(paste0('m2_',mutect2.features),
                   paste0('f_',freebayes.features),
                   paste0('vs_',varscan.features),
                   paste0('vd_',vardict.features))
  
  parse_snv <- meta_data[,c("seqnames","start","end","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                            "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
                            names.features)]
  colnames(parse_snv) <- c("Chr","START_POS_REF","END_POS_REF","REF","ALT","REF_MFVdVs","ALT_MFVdVs","Sample_Name",
                           "FILTER_Mutect2","FILTER_Freebayes","FILTER_Vardict","FILTER_Varscan",
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
  
  write.arff(parse_snv_with_truth, file = paste0(getwd(), "/", folder, "test.arff"), eol = "\n")
  return(parse_snv_with_truth)
}

setwd("C:/Users/mrlim/OneDrive/Desktop/NUS/Y2S2/CS4220/Data2")
parse_snv_with_truth = writeArff("real1")
 