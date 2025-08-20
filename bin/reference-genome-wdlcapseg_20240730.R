#################################################################
##                      READ IN LIBRARIES                      ##
#################################################################
library(data.table)
library(parallel)
library(optparse)
library(tidyverse)
library(dplyr)

##################################################################
##                DEFINE INPUT OPTIONS AND FLAGS                ##
##################################################################
option_list <- list(
  make_option(c("-s", "--segfile"), type = "character", default = NA,
              help = "Path to iconicc segfile output"),
  make_option(c("-c", "--processed_cts"), type = "character", default = NA,
              help = "Path to processed allelic counts file"),
  make_option(c("-i", "--participant_id"), type = "character", default = NA,
              help = "Unique sample identiifer")
  
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

iconicc_seg <- fread(opt$segfile)
processed_cts <- fread(opt$processed_cts)
id <- opt$participant_id


canonical_contigs <- c(paste0("chr",1:22),"chrX","chrY")
assign_group <- function(itter,dt,segments_CNA,seg.source){
  row_c <- dt[itter,] 
  chr <- row_c$CONTIG 
  segments_info <- subset(segments_CNA,chrom==chr)
  seg_group <- segments_info[which((as.numeric(row_c$START)>=loc.start)&(as.numeric(row_c$END)<=(loc.end+24999))),]$SegmentID
  if (identical(integer(0),seg_group)){
    seg_group <- NA
  }
  return(row_c[,paste0(seg.source,'_segment'):=seg_group])
}

iconicc_to_capseg <- function(segfile,full_dt){

  segfile[,SegmentID := 1:nrow(segfile)]
  alleliccapseg <- segfile
  full_dt <- full_dt[order(factor(CONTIG,levels = canonical_contigs),as.numeric(START))]
  full_dt <- rbindlist(mclapply(1:nrow(full_dt),assign_group,full_dt,segfile,seg.source='final',mc.cores = 7))
  full_dt[,tangent_b10 := (2^log2_tangent)]
  alleliccapseg[,tau:=2*(2^seg.mean)]
  sd.tau <- aggregate(full_dt$tangent_b10,by=list(full_dt$final_segment),FUN = function(x){return (sd(x,na.rm=T))})
  sd.tau <- as.data.table(sd.tau)
  colnames(sd.tau) <- c('SegmentID','sigma.tau')
  alleliccapseg <- merge(alleliccapseg,sd.tau,by='SegmentID')
  
  ##figure out threshold for setting f to 0.5;
  ###I'll first determine the 90% percentile of AFMIN for the rows in full_dt that we could have used to call the segments
  maf90threshold <- quantile(full_dt[snp_count>5,]$AFMIN,0.9) #shoudl I set this value or let it change each time
  ###what if we had a sample that was riddled with allelic imbalance?
  #correct artificial 0 and 1 for AFMIN/AFMAX so that calculations here are not skewed
  #full_dt[,AFMIN:=ifelse(snp_count==0,NA,AFMIN)]
  #full_dt[,AFMAX:=ifelse(snp_count==0,NA,AFMAX)]
  #or instead correct all things with snp_count < 5 to consider only nodes that were incorporated into defining the segments
  full_dt[,AFMIN:=ifelse(snp_count<5,NA,AFMIN)]
  full_dt[,AFMAX:=ifelse(snp_count<5,NA,AFMAX)]
  
  min.af <- aggregate(full_dt$AFMIN,by=list(full_dt$final_segment),FUN = function(x){return (c(afmin.avg = mean(x,na.rm=T), sigma.min=sd(x,na.rm=T)))})
  min.af <- as.data.table(min.af)
  colnames(min.af) <- c('SegmentID','afmin.avg','sigma.min')
  
  ####do the same for the major allele but only retain sd info since the avg will be 1-afmin
  maj.af <- aggregate(full_dt$AFMAX,by=list(full_dt$final_segment),FUN = function(x){return (sd(x,na.rm=T))})
  maj.af <- as.data.table(maj.af)
  colnames(maj.af) <- c('SegmentID','sigma.maj')
  
  alleliccapseg <- Reduce(merge,list(alleliccapseg,min.af,maj.af))
  alleliccapseg[,f:=ifelse(afmin.avg>=maf90threshold,0.5,afmin.avg)]
  alleliccapseg[,mu.minor:= f*tau]
  #alleliccapseg[,sigma.minor:= mu.minor*(sqrt((sigma.tau/tau)^2 + (sigma.min/f)^2))] 
  alleliccapseg[,sigma.minor:= (sqrt((sigma.tau)^2 + (sigma.min)^2))] #following the error propagation on the github
  alleliccapseg[,mu.major:= (1-f)*tau]
  #alleliccapseg[,sigma.major:= mu.major*(sqrt((sigma.tau/tau)^2 + (sigma.maj/1-f)^2))]
  alleliccapseg[,sigma.major:= (sqrt((sigma.tau)^2 + (sigma.maj)^2))]
  
  ###collect no. of snps
  snps.agg <- aggregate(full_dt$snp_count,by=list(full_dt$final_segment),FUN = sum)
  snps.agg <- as.data.table(snps.agg)
  colnames(snps.agg) <- c('SegmentID','n_hets')
  alleliccapseg <- merge(alleliccapseg,snps.agg,by='SegmentID')
  
  #tau=2&afmax=1 for CNLOH
  #user determined thresholds for afmin and tau shift to call LOH
  alleliccapseg[,SegLabelCNLOH:=ifelse((tau>=1.8&tau<=2.2)&afmin.avg<=0.05,1,0)] #not necessary for absolute so a guess here is fine too
  alleliccapseg[,length:=loc.end-loc.start]
  
  
  #select relavant columns
  alleliccapseg <- alleliccapseg[,c('chrom','loc.start','loc.end','num.mark','length','n_hets',
                                    'f','tau','sigma.tau','mu.minor','sigma.minor','mu.major',
                                    'sigma.major','SegLabelCNLOH')]
  colnames(alleliccapseg) <- c("Chromosome","Start.bp","End.bp","n_probes","length","n_hets","f","tau",          
                               "sigma.tau","mu.minor","sigma.minor","mu.major","sigma.major","SegLabelCNLOH")
  alleliccapseg[,Chromosome:=gsub('chr','',Chromosome)]
  alleliccapseg <- alleliccapseg %>% drop_na(tau) #drop rows that segment info is essentially missing
  
  write.table(alleliccapseg, file = paste0(id,'.capseg.txt'),sep = '\t',col.names = T,row.names = F,quote = F)
  return(alleliccapseg)
}

iconicc_to_capseg(iconicc_seg,processed_cts)
