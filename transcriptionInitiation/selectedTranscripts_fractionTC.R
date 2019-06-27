
theme_ameres <- function (type) { 
  
  types_plot = c("boxplot","barplot","violin")
  if(type %in% types_plot ==T)
  {
    
    if(type == "boxplot"){
      theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt"),size = 15),axis.text.y = element_text(margin=margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    if(type == "barplot"){
      theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt"),size = 15, hjust = 1),axis.text.y = element_text(margin=margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    
  }
  
}


sampleInfo = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt",sep="\t",stringsAsFactors = F)
# errorRates = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019//errorRates_predicted_observed.txt",header = T,stringsAsFactors = F)
# errorRates_incubation = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019//errorRates_predicted_observed_incubation.txt",header = T,stringsAsFactors = F)


#Inj_R1 = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/count_revised/combinedFile_CAAGCA.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_tcount.tsv",sep="\t",stringsAsFactors = F,header = T)
Inj_R2 = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019//count_revised///combinedFile_TAGGCT.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_tcount.tsv",sep="\t",stringsAsFactors = F,header = T)
Inj_R3 = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/count_revised/combinedFile_ACGTCT.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_tcount.tsv",sep="\t",stringsAsFactors = F,header = T)

#Inc_R1 = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/count_revised/combinedFile_GTCAGG.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_tcount.tsv",sep="\t",stringsAsFactors = F,header = T)
Inc_R2 = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/count_revised/combinedFile_CGCAAC.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_tcount.tsv",sep="\t",stringsAsFactors = F,header = T)
Inc_R3 = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/count_revised/combinedFile_AATGAA.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_tcount.tsv",sep="\t",stringsAsFactors = F,header = T)

UntreatedSample = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/count_revised//combinedFile_TCAGGA.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_tcount.tsv",sep="\t",stringsAsFactors = F,header = T)


Inc = list(Inc_R2,Inc_R3)
Inj = list(Inj_R2,Inj_R3)


returnFractionTC = function(list_samples,nameGene){
  Inj_ReadCount = do.call(rbind,lapply(list_samples,function(x) x %>%  dplyr::filter(Name == nameGene) %>% dplyr::select(ReadCount)))
  Inj_TCReadCount = do.call(rbind,lapply(list_samples,function(x) x %>%  dplyr::filter(Name == nameGene) %>% dplyr::select(TcReadCount)))
  fractionTC = Inj_TCReadCount$TcReadCount/Inj_ReadCount$ReadCount
  fractionTc_df = data.frame(Name = nameGene, replicate = c("R2","R3"),fractionTC = fractionTC)
  return(fractionTc_df)
}



zygoticGenes = c("vox","vent","mir430_site1","krt18")
zygoticGenes_fraction_Inj = vector("list",length(zygoticGenes))
zygoticGenes_fraction_Inc = vector("list",length(zygoticGenes))



for(i in 1:length(zygoticGenes)){
  zygoticGenes_fraction_Inj[[i]] = returnFractionTC(Inj,zygoticGenes[i])
  zygoticGenes_fraction_Inc[[i]] = returnFractionTC(Inc,zygoticGenes[i])
  }

zygoticGenes_Inj = do.call(rbind.data.frame,zygoticGenes_fraction_Inj)
zygoticGenes_Inj = zygoticGenes_Inj[-which(zygoticGenes_Inj$fractionTC == "NaN"),]

zygoticGenes_Inc = do.call(rbind.data.frame,zygoticGenes_fraction_Inc)
zygoticGenes_Inc = zygoticGenes_Inc[-which(zygoticGenes_Inc$fractionTC == "NaN"),]

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-plyr::ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}


DataSum_Zygotic_Inj = data_summary(data = zygoticGenes_Inj,varname = "fractionTC",groupnames = c("Name") ) %>% dplyr::mutate(type = "Injection")
DataSum_Zygotic_Inc = data_summary(data = zygoticGenes_Inc,varname = "fractionTC",groupnames = c("Name") ) %>% dplyr::mutate(type = "Incubation") 
UntreatedSample_zygotic = UntreatedSample[UntreatedSample$Name %in% zygoticGenes,] %>% dplyr::filter(ReadsCPM !=0) %>% dplyr::mutate(fractionTC = TcReadCount/ReadCount) %>%dplyr::select(Name,fractionTC) %>%
  dplyr::mutate(sd =0) %>% dplyr::mutate(type = "-4SU")

zygoticGenes_Inj = zygoticGenes_Inj %>% dplyr::mutate(type = "Injection")
zygoticGenes_Inc = zygoticGenes_Inc %>% dplyr::mutate(type = "Incubation")
UntreatedSample_zygotic = UntreatedSample_zygotic %>% dplyr::mutate(type = "Untreated",replicate = "R1") %>% dplyr::select(- 'sd') %>% dplyr::select(c('Name','replicate','fractionTC','type'))
totalData = rbind(zygoticGenes_Inj,zygoticGenes_Inc,UntreatedSample_zygotic)

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/transcriptionInitiation/exampleGenes_fractionTC.pdf",height = 4)

p = ggpubr::ggscatter(totalData,x='Name',y='fractionTC',color = 'type',palette = 'Dark2',size=4) + theme_ameres(type = "barplot")  + ylim(c(0,1))


#### same for maternal 



maternalGenes = c("eml2","smurf1","ccna1","cldnd")
maternalGenes_fraction_Inj = vector("list",length(maternalGenes))
maternalGenes_fraction_Inc = vector("list",length(maternalGenes))



for(i in 1:length(maternalGenes)){
  maternalGenes_fraction_Inj[[i]] = returnFractionTC(Inj,maternalGenes[i])
  maternalGenes_fraction_Inc[[i]] = returnFractionTC(Inc,maternalGenes[i])
}

maternalGenes_Inj = do.call(rbind.data.frame,maternalGenes_fraction_Inj)

maternalGenes_Inc = do.call(rbind.data.frame,maternalGenes_fraction_Inc)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-plyr::ddply(data, groupnames, .fun=summary_func,
                        varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}


DataSum_Maternal_Inj = data_summary(data = maternalGenes_Inj,varname = "fractionTC",groupnames = c("Name") ) %>% dplyr::mutate(type = "Injection")
DataSum_Maternal_Inc = data_summary(data = maternalGenes_Inc,varname = "fractionTC",groupnames = c("Name") ) %>% dplyr::mutate(type = "Incubation") 
UntreatedSample_maternal = UntreatedSample[UntreatedSample$Name %in% maternalGenes,] %>% dplyr::filter(ReadsCPM !=0) %>% dplyr::mutate(fractionTC = TcReadCount/ReadCount) %>%dplyr::select(Name,fractionTC) %>%
  dplyr::mutate(sd =0) %>% dplyr::mutate(type = "-4SU")


maternalGenes_Inj = maternalGenes_Inj %>% dplyr::mutate(type = "Injection")
maternalGenes_Inc = maternalGenes_Inc %>% dplyr::mutate(type = "Incubation")
UntreatedSample_maternal = UntreatedSample_maternal %>% dplyr::mutate(type = "Untreated",replicate = "R1") %>% dplyr::select(- 'sd') %>% dplyr::select(c('Name','replicate','fractionTC','type'))
totalData = rbind(maternalGenes_Inj,maternalGenes_Inc,UntreatedSample_maternal)


q = ggpubr::ggscatter(totalData,x='Name',y='fractionTC',color = 'type',palette = 'Dark2',size=4) + theme_ameres(type = "barplot")   + ylim(c(0,1))



# RColorBrewer::brewer.pal(n = 3,name = "Dark2")
# 
# p = ggplot(totalData, aes(x=Name, y=fractionTC,fill=type)) + ylim(c(0,1))+
#   geom_bar(stat="identity", position=position_dodge() )+
#   geom_errorbar(aes(ymin=fractionTC-sd, ymax=fractionTC+sd), width=.2,
#                 position=position_dodge(.9)) + theme_bw() + theme_ameres(type = "barplot")
# p + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                          panel.background = element_blank(), axis.line = element_line(colour = "black")) 

print(p)
print(q)


dev.off()
