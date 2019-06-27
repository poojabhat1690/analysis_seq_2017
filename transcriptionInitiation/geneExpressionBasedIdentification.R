#### 

##########################
### transcripts that are continuously expressed from stage 1 -TP1
### transcripts that are continuously expressed from stage 2 -TP2 .... to stage 9. 


library(ggplot2)
library(reshape)
library(dplyr)
library(biomaRt)
library(topGO)
library(ggrepel)

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

allRPMs = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/dataTables_expandedCounting/RPM_allCws.txt",
                     sep= "\t",stringsAsFactors = F,header = T)

TCconversions = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/dataTables_expandedCounting/ConversionRates_allCws.txt",
                           sep= "\t",stringsAsFactors = F,header = T)


splitReplicates = function(dataFrameToSplit,condition,metadata_add){
  dataFrameToSplit_condition = dataFrameToSplit %>% select_if(grepl(condition,names(.)))
  dataFrameToSplit_condition_R1 = dataFrameToSplit_condition %>% select_if(grepl("R1",names(.)))
  dataFrameToSplit_condition_R2 = dataFrameToSplit_condition %>% select_if(grepl("R2",names(.)))
  dataFrameToSplit_condition_R3 = dataFrameToSplit_condition %>% select_if(grepl("R3",names(.)))
  mean_repl = (dataFrameToSplit_condition_R1+dataFrameToSplit_condition_R2+dataFrameToSplit_condition_R3)/3
  dataFrameToSplit_condition_R1 = cbind.data.frame(dataFrameToSplit_condition_R1,metadata_add)
  dataFrameToSplit_condition_R2 = cbind.data.frame(dataFrameToSplit_condition_R2,metadata_add)
  dataFrameToSplit_condition_R3 = cbind.data.frame(dataFrameToSplit_condition_R3,metadata_add)
  mean_repl = cbind.data.frame(mean_repl,metadata_add)
  splitReplicates = list(dataFrameToSplit_condition_R1,dataFrameToSplit_condition_R2,dataFrameToSplit_condition_R3,mean_repl)
  names(splitReplicates) = c("R1","R2","R3","mean")
  return(splitReplicates)
}

annotation_data = allRPMs %>% dplyr::select(chr,start,end,name,strand )
splitData_RPMs = splitReplicates(dataFrameToSplit =allRPMs,condition = "Inj",metadata_add = annotation_data )
splitData_conversions = splitReplicates(dataFrameToSplit =TCconversions,condition = "Inj",metadata_add = annotation_data )

splitData_RPMs_means =splitData_RPMs$mean
splitData_conversions_means = splitData_conversions$mean %>% dplyr::select( dplyr::contains('Inj'))
splitData_conversions_untreated = TCconversions %>% dplyr::select( dplyr::contains('Untreated'))
conversions_backgroundSubtracted = as.matrix(splitData_conversions_means - splitData_conversions_untreated )
conversions_backgroundSubtracted[which(conversions_backgroundSubtracted <0)] <- 0
conversions_backgroundSubtracted = as.data.frame(conversions_backgroundSubtracted)

colnames(conversions_backgroundSubtracted) = paste0("conversoins_",colnames(conversions_backgroundSubtracted))
conversionAndRPMs = cbind.data.frame(splitData_RPMs_means,conversions_backgroundSubtracted)


classifiedTranscripts = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11//dataTables/countingWindows_classified_1.txt",sep="\t",stringsAsFactors = F,header = T)
classifiedTranscripts = classifiedTranscripts %>% dplyr::select(V1,V2,V3,V4,V5,V6,category) %>% dplyr::mutate(chr = V1, start = V2, 
                                                                                                              end = V3 , name = V4, strand = V6)   %>%   dplyr::select(-c(V1,V2,V3,V4,V5,V6))
classifiedTranscripts = classifiedTranscripts %>% dplyr::filter(category == "Z" | category == "MZ")
conversionAndRPMs = conversionAndRPMs %>%  dplyr::filter(Inj_R1_TP9 > 5 | Inj_R1_TP8 > 5 | Inj_R1_TP7 > 5 | Inj_R1_TP6 >5 | Inj_R1_TP5 > 5 | Inj_R1_TP4 > 5 | Inj_R1_TP3 > 5 | Inj_R1_TP2 > 5 | Inj_R1_TP1 >5)

conversionAndRPMs = plyr::join(classifiedTranscripts,conversionAndRPMs)

##### getting the transcripts expressed at each stage


expressedFromTP1 = conversionAndRPMs %>% dplyr::filter(Inj_R1_TP1 > 0 & Inj_R1_TP7 > 0 & Inj_R1_TP8 > 0 & Inj_R1_TP9 > 0 ) %>% dplyr::mutate(geneXpression = "TP1")
expressedFromTP2 = conversionAndRPMs %>% dplyr::filter(Inj_R1_TP1 == 0 & Inj_R1_TP2 > 0  & Inj_R1_TP7 > 0 & Inj_R1_TP8 > 0 & Inj_R1_TP9 > 0 ) %>% dplyr::mutate(geneXpression = "TP2")
expressedFromTP3 = conversionAndRPMs %>% dplyr::filter(Inj_R1_TP1 == 0 & Inj_R1_TP2 == 0 & Inj_R1_TP3 > 0  & Inj_R1_TP7 > 0 & Inj_R1_TP8 > 0 & Inj_R1_TP9 > 0 ) %>% dplyr::mutate(geneXpression = "TP3")
expressedFromTP4 = conversionAndRPMs %>% dplyr::filter(Inj_R1_TP1 == 0 & Inj_R1_TP2 == 0 & Inj_R1_TP3 == 0 & Inj_R1_TP4 > 0  & Inj_R1_TP7 > 0 & Inj_R1_TP8 > 0 & Inj_R1_TP9 > 0 ) %>% dplyr::mutate(geneXpression = "TP4")
expressedFromTP5 = conversionAndRPMs %>% dplyr::filter(Inj_R1_TP1 == 0 & Inj_R1_TP2 == 0 & Inj_R1_TP3 == 0 & Inj_R1_TP4 == 0 & Inj_R1_TP5 > 0 & Inj_R1_TP7 > 0 & Inj_R1_TP8 > 0 & Inj_R1_TP9 > 0 ) %>% dplyr::mutate(geneXpression = "TP5")
expressedFromTP6 = conversionAndRPMs %>% dplyr::filter(Inj_R1_TP1 == 0 & Inj_R1_TP2 == 0 & Inj_R1_TP3 == 0 & Inj_R1_TP4 == 0 & Inj_R1_TP5 == 0 & Inj_R1_TP6 > 0 & Inj_R1_TP7 > 0 & Inj_R1_TP8 > 0 & Inj_R1_TP9 > 0 ) %>% dplyr::mutate(geneXpression = "TP6")
expressedFromTP7 = conversionAndRPMs %>% dplyr::filter(Inj_R1_TP1 == 0 & Inj_R1_TP2 == 0 & Inj_R1_TP3 == 0 & Inj_R1_TP4 == 0 & Inj_R1_TP5 == 0 & Inj_R1_TP6 == 0 & Inj_R1_TP7 > 0 & Inj_R1_TP8 > 0 & Inj_R1_TP9 > 0 ) %>% dplyr::mutate(geneXpression = "TP7")
expressedFromTP8 = conversionAndRPMs %>% dplyr::filter(Inj_R1_TP1 == 0 & Inj_R1_TP2 == 0 & Inj_R1_TP3 == 0 & Inj_R1_TP4 == 0 & Inj_R1_TP5 == 0 & Inj_R1_TP6 == 0 & Inj_R1_TP7 == 0 & Inj_R1_TP8 > 0 & Inj_R1_TP9 > 0 ) %>% dplyr::mutate(geneXpression = "TP8")
expressedFromTP9 = conversionAndRPMs %>% dplyr::filter(Inj_R1_TP1 == 0 & Inj_R1_TP2 == 0 & Inj_R1_TP3 == 0 & Inj_R1_TP4 == 0 & Inj_R1_TP5 == 0 & Inj_R1_TP6 == 0 & Inj_R1_TP7 == 0 & Inj_R1_TP8 == 0 & Inj_R1_TP9 > 0 ) %>% dplyr::mutate(geneXpression = "TP9")

##### but these show de novo transcription from different time points... 
totalCapturedConversions = rbind.data.frame(expressedFromTP1,expressedFromTP2,expressedFromTP3,expressedFromTP4,expressedFromTP5,expressedFromTP6,expressedFromTP7,expressedFromTP8,expressedFromTP9)

###### now i want to add infos about when TC conversions are detected for these transcripts. 
TCconversionsAtTP1 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 > 0) %>% dplyr::mutate(TCconversionTime = "TP1")
TCconversionsAtTP2 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 == 0 & conversoins_Inj_R1_TP2 > 0) %>% dplyr::mutate(TCconversionTime = "TP2")
TCconversionsAtTP3 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 == 0 & conversoins_Inj_R1_TP2 == 0  & conversoins_Inj_R1_TP3> 0 ) %>% 
  dplyr::mutate(TCconversionTime = "TP3")
TCconversionsAtTP4 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 == 0 & conversoins_Inj_R1_TP2 == 0  & conversoins_Inj_R1_TP3 == 0 & conversoins_Inj_R1_TP4 > 0 ) %>% 
  dplyr::mutate(TCconversionTime = "TP4")
TCconversionsAtTP5 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 == 0 & conversoins_Inj_R1_TP2 == 0  & conversoins_Inj_R1_TP3 == 0 & conversoins_Inj_R1_TP4 == 0 & conversoins_Inj_R1_TP5 > 0 ) %>% 
  dplyr::mutate(TCconversionTime = "TP5")
TCconversionsAtTP6 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 == 0 & conversoins_Inj_R1_TP2 == 0  & conversoins_Inj_R1_TP3 == 0 & conversoins_Inj_R1_TP4 == 0 & conversoins_Inj_R1_TP5 == 0  & conversoins_Inj_R1_TP6 > 0 ) %>% 
  dplyr::mutate(TCconversionTime = "TP6")
TCconversionsAtTP7 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 == 0 & conversoins_Inj_R1_TP2 == 0  & conversoins_Inj_R1_TP3 == 0 & conversoins_Inj_R1_TP4 == 0 & conversoins_Inj_R1_TP5 == 0  & conversoins_Inj_R1_TP6 == 0 & conversoins_Inj_R1_TP7 > 0 ) %>% 
  dplyr::mutate(TCconversionTime = "TP7")
TCconversionsAtTP8 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 == 0 & conversoins_Inj_R1_TP2 == 0  & conversoins_Inj_R1_TP3 == 0 & conversoins_Inj_R1_TP4 == 0 & conversoins_Inj_R1_TP5 == 0  & conversoins_Inj_R1_TP6 == 0 & conversoins_Inj_R1_TP7 == 0 & conversoins_Inj_R1_TP8 > 0 ) %>% 
  dplyr::mutate(TCconversionTime = "TP8")
TCconversionsAtTP9 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 == 0 & conversoins_Inj_R1_TP2 == 0  & conversoins_Inj_R1_TP3 == 0 & conversoins_Inj_R1_TP4 == 0 & conversoins_Inj_R1_TP5 == 0  & conversoins_Inj_R1_TP6 == 0 & conversoins_Inj_R1_TP7 == 0 & conversoins_Inj_R1_TP8 == 0 & conversoins_Inj_R1_TP9 > 0 ) %>% 
  dplyr::mutate(TCconversionTime = "TP9")


expression_TCtiming = rbind.data.frame(TCconversionsAtTP1,TCconversionsAtTP2,TCconversionsAtTP3,TCconversionsAtTP4,TCconversionsAtTP5,TCconversionsAtTP6,TCconversionsAtTP7,TCconversionsAtTP8,TCconversionsAtTP9)
expression_TCtiming_Z = expression_TCtiming %>% dplyr::filter(category == "Z")
stageSpecificExpression_Z = as.data.frame(table(expression_TCtiming_Z$TCconversionTime)) %>% dplyr::mutate(type = "Z") 

expression_TCtiming_MZ = expression_TCtiming %>% dplyr::filter(category == "MZ")
stageSpecificExpression_MZ = as.data.frame(table(expression_TCtiming_MZ$TCconversionTime)) %>% dplyr::mutate(type = "MZ") 
expression_TCtiming_Z_write = expression_TCtiming_Z %>% dplyr::select(-'category')

#write.table(expression_TCtiming_Z_write,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/countinWindows_expressedAccordingToTime.bed",sep="\t",quote = F,row.names = F,col.names = F)

totalIncreaseOfTCconversions = rbind.data.frame(stageSpecificExpression_Z,stageSpecificExpression_MZ)
p = ggpubr::ggbarplot(totalIncreaseOfTCconversions,x='Var1',y='Freq',fill='type',palette = 'grey',position = position_dodge(0.9),label = 'Freq',xlab = 'Timepoint',ylab='number of transcripts')  + theme_ameres(type = 'barplot')

#### now i want to check how many transcripts actually have continuous SLAMseq signal. 



conversionAndRPMs_filterConversions_atTP1 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 > 0 & conversoins_Inj_R1_TP2>0 & conversoins_Inj_R1_TP3>0 & conversoins_Inj_R1_TP4>0 & conversoins_Inj_R1_TP5>0 & conversoins_Inj_R1_TP6>0 & conversoins_Inj_R1_TP7>0 & conversoins_Inj_R1_TP8>0 & conversoins_Inj_R1_TP9>0)  %>% dplyr::mutate(timepoint = "TP1")
conversionAndRPMs_filterConversions_atTP2 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 == 0 & conversoins_Inj_R1_TP2>0 & conversoins_Inj_R1_TP3>0 & conversoins_Inj_R1_TP4>0 & conversoins_Inj_R1_TP5>0 & conversoins_Inj_R1_TP6>0 & conversoins_Inj_R1_TP7>0 & conversoins_Inj_R1_TP8>0 & conversoins_Inj_R1_TP9>0)  %>% dplyr::mutate(timepoint = "TP2")
conversionAndRPMs_filterConversions_atTP3 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 == 0 & conversoins_Inj_R1_TP2 == 0 & conversoins_Inj_R1_TP3>0 & conversoins_Inj_R1_TP4>0 & conversoins_Inj_R1_TP5>0 & conversoins_Inj_R1_TP6>0 & conversoins_Inj_R1_TP7>0 & conversoins_Inj_R1_TP8>0 & conversoins_Inj_R1_TP9>0)  %>% dplyr::mutate(timepoint = "TP3")
conversionAndRPMs_filterConversions_atTP4 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 == 0 & conversoins_Inj_R1_TP2 == 0 & conversoins_Inj_R1_TP3==0 & conversoins_Inj_R1_TP4>0 & conversoins_Inj_R1_TP5>0 & conversoins_Inj_R1_TP6>0 & conversoins_Inj_R1_TP7>0 & conversoins_Inj_R1_TP8>0 & conversoins_Inj_R1_TP9>0)  %>% dplyr::mutate(timepoint = "TP4")
conversionAndRPMs_filterConversions_atTP5 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 == 0 & conversoins_Inj_R1_TP2 == 0 & conversoins_Inj_R1_TP3==0 & conversoins_Inj_R1_TP4==0 & conversoins_Inj_R1_TP5>0  & conversoins_Inj_R1_TP6>0 & conversoins_Inj_R1_TP7>0 & conversoins_Inj_R1_TP8>0 & conversoins_Inj_R1_TP9>0)  %>% dplyr::mutate(timepoint = "TP5")
conversionAndRPMs_filterConversions_atTP6 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 == 0 & conversoins_Inj_R1_TP2 == 0 & conversoins_Inj_R1_TP3==0 & conversoins_Inj_R1_TP4==0 & conversoins_Inj_R1_TP5==0 & conversoins_Inj_R1_TP6>0 & conversoins_Inj_R1_TP7>0 & conversoins_Inj_R1_TP8>0 & conversoins_Inj_R1_TP9>0)  %>% dplyr::mutate(timepoint = "TP6")
conversionAndRPMs_filterConversions_atTP7 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 == 0 & conversoins_Inj_R1_TP2 == 0 & conversoins_Inj_R1_TP3==0 & conversoins_Inj_R1_TP4==0 & conversoins_Inj_R1_TP5==0 & conversoins_Inj_R1_TP6==0 & conversoins_Inj_R1_TP7>0 & conversoins_Inj_R1_TP8>0 & conversoins_Inj_R1_TP9>0) %>% dplyr::mutate(timepoint = "TP7")
conversionAndRPMs_filterConversions_atTP8 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 == 0 & conversoins_Inj_R1_TP2 == 0 & conversoins_Inj_R1_TP3==0 & conversoins_Inj_R1_TP4==0 & conversoins_Inj_R1_TP5==0 & conversoins_Inj_R1_TP6==0 & conversoins_Inj_R1_TP7==0 & conversoins_Inj_R1_TP8>0 & conversoins_Inj_R1_TP9>0)  %>% dplyr::mutate(timepoint = "TP8")
conversionAndRPMs_filterConversions_atTP9 = totalCapturedConversions %>% dplyr::filter(conversoins_Inj_R1_TP1 == 0 & conversoins_Inj_R1_TP2 == 0 & conversoins_Inj_R1_TP3==0 & conversoins_Inj_R1_TP4==0 & conversoins_Inj_R1_TP5==0  & conversoins_Inj_R1_TP6==0 & conversoins_Inj_R1_TP7==0 & conversoins_Inj_R1_TP8==0 & conversoins_Inj_R1_TP9>0)  %>% dplyr::mutate(timepoint = "TP9")

totalCapturedConversions = rbind.data.frame(conversionAndRPMs_filterConversions_atTP9,conversionAndRPMs_filterConversions_atTP8,conversionAndRPMs_filterConversions_atTP7,
                                            conversionAndRPMs_filterConversions_atTP6,conversionAndRPMs_filterConversions_atTP5,conversionAndRPMs_filterConversions_atTP4,
                                            conversionAndRPMs_filterConversions_atTP3,conversionAndRPMs_filterConversions_atTP2,conversionAndRPMs_filterConversions_atTP1)


Continuedexpression_TCtiming_Z = totalCapturedConversions %>% dplyr::filter(category == "Z")
ContinuedstageSpecificExpression_Z = as.data.frame(table(Continuedexpression_TCtiming_Z$timepoint)) %>% dplyr::mutate(type = "Z") 

Continuedexpression_TCtiming_MZ = totalCapturedConversions %>% dplyr::filter(category == "MZ")
ContinuedstageSpecificExpression_MZ = as.data.frame(table(Continuedexpression_TCtiming_MZ$timepoint)) %>% dplyr::mutate(type = "MZ") 

ContinuoustotalIncreaseOfTCconversions = rbind.data.frame(ContinuedstageSpecificExpression_Z,ContinuedstageSpecificExpression_MZ)
q = ggpubr::ggbarplot(ContinuoustotalIncreaseOfTCconversions,x='Var1',y='Freq',fill='type',palette = 'grey',position = position_dodge(0.9),label = 'Freq',xlab = 'Timepoint',ylab='number of transcripts')  + theme_ameres(type = 'barplot')

pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/transcriptionInitiation/numberOfGenes_eachStage.pdf",height = 4,width = 5)
  print(p)
  print(q)
dev.off()


##### plotting TC conversion rates for the 7 purely zygotic transcripts detected at TP1 and TP2. 

zygoticTranscripts_TP1 = expression_TCtiming_Z %>% dplyr::filter(TCconversionTime == "TP1") 
zygoticTranscripts_TP2 = expression_TCtiming_Z %>% dplyr::filter(TCconversionTime == "TP2")


#####################################################################################################################################################################

############### Checking H3K27Ac signal and ATACseq signal at purely zygotic transcripts detected at TP1 and TP2 ###############

#####################################################################################################################################################################

ATACseqSignal_1 = read.table("/Volumes/clustertmp/pooja/ATACseq_fastqDump/countDir/SRR6476749_counts1kb.bed")
ATACseqSignal_2 = read.table("/Volumes/clustertmp/pooja/ATACseq_fastqDump/countDir/SRR6476750_counts1kb.bed")
ATACseqSignal_3 = read.table("/Volumes/clustertmp/pooja/ATACseq_fastqDump/countDir/SRR6476751_counts1kb.bed")
ATACseqSignal_1$RPM = ATACseqSignal_1$V5 * 1000000 /sum(ATACseqSignal_1$V5)
ATACseqSignal_2$RPM = ATACseqSignal_2$V5 * 1000000 /sum(ATACseqSignal_2$V5)
ATACseqSignal_3$RPM = ATACseqSignal_3$V5 * 1000000 /sum(ATACseqSignal_3$V5)


geneAnnotation = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/dr11/ensembl_dr11_Ensembl_Genes_93/transcriptStartsAndEnds_all.txt",sep="\t",header = T,stringsAsFactors = F)
intersect(geneAnnotation$external_gene_name, zygoticTranscripts_TP1$Name)


geneAnnotation_samples_TP1 = geneAnnotation[geneAnnotation$external_gene_name %in% zygoticTranscripts_TP1$Name,]
geneAnnotation_samples_TP1 = geneAnnotation_samples_TP1[!duplicated(geneAnnotation_samples_TP1$external_gene_name) ,]


geneAnnotation_samples_TP2 = geneAnnotation[geneAnnotation$external_gene_name %in% zygoticTranscripts_TP2$Name,]
geneAnnotation_samples_TP2 = geneAnnotation_samples_TP2[!duplicated(geneAnnotation_samples_TP2$external_gene_name) ,]



overlapWithSeq = function(geneAnnotation,seqData){
  geneAnnotation_samples_Granges= makeGRangesFromDataFrame(df = geneAnnotation,keep.extra.columns = T,ignore.strand = F,seqnames.field = 'chromosome_name',start.field = 'transcript_start',end.field = 'transcript_end')
  seqData_Granges  = makeGRangesFromDataFrame(seqData,keep.extra.columns = T,ignore.strand = T,start.field = 'V2',end.field = 'V3',seqnames.field = 'V1')
  geneAnnotationOverlapWithATACseq = findOverlaps(query = geneAnnotation_samples_Granges,subject = seqData_Granges)
  ATACseqOverlap = seqData[subjectHits(geneAnnotationOverlapWithATACseq),]
  geneAnnotationOverlap = geneAnnotation[queryHits(geneAnnotationOverlapWithATACseq),]
  subjectQueryOverlap = cbind.data.frame(ATACseqOverlap,geneAnnotationOverlap)
  a = subjectQueryOverlap %>% dplyr::group_by(external_gene_name) %>% mutate(meanSignal = mean(RPM))
  a = a[!duplicated(a[,c("external_gene_name","meanSignal")]),]
  
  return(a)
}



## TP1
ATACseqOverlap_1 = overlapWithSeq(geneAnnotation = geneAnnotation_samples_TP1,seqData = ATACseqSignal_1)
ATACseqOverlap_2 = overlapWithSeq(geneAnnotation = geneAnnotation_samples_TP1,seqData = ATACseqSignal_2)
ATACseqOverlap_3 = overlapWithSeq(geneAnnotation = geneAnnotation_samples_TP1,seqData = ATACseqSignal_3)

ATACseqOverlap_1$RPM_2 = ATACseqOverlap_2$meanSignal
ATACseqOverlap_1$RPM_3 = ATACseqOverlap_3$meanSignal
ATACseqOverlap_1$RPM_1 = ATACseqOverlap_1$meanSignal
ATACseqOverlap_1$meanRPM = rowMeans(ATACseqOverlap_1[,c("RPM_1","RPM_2","RPM_3")])
#ATACseqOverlap_1 = ATACseqOverlap_1[order(ATACseqOverlap_1$meanRPM),]

ATACseqOverlap_1$geneNumber = c(1:nrow(ATACseqOverlap_1))
ATACseqOverlap_1_rep1 = ATACseqOverlap_1
ggpubr::ggscatter(ATACseqOverlap_1,x='geneNumber',y='meanRPM',label = 'external_gene_name')

### TP2

ATACseqOverlap_1 = overlapWithSeq(geneAnnotation = geneAnnotation_samples_TP2,seqData = ATACseqSignal_1)
ATACseqOverlap_2 = overlapWithSeq(geneAnnotation = geneAnnotation_samples_TP2,seqData = ATACseqSignal_2)
ATACseqOverlap_3 = overlapWithSeq(geneAnnotation = geneAnnotation_samples_TP2,seqData = ATACseqSignal_3)

ATACseqOverlap_1$RPM_2 = ATACseqOverlap_2$meanSignal
ATACseqOverlap_1$RPM_3 = ATACseqOverlap_3$meanSignal
ATACseqOverlap_1$RPM_1 = ATACseqOverlap_1$meanSignal
ATACseqOverlap_1$meanRPM = rowMeans(ATACseqOverlap_1[,c("RPM_1","RPM_2","RPM_3")])
ATACseqOverlap_1 = ATACseqOverlap_1[order(ATACseqOverlap_1$meanRPM),]

ATACseqOverlap_1$geneNumber = c(1:nrow(ATACseqOverlap_1))
ggpubr::ggscatter(ATACseqOverlap_1,x='geneNumber',y='meanRPM',label = 'external_gene_name')


################# i want to do the same with H3k27AC


h3k27ac_1 = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/CHIPseq/coverage/SRR7235547_counts1kb.bed")
h3k27ac_2 = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/CHIPseq/coverage/SRR7235548_counts1kb.bed")

h3k27ac_1$RPM = h3k27ac_1$V5 * 1000000 /sum(h3k27ac_1$V5)
h3k27ac_2$RPM = h3k27ac_2$V5 * 1000000 /sum(h3k27ac_2$V5)


## TP1
h3k27acOverlap_1 = overlapWithSeq(geneAnnotation = geneAnnotation_samples_TP1,seqData = h3k27ac_1)
h3k27acOverlap_2 = overlapWithSeq(geneAnnotation = geneAnnotation_samples_TP1,seqData = h3k27ac_2)
h3k27acOverlap_1$RPM_2 = h3k27acOverlap_2$meanSignal
h3k27acOverlap_1$RPM_1 = h3k27acOverlap_1$meanSignal

h3k27acOverlap_1$meanRPM = rowMeans(h3k27acOverlap_1[,c("RPM_1","RPM_2")])
#h3k27acOverlap_1 = h3k27acOverlap_1[order(h3k27acOverlap_1$meanRPM),]

h3k27acOverlap_1$geneNumber = c(1:nrow(h3k27acOverlap_1))
h3k27acOverlap_1_rep1 = h3k27acOverlap_1
ggpubr::ggscatter(h3k27acOverlap_1,x='geneNumber',y='meanRPM',label = 'external_gene_name')

## TP2

h3k27acOverlap_1 = overlapWithSeq(geneAnnotation = geneAnnotation_samples_TP2,seqData = h3k27ac_1)
h3k27acOverlap_2 = overlapWithSeq(geneAnnotation = geneAnnotation_samples_TP2,seqData = h3k27ac_2)
h3k27acOverlap_1$RPM_2 = h3k27acOverlap_2$meanSignal
h3k27acOverlap_1$RPM_1 = h3k27acOverlap_1$meanSignal

h3k27acOverlap_1$meanRPM = rowMeans(h3k27acOverlap_1[,c("RPM_1","RPM_2")])
#h3k27acOverlap_1 = h3k27acOverlap_1[order(h3k27acOverlap_1$meanRPM),]

h3k27acOverlap_1$geneNumber = c(1:nrow(h3k27acOverlap_1))
ggpubr::ggscatter(h3k27acOverlap_1,x='geneNumber',y='meanRPM',label = 'external_gene_name',cor.coef = T) + theme_ameres(type = "barplot")



pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/transcriptionInitiation/H3K27Ac_ATACseq_TP1_TP2.pdf")
      ac_openChr = data.frame(H3K27Ac = h3k27acOverlap_1_rep1$meanRPM,ATACseq = ATACseqOverlap_1_rep1$meanRPM,geneName_1 = ATACseqOverlap_1_rep1$external_gene_name,geneName_2 = h3k27acOverlap_1_rep1$external_gene_name)
      p = ggpubr::ggscatter(ac_openChr,x='H3K27Ac',y = 'ATACseq',cor.coef = T, xlab = 'H3K27Ac (4cell stage) RPM',ylab = 'ATACseq (64 cell stage) RPM',title = '2 cell Stage') + theme_ameres(type = "barplot")+ 
        geom_text_repel(data=subset(ac_openChr,  H3K27Ac > 0 | ATACseq > 0),aes(label = geneName_1)) 
      
      ac_openChr_TP2_h3k27 = data.frame(H3K27Ac = h3k27acOverlap_1$meanRPM,geneName_2 = h3k27acOverlap_1$external_gene_name)
      ac_openChr_TP2_ATACseq = data.frame(ATACseq = ATACseqOverlap_1$meanRPM,geneName_2 = ATACseqOverlap_1$external_gene_name)
      ac_openChr_TP2 = plyr::join(ac_openChr_TP2_h3k27,ac_openChr_TP2_ATACseq)
      q = ggpubr::ggscatter(ac_openChr_TP2,x='H3K27Ac',y = 'ATACseq',cor.coef = T , xlab = 'H3K27Ac (4 cell stage) RPM',ylab = 'ATACseq (64 cell stage) RPM',title = '64 cell Stage') + theme_ameres(type = "barplot") + 
        geom_text_repel(data=subset(ac_openChr_TP2,  H3K27Ac > 0 | ATACseq > 0),aes(label = geneName_2)) 
      print(p)
      print(q)
dev.off()

############### seem like quite a number of genes that are detected by SLAMseq also have H3K27Ac signal and habe open chromatin. ##############################

### combine three replicates for this analysis below.... to get the value of the fraction. 

Continuedexpression_TCtiming_Z_TP1 = Continuedexpression_TCtiming_Z %>% dplyr::filter(timepoint == "TP1")
Continuedexpression_TCtiming_Z_TP2 = Continuedexpression_TCtiming_Z %>% dplyr::filter(timepoint == "TP2")

Continuedexpression_TCtiming_Z_TP1_melt = melt(Continuedexpression_TCtiming_Z_TP1 %>% dplyr::select(dplyr::contains("conver")))
Continuedexpression_TCtiming_Z_TP1_melt$names = Continuedexpression_TCtiming_Z_TP1$Name
Continuedexpression_TCtiming_Z_TP1_melt$time = rep(c(0.75,2,2.5,3,3.5,4,4.5,5,5.5),each=nrow(Continuedexpression_TCtiming_Z_TP1))

ggpubr::ggline(Continuedexpression_TCtiming_Z_TP1_melt,x='time',y='value',group = 'names',color = 'names',facet.by = 'names')

## TP2
Continuedexpression_TCtiming_Z_TP2 = Continuedexpression_TCtiming_Z_TP2[-which(Continuedexpression_TCtiming_Z_TP2$Name == "ENSDART00000160631"),]
Continuedexpression_TCtiming_Z_TP2_melt = melt(Continuedexpression_TCtiming_Z_TP2 %>% dplyr::select(dplyr::contains("conver")))
Continuedexpression_TCtiming_Z_TP2_melt$names = Continuedexpression_TCtiming_Z_TP2$Name
Continuedexpression_TCtiming_Z_TP2_melt$time = rep(c(0.75,2,2.5,3,3.5,4,4.5,5,5.5),each=nrow(Continuedexpression_TCtiming_Z_TP2))

ggpubr::ggline(Continuedexpression_TCtiming_Z_TP2_melt,x='time',y='value',group = 'names',color = 'names',facet.by = 'names')


##### creating a data table for analysis of TC read counts.... - this table will contain - number of reads, number of TC reads, number of TC reads with 2 conversions, number of Ts converted, numbeer of Ts covered

numberOfreads = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/dataTables/reads_allCountingWindows_allSNPsRemoved_slamdunk3.4.txt",sep="\t",stringsAsFactors = F,header=T)
numberOfTCreads = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/dataTables/TCreads_allCountingWindows_allSNPsRemoved_slamdunk3.4.txt",sep="\t",stringsAsFactors = F,header=T)
numberOfTCreads_with2TC = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/dataTables/TCreads_allCountingWindows_2tc.txt",sep="\t",stringsAsFactors = F,header=T)
TsCovered = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/dataTables/TsCovered_allCountingWindows_allSNPsRemoved_slamdunk3.4.txt",sep="\t",stringsAsFactors = F,header=T)
TsConverted = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/dataTables/TsConverted_allCountingWindows_allSNPsRemoved_slamdunk3.4.txt",sep="\t",stringsAsFactors = F,header=T)


#### combining the replicates for each of these... 

splitReplicates_combine = function(dataFrameToSplit,condition,metadata_add){
  dataFrameToSplit_condition = dataFrameToSplit %>% select_if(grepl(condition,names(.)))
  dataFrameToSplit_condition_R1 = dataFrameToSplit_condition %>% select_if(grepl("R1",names(.)))
  dataFrameToSplit_condition_R2 = dataFrameToSplit_condition %>% select_if(grepl("R2",names(.)))
  dataFrameToSplit_condition_R3 = dataFrameToSplit_condition %>% select_if(grepl("R3",names(.)))
  mean_repl = (dataFrameToSplit_condition_R1+dataFrameToSplit_condition_R2+dataFrameToSplit_condition_R3)
  dataFrameToSplit_condition_R1 = cbind.data.frame(dataFrameToSplit_condition_R1,metadata_add)
  dataFrameToSplit_condition_R2 = cbind.data.frame(dataFrameToSplit_condition_R2,metadata_add)
  dataFrameToSplit_condition_R3 = cbind.data.frame(dataFrameToSplit_condition_R3,metadata_add)
  mean_repl = cbind.data.frame(mean_repl,metadata_add)
  splitReplicates = list(dataFrameToSplit_condition_R1,dataFrameToSplit_condition_R2,dataFrameToSplit_condition_R3,mean_repl)
  names(splitReplicates) = c("R1","R2","R3","mean")
  return(splitReplicates)
}

numberOfreads = splitReplicates_combine(dataFrameToSplit = numberOfreads,condition = "Inj",metadata_add = numberOfreads[,c(1:6)])
numberOfTCreads = splitReplicates_combine(dataFrameToSplit = numberOfTCreads,condition = "Inj",metadata_add = numberOfTCreads[,c(1:6)])
numberOfTCreads_with2TC = splitReplicates_combine(dataFrameToSplit = numberOfTCreads_with2TC,condition = "Inj",metadata_add = numberOfTCreads_with2TC[,c(1:6)])
TsCovered = splitReplicates_combine(dataFrameToSplit = TsCovered,condition = "Inj",metadata_add = TsCovered[,c(1:6)])
TsConverted = splitReplicates_combine(dataFrameToSplit = TsConverted,condition = "Inj",metadata_add = TsConverted[,c(1:6)])

numberOfreads = numberOfreads$mean
numberOfTCreads = numberOfTCreads$mean
numberOfTCreads_with2TC = numberOfTCreads_with2TC$mean
TsCovered = TsCovered$mean
TsConverted = TsConverted$mean

allDfsTOgether = list(numberOfreads,numberOfTCreads,numberOfTCreads_with2TC,TsCovered,TsConverted)
names(allDfsTOgether) = c("numberOfReads","numberOfTCreads","numberOf2TCreads","TsCovered","TsConverted")
TP1Data = cbind.data.frame(do.call(cbind.data.frame,lapply(allDfsTOgether,function(x) x[,1])),numberOfTCreads[,c("Chromosome","Start","End","Name","Length","Strand")]) %>% dplyr::mutate(conversionScore = (TsConverted)/ (TsCovered * numberOfReads )) %>%
  dplyr::mutate(TCweighting = numberOf2TCreads/numberOfTCreads)%>% dplyr::mutate(finalScore = conversionScore * TCweighting) %>% dplyr::mutate(label_add = paste0("co=",TsConverted,",","cv=",TsCovered,",","nRead=",numberOfReads,",","2tc=",numberOf2TCreads,",","tcReads=",numberOfTCreads))
TP2Data = cbind.data.frame(do.call(cbind.data.frame,lapply(allDfsTOgether,function(x) x[,2])),numberOfTCreads[,c("Chromosome","Start","End","Name","Length","Strand")]) %>% dplyr::mutate(conversionScore = (TsConverted )/ (TsCovered * numberOfReads )) %>% 
  dplyr::mutate(TCweighting = numberOf2TCreads/numberOfTCreads)%>% dplyr::mutate(finalScore = conversionScore * TCweighting) %>%  dplyr::mutate(label_add = paste0("co=",TsConverted,",","cv=",TsCovered,",","nRead=",numberOfReads,",","2tc=",numberOf2TCreads,",","tcReads=",numberOfTCreads))
TP3Data = cbind.data.frame(do.call(cbind.data.frame,lapply(allDfsTOgether,function(x) x[,3])),numberOfTCreads[,c("Chromosome","Start","End","Name","Length","Strand")]) %>% dplyr::mutate(conversionScore = (TsConverted )/ (TsCovered * numberOfReads )) %>% 
  dplyr::mutate(TCweighting = numberOf2TCreads/numberOfTCreads)%>% dplyr::mutate(finalScore = conversionScore * TCweighting) %>%  dplyr::mutate(label_add = paste0("co=",TsConverted,",","cv=",TsCovered,",","nRead=",numberOfReads,",","2tc=",numberOf2TCreads,",","tcReads=",numberOfTCreads))
TP4Data = cbind.data.frame(do.call(cbind.data.frame,lapply(allDfsTOgether,function(x) x[,4])),numberOfTCreads[,c("Chromosome","Start","End","Name","Length","Strand")]) %>% dplyr::mutate(conversionScore = (TsConverted)/ (TsCovered * numberOfReads )) %>% 
  dplyr::mutate(TCweighting = numberOf2TCreads/numberOfTCreads)%>% dplyr::mutate(finalScore = conversionScore * TCweighting) %>%  dplyr::mutate(label_add = paste0("co=",TsConverted,",","cv=",TsCovered,",","nRead=",numberOfReads,",","2tc=",numberOf2TCreads,",","tcReads=",numberOfTCreads))
TP5Data = cbind.data.frame(do.call(cbind.data.frame,lapply(allDfsTOgether,function(x) x[,5])),numberOfTCreads[,c("Chromosome","Start","End","Name","Length","Strand")]) %>% dplyr::mutate(conversionScore = (TsConverted )/ (TsCovered * numberOfReads )) %>% 
  dplyr::mutate(TCweighting = numberOf2TCreads/numberOfTCreads)%>% dplyr::mutate(finalScore = conversionScore * TCweighting) %>%  dplyr::mutate(label_add = paste0("co=",TsConverted,",","cv=",TsCovered,",","nRead=",numberOfReads,",","2tc=",numberOf2TCreads,",","tcReads=",numberOfTCreads))
TP6Data = cbind.data.frame(do.call(cbind.data.frame,lapply(allDfsTOgether,function(x) x[,6])),numberOfTCreads[,c("Chromosome","Start","End","Name","Length","Strand")]) %>% dplyr::mutate(conversionScore = (TsConverted )/ (TsCovered * numberOfReads )) %>% 
  dplyr::mutate(TCweighting = numberOf2TCreads/numberOfTCreads)%>% dplyr::mutate(finalScore = conversionScore * TCweighting) %>%  dplyr::mutate(label_add = paste0("co=",TsConverted,",","cv=",TsCovered,",","nRead=",numberOfReads,",","2tc=",numberOf2TCreads,",","tcReads=",numberOfTCreads))
TP7Data = cbind.data.frame(do.call(cbind.data.frame,lapply(allDfsTOgether,function(x) x[,7])),numberOfTCreads[,c("Chromosome","Start","End","Name","Length","Strand")]) %>% dplyr::mutate(conversionScore = (TsConverted)/ (TsCovered * numberOfReads )) %>%
  dplyr::mutate(TCweighting = numberOf2TCreads/numberOfTCreads)%>% dplyr::mutate(finalScore = conversionScore * TCweighting) %>%   dplyr::mutate(label_add = paste0("co=",TsConverted,",","cv=",TsCovered,",","nRead=",numberOfReads,",","2tc=",numberOf2TCreads,",","tcReads=",numberOfTCreads))
TP8Data = cbind.data.frame(do.call(cbind.data.frame,lapply(allDfsTOgether,function(x) x[,8])),numberOfTCreads[,c("Chromosome","Start","End","Name","Length","Strand")]) %>% dplyr::mutate(conversionScore = (TsConverted )/ (TsCovered * numberOfReads )) %>%
  dplyr::mutate(TCweighting = numberOf2TCreads/numberOfTCreads)%>% dplyr::mutate(finalScore = conversionScore * TCweighting) %>% dplyr::mutate(label_add = paste0("co=",TsConverted,",","cv=",TsCovered,",","nRead=",numberOfReads,",","2tc=",numberOf2TCreads,",","tcReads=",numberOfTCreads))
TP9Data = cbind.data.frame(do.call(cbind.data.frame,lapply(allDfsTOgether,function(x) x[,9])),numberOfTCreads[,c("Chromosome","Start","End","Name","Length","Strand")]) %>% dplyr::mutate(conversionScore = (TsConverted )/ (TsCovered * numberOfReads )) %>% 
  dplyr::mutate(TCweighting = numberOf2TCreads/numberOfTCreads)%>% dplyr::mutate(finalScore = conversionScore * TCweighting) %>% dplyr::mutate(label_add = paste0("co=",TsConverted,",","cv=",TsCovered,",","nRead=",numberOfReads,",","2tc=",numberOf2TCreads,",","tcReads=",numberOfTCreads))

###################### 
allData = list(TP1Data,TP2Data,TP3Data,TP4Data,TP5Data,TP6Data,TP7Data,TP8Data,TP9Data)


#### i want to plot these parameters for purely zygotic traqnscripts defined by RNAseq using :
  ### treated samples
  ### untreated samples

pureZtranscripts = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/microsymosium2018/data/zygoticGenes.txt",sep="\t")



TP1Data_subset = lapply(allData,function(x) plyr::join(Continuedexpression_TCtiming_Z_TP1,as.data.frame(x)))
TP2Data_subset = lapply(allData,function(x) plyr::join(Continuedexpression_TCtiming_Z_TP2,as.data.frame(x)))

TP1Data_subset = lapply(TP1Data_subset,function(x) x[,c("Name","label_add")])
TP2Data_subset = lapply(TP2Data_subset,function(x) x[,c("Name","label_add")])
names(TP1Data_subset) = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)
names(TP2Data_subset) =  c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)

TP1Data_subset_melt = melt(TP1Data_subset) %>% mutate(names = Name,time=L1)
TP2Data_subset_melt = melt(TP2Data_subset) %>% mutate(names = Name,time=L1)

Continuedexpression_TCtiming_Z_TP2_melt = plyr::join(Continuedexpression_TCtiming_Z_TP2_melt,TP2Data_subset_melt)

ggpubr::ggline(Continuedexpression_TCtiming_Z_TP2_melt,x='time',y='value',group = 'names',color = 'names',facet.by = 'names',label = "label_add",repel=T,ylim = c(0,0.25))
