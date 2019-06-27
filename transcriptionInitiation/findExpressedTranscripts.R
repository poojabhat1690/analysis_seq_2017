######## identification of transcripts that are expressed before the main wave of ZGA.. 
#### the idea is to use 3 parameters to find these transcripts... 
#### should be greater than the background conversion rate
#### should have percentage TC conversion rate greater than a set of purely zygotic genes. 
#### should have multiple TC conversions..

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
splitReplicates = function(dataFrameToSplit,condition,metadata_add){
  dataFrameToSplit_condition = dataFrameToSplit[,grep(condition,colnames(dataFrameToSplit ))]
  # dataFrameToSplit_condition_R1 = dataFrameToSplit_condition[,grep("R1",colnames(dataFrameToSplit_condition ))]
  dataFrameToSplit_condition_R2 = dataFrameToSplit_condition[,grep("R2",colnames(dataFrameToSplit_condition ))]
  dataFrameToSplit_condition_R3 = dataFrameToSplit_condition[,grep("R3",colnames(dataFrameToSplit_condition ))]
  mean_repl = (dataFrameToSplit_condition_R2+dataFrameToSplit_condition_R3)/2
  #dataFrameToSplit_condition_R1 = cbind.data.frame(dataFrameToSplit_condition_R1,metadata_add)
  dataFrameToSplit_condition_R2 = cbind.data.frame(dataFrameToSplit_condition_R2,metadata_add)
  dataFrameToSplit_condition_R3 = cbind.data.frame(dataFrameToSplit_condition_R3,metadata_add)
  mean_repl = cbind.data.frame(mean_repl,metadata_add)
  splitReplicates = list(dataFrameToSplit_condition_R2,dataFrameToSplit_condition_R3,mean_repl)
  names(splitReplicates) = c("R2","R3","mean")
  return(splitReplicates)
}

##### finding genes that have conversions at the first two time points and subtracting stage specific background from them

classifiedTranscripts = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/dataTables_expandedCounting/countingWindows_classified_conversions.txt",sep="\t",stringsAsFactors = F, header = T)
classifiedTranscripts = classifiedTranscripts %>% dplyr::filter(class == "Z" | class == "MZ") %>% dplyr::select(gene)


#### i need to get the background conversions per replicate and subtract from the conversion rate. 
errorRates = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019//errorRates_predicted_observed.txt",header = T)

conversionRates_all = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017//analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019//dataTables_expandedCounting//perGenes_conversion.txt",sep="\t",stringsAsFactors = F,header=T)

conversionRates_all[is.na(conversionRates_all)] <- 0
conversionRates_all = as.data.frame(conversionRates_all)
conversionRates_all_injection = conversionRates_all %>% dplyr::select(-starts_with("Inc"))
conversionRates_all_injection = conversionRates_all_injection %>% dplyr::select(starts_with("Inj"))

for(i in 1:ncol(conversionRates_all_injection)){
  conversionRates_all_injection[,i] = conversionRates_all_injection[,i] - errorRates$predicted[i] 
}
conversionRates_all_injection = as.matrix(conversionRates_all_injection)
conversionRates_all_injection[which(conversionRates_all_injection <0 )] <- 0
conversionRates_all_injection[grep("Inf",conversionRates_all_injection)]<-0
conversionRates_all_injection = as.data.frame(conversionRates_all_injection)
conversionRates_all_injection$name = conversionRates_all$name

##### splitting into replicates                    
conversionRates_all_injection_split = splitReplicates(dataFrameToSplit = conversionRates_all_injection,condition = "Inj",
                                                      metadata_add = conversionRates_all_injection$name)
####### now getting the MZ and Z genes from these and checking if these are expressed at stages 1 and 2.                    
conversionRates_all_injection_split = lapply(conversionRates_all_injection_split,function(x) 
  x[x$metadata_add %in% classifiedTranscripts$gene,])           

conversionRates_all_injection_split = lapply(conversionRates_all_injection_split,function(x)   x %>%
                                               dplyr::mutate(hpr_0.75 = if_else(.[,1]>0,true = T,false = F)) %>% 
                                               dplyr::mutate(hpr_0.2 = if_else(.[,2]>0,true = T,false = F)) )


######## now i want to get the percentage TC conversions... i.e (#TC converted * 100 * totalNumberOf reads)/ (#TC convered)
##### i.e conversion rate * total number of reads. 

#### reading in the total number of reads... 

readCounts = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019//dataTables_expandedCounting/totalCounts_allCws.txt",
                        stringsAsFactors = F,sep="\t",header = T)
readCounts = readCounts %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE) ### per gene read count
readCounts = readCounts[readCounts$name %in% conversionRates_all_injection_split$R2$metadata_add,]
### now splitting the read counts per replicate. 
readCounts_split = splitReplicates(dataFrameToSplit = readCounts,condition = "Inj",
                                   metadata_add = readCounts$name)


#### creating the perentage TC score by multiplying this conversion TC with the read counts. 
percentageTCscore = conversionRates_all_injection_split
for(i in 1:length(readCounts_split)){
  percentageTCscore[[i]][,c(1:9)] = (conversionRates_all_injection_split[[i]][,c(1:9)] * 100) / readCounts_split[[i]][,c(1:9)]
}  

########### getting the multiple TC score... 

readsWithMultipleTC = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019//dataTables_expandedCounting/numberOfreadsWithMultipleTC.txt",
                                 stringsAsFactors = F, sep="\t",header = T)                      


readsWithMultipleTC = readsWithMultipleTC %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE) ### per gene read count
readsWithMultipleTC = readsWithMultipleTC[readsWithMultipleTC$name %in% conversionRates_all_injection_split$R2$metadata_add,]
readsWithMultipleTC_split = splitReplicates(dataFrameToSplit = readsWithMultipleTC,condition = "Inj",
                                            metadata_add = readsWithMultipleTC$name)




### reads with single TC

readsWithTC = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019//dataTables_expandedCounting/numberOfreadsWithTC.txt",
                         stringsAsFactors = F, sep="\t",header = T)                      


readsWithTC = readsWithTC %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE) ### per gene read count
readsWithTC = readsWithTC[readsWithTC$name %in% conversionRates_all_injection_split$R2$metadata_add,]
readsWithTC_split = splitReplicates(dataFrameToSplit = readsWithTC,condition = "Inj",
                                    metadata_add = readsWithTC$name)

### fraction of multipleTC reads... 


#### i also want to know the number of Ts covered ... i .e the number of Ts in each counting window..

slamdunkOp = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/count/combinedFile_AACGCC.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_tcount.tsv",sep="\t",stringsAsFactors = F,header = T)
slamdunkOp = slamdunkOp %>% dplyr::group_by(Name)  %>%   dplyr::summarise_at(vars(Tcontent), sum, na.rm = TRUE) ### per gene read count
slamdunkOp = slamdunkOp[slamdunkOp$Name %in% conversionRates_all_injection_split$R2$metadata_add,]

fractionMultipleTC = readsWithTC_split
for(i in 1:3){
  fractionMultipleTC[[i]][,c(1:9)] = (readsWithMultipleTC_split[[i]][,c(1:9)]*100)/ (readsWithTC_split[[i]][,c(1:9)] * slamdunkOp$Tcontent)
  
}


##### getting these parameters for purely zygotic transcripts - defined previously by RNAseq. 
zygoticGenes = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/microsymosium2018/data/zygoticGenes.txt",sep="\t",stringsAsFactors = F)
zygoticGeneMultipleTC = lapply(fractionMultipleTC,function(x) x[x$metadata_add %in% zygoticGenes$external_gene_name,] )
zygoticGenepercentageTCscore = lapply(percentageTCscore,function(x) x[x$metadata_add %in% zygoticGenes$external_gene_name,] )
zygoticGeneMultipleTC_conf = unlist(lapply(zygoticGeneMultipleTC,function(x) t.test(x[,9])$conf.int[1])) ##3 getting the confidence interval
zygoticGenepercentageTCscore_conf = unlist(lapply(zygoticGenepercentageTCscore,function(x) t.test(x[,9])$conf.int[1])) ##3 getting the confidence interval

### comparing multiple TC with the considence interval...

percentageTCscore$R2  =  percentageTCscore$R2 %>% dplyr::mutate(precentageTC_TP1 = if_else(Inj_R2_TP1 > zygoticGenepercentageTCscore_conf[1],T,F)) %>%
  dplyr::mutate(precentageTC_TP2 = if_else(Inj_R2_TP2 > zygoticGenepercentageTCscore_conf[1],T,F))
percentageTCscore$R3  =  percentageTCscore$R3 %>% dplyr::mutate(precentageTC_TP1 = if_else(Inj_R3_TP2 > zygoticGenepercentageTCscore_conf[2],T,F)) %>%
  dplyr::mutate(precentageTC_TP2 = if_else(Inj_R3_TP2 > zygoticGenepercentageTCscore_conf[2],T,F))

fractionMultipleTC$R2  =  fractionMultipleTC$R2 %>% dplyr::mutate(MultipleTC_TP1 = if_else(Inj_R2_TP1 > zygoticGeneMultipleTC_conf[1],T,F)) %>%
  dplyr::mutate(MultipleTC_TP2 = if_else(Inj_R2_TP2 > zygoticGeneMultipleTC_conf[1],T,F))
fractionMultipleTC$R3  =  fractionMultipleTC$R3 %>% dplyr::mutate(MultipleTC_TP1 = if_else(Inj_R3_TP2 > zygoticGeneMultipleTC_conf[2],T,F)) %>%
  dplyr::mutate(MultipleTC_TP2 = if_else(Inj_R3_TP2 > zygoticGeneMultipleTC_conf[2],T,F))

total_R2 = cbind.data.frame(fractionMultipleTC$R2 %>% dplyr::select(c(MultipleTC_TP1,MultipleTC_TP2)), percentageTCscore$R2 %>% 
                              dplyr::select(c(precentageTC_TP1,precentageTC_TP2,hpr_0.75,hpr_0.2, metadata_add)))
R2_TP1 = total_R2 %>% dplyr::filter(hpr_0.75 == T & (precentageTC_TP1 == T |  MultipleTC_TP1 == T)) %>% dplyr::mutate(replicate = "R2")  %>%
  dplyr::mutate(time = 0.75)
R2_TP2 = total_R2 %>% dplyr::filter(hpr_0.2 == T & (precentageTC_TP2 == T |  MultipleTC_TP2 == T)) %>% dplyr::mutate(replicate = "R2")  %>%
  dplyr::mutate(time = 2)


total_R3 = cbind.data.frame(fractionMultipleTC$R3 %>% dplyr::select(c(MultipleTC_TP1,MultipleTC_TP2)), percentageTCscore$R3 %>% 
                              dplyr::select(c(precentageTC_TP1,precentageTC_TP2,hpr_0.75,hpr_0.2, metadata_add)))
R3_TP1 = total_R3 %>% dplyr::filter(hpr_0.75 == T & (precentageTC_TP1 == T |  MultipleTC_TP1 == T)) %>% dplyr::mutate(replicate = "R3")  %>%
  dplyr::mutate(time = 0.75)
R3_TP2 = total_R3 %>% dplyr::filter(hpr_0.2 == T & (precentageTC_TP2 == T |  MultipleTC_TP2 == T)) %>% dplyr::mutate(replicate = "R3")  %>%
  dplyr::mutate(time = 2)

#### total time points... 

R3_combined = rbind.data.frame(R3_TP2,R3_TP1)
R2_combined = rbind.data.frame(R2_TP2,R2_TP1)
totalGenes = rbind.data.frame(R2_combined,R3_combined)
#totalGenes = totalGenes[!duplicated(totalGenes$metadata_add),]
dir.create("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/transcriptionInitiation/")
write.table(totalGenes,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/transcriptionInitiation/genes_T1_T2.txt",
            quote = F,row.names = F,sep="\t")


###### plotting this data... 
classifiedTranscripts = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/dataTables_expandedCounting/countingWindows_classified_conversions.txt",sep="\t",stringsAsFactors = F, header = T)

classifiedTranscripts_intotalGenes = classifiedTranscripts[classifiedTranscripts$gene %in% totalGenes$metadata_add,]

totalGenes$metadata_add = as.character(totalGenes$metadata_add)

TP1 = totalGenes %>% dplyr::filter(time == 0.75)  %>% dplyr::mutate(gene = metadata_add) %>% plyr::join(.,classifiedTranscripts_intotalGenes) %>% dplyr::distinct(gene,.keep_all = TRUE)
TP2 = totalGenes %>% dplyr::filter(time == 2) %>% dplyr::mutate(gene = metadata_add) %>% plyr::join(.,classifiedTranscripts_intotalGenes)%>% dplyr::distinct(gene,.keep_all = TRUE)

TP1_table = data.frame(table(TP1$class), time = 0.75) 
TP2_table = data.frame(table(TP2$class), time = 2) 

totalData = rbind.data.frame(TP1_table,TP2_table)

p = ggpubr::ggbarplot(totalData, x='time',y='Freq',col='Var1', fill='Var1',xlab = 'Time (hpf)', ylab = 'Number of genes',palette = 'Set1') + theme_ameres(type = 'barplot') 
p


#ggsave("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/transcriptionInitiation/genes_expressedTP1_TP2.pdf",p,height = 5,width=3)

conversions_samples = melt(classifiedTranscripts_intotalGenes %>% dplyr::select(c("Inj_R2_TP1","Inj_R2_TP2")))
p = ggpubr::ggboxplot(conversions_samples, x='variable',y='value',ylim = c(0,0.075),fill='red',ylab = 'TC conversion rate',xlab = 'Time (hpf)') + theme_ameres(type = "barplot")
p
ggsave("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/transcriptionInitiation/conversions_TP1_TP2.pdf",p,height = 5,width=3)



#### creating an alluvial plot... 
TP1 = TP1 %>% dplyr::select(c('class','time','metadata_add'))
TP2 = TP2 %>% dplyr::select(c('class','time','metadata_add'))
totalTPs = rbind.data.frame(TP1,TP2)
totalTPs$time = as.factor(totalTPs$time)

p = ggplot(totalTPs,
       aes(x = time, stratum = class, alluvium = metadata_add,
           fill = class, label = class)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft") +
  geom_stratum() +
  theme(legend.position = "bottom") +
  theme_ameres(type = 'barplot') + geom_label(stat = "stratum") + xlab('Study')
p
ggsave("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/transcriptionInitiation/genes_expressedTP1_TP2.pdf",p,height = 5,width=3)
