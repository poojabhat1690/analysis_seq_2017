###### making TPMs and count data from SLAMseq samples ... 

##### first from ENSEMBL annotations.... 

countDataSets_quantSeq = list.files(path = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/initialEvaluationOfData/data/countData_6thFeb/",pattern = ".tsv")
countDataSets_quantSeq_paths = paste0( "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/initialEvaluationOfData/data/countData_6thFeb/",countDataSets_quantSeq)
# countDataSets_quantSeq = list.files(path = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/count//",pattern = ".tsv")
# countDataSets_quantSeq_paths = paste0( "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific/count/",countDataSets_quantSeq)

countDataSets_quantSeq_data = lapply(countDataSets_quantSeq_paths,function(x) read.table(x,sep="\t",header=T))

countData_timeCourse_datasets = countDataSets_quantSeq_data
splitNames = unlist(lapply(strsplit(countDataSets_quantSeq,"_",T),function(x) x[2]))
barcodes = unlist(lapply(strsplit(splitNames,".",T),function(x) x[1]))
names(countData_timeCourse_datasets) = barcodes

### reading in the sample information file

sampleInfo = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt",sep="\t")
countData_timeCourse_datasets_names = as.data.frame(names(countData_timeCourse_datasets))
colnames(countData_timeCourse_datasets_names) = "V2"

countData_timeCourse_datasets_names = plyr::join(countData_timeCourse_datasets_names,sampleInfo)
names(countData_timeCourse_datasets) = countData_timeCourse_datasets_names$V3

countDataSets_quantSeq_data = countData_timeCourse_datasets


timepoints = paste0("TP",c(1:9))

reads_samples = do.call(cbind.data.frame,lapply(countDataSets_quantSeq_data,function(x) x$ReadCount))
orderAll = (c(c(paste0("Inj_","R1_",timepoints),paste0("Inj_","R2_",timepoints),paste0("Inj_","R3_",timepoints)),c(paste0("Inc_","R1_",timepoints),paste0("Inc_","R2_",timepoints),paste0("Inc_","R3_",timepoints)),c(paste0("Untreated_",timepoints))))
reads_samples = reads_samples[,orderAll]
reads_samples$names = countDataSets_quantSeq_data$Inj_R1_TP8$Name
write.table(reads_samples,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/microsymosium2018/data/countData.txt",sep="\t",quote = F,row.names = F)

TCreads_samples = do.call(cbind.data.frame,lapply(countDataSets_quantSeq_data,function(x) x$TcReadCount))
orderAll = (c(c(paste0("Inj_","R1_",timepoints),paste0("Inj_","R2_",timepoints),paste0("Inj_","R3_",timepoints)),c(paste0("Inc_","R1_",timepoints),paste0("Inc_","R2_",timepoints),paste0("Inc_","R3_",timepoints)),c(paste0("Untreated_",timepoints))))
TCreads_samples = TCreads_samples[,orderAll]
TCreads_samples$names = countDataSets_quantSeq_data$Inj_R1_TP8$Name
write.table(TCreads_samples,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/microsymosium2018/data/TCcountData.txt",sep="\t",quote = F,row.names = F)


cpm_samples = do.call(cbind.data.frame,lapply(countDataSets_quantSeq_data,function(x) x$ReadsCPM))
orderAll = (c(c(paste0("Inj_","R1_",timepoints),paste0("Inj_","R2_",timepoints),paste0("Inj_","R3_",timepoints)),c(paste0("Inc_","R1_",timepoints),paste0("Inc_","R2_",timepoints),paste0("Inc_","R3_",timepoints)),c(paste0("Untreated_",timepoints))))
cpm_samples = cpm_samples[,orderAll]
cpm_samples$names = countDataSets_quantSeq_data$Inj_R1_TP8$Name
write.table(cpm_samples,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/microsymosium2018/data/RPMData.txt",sep="\t",quote = F,row.names = F)

#### cpms of T>C reads 

cpm_TCsamples = do.call(cbind.data.frame,lapply(countDataSets_quantSeq_data,function(x) x$ConversionRate))
orderAll = (c(c(paste0("Inj_","R1_",timepoints),paste0("Inj_","R2_",timepoints),paste0("Inj_","R3_",timepoints)),c(paste0("Inc_","R1_",timepoints),paste0("Inc_","R2_",timepoints),paste0("Inc_","R3_",timepoints)),c(paste0("Untreated_",timepoints))))
cpm_TCsamples = cpm_TCsamples[,orderAll]
cpm_TCsamples$names = countDataSets_quantSeq_data$Inj_R1_TP8$Name
write.table(cpm_TCsamples,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/microsymosium2018/data/TCConversionRate.txt",sep="\t",quote = F,row.names = F)



######################################################################################################

######### from custom annotation 

######################################################################################################




###### making TPMs and count data from SLAMseq samples ... 
countDataSets_quantSeq = list.files(path = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11//count/",pattern = ".tsv")
countDataSets_quantSeq_paths = paste0( "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11//count/",countDataSets_quantSeq)
countDataSets_quantSeq_data = lapply(countDataSets_quantSeq_paths,function(x) read.table(x,sep="\t",header=T))

countData_timeCourse_datasets = countDataSets_quantSeq_data
splitNames = unlist(lapply(strsplit(countDataSets_quantSeq,"_",T),function(x) x[2]))
barcodes = unlist(lapply(strsplit(splitNames,".",T),function(x) x[1]))
names(countData_timeCourse_datasets) = barcodes

### reading in the sample information file

sampleInfo = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt",sep="\t")
countData_timeCourse_datasets_names = as.data.frame(names(countData_timeCourse_datasets))
colnames(countData_timeCourse_datasets_names) = "V2"

countData_timeCourse_datasets_names = plyr::join(countData_timeCourse_datasets_names,sampleInfo)
names(countData_timeCourse_datasets) = countData_timeCourse_datasets_names$V3

countDataSets_quantSeq_data = countData_timeCourse_datasets


timepoints = paste0("TP",c(1:9))

reads_samples = do.call(cbind.data.frame,lapply(countDataSets_quantSeq_data,function(x) x$ReadCount))
orderAll = (c(c(paste0("Inj_","R1_",timepoints),paste0("Inj_","R2_",timepoints),paste0("Inj_","R3_",timepoints)),c(paste0("Inc_","R1_",timepoints),paste0("Inc_","R2_",timepoints),paste0("Inc_","R3_",timepoints)),c(paste0("Untreated_",timepoints))))
reads_samples = reads_samples[,orderAll]
reads_samples$names = countDataSets_quantSeq_data$Inj_R1_TP8$Name
write.table(reads_samples,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/microsymosium2018/data/countData_customAnnotation.txt",sep="\t",quote = F,row.names = F)


TCreads_samples = do.call(cbind.data.frame,lapply(countDataSets_quantSeq_data,function(x) x$TcReadCount))
orderAll = (c(c(paste0("Inj_","R1_",timepoints),paste0("Inj_","R2_",timepoints),paste0("Inj_","R3_",timepoints)),c(paste0("Inc_","R1_",timepoints),paste0("Inc_","R2_",timepoints),paste0("Inc_","R3_",timepoints)),c(paste0("Untreated_",timepoints))))
TCreads_samples = TCreads_samples[,orderAll]
TCreads_samples$names = countDataSets_quantSeq_data$Inj_R1_TP8$Name
write.table(TCreads_samples,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/microsymosium2018/data/TCcountData_customAnnotation.txt",sep="\t",quote = F,row.names = F)


cpm_samples = do.call(cbind.data.frame,lapply(countDataSets_quantSeq_data,function(x) x$ReadsCPM))
orderAll = (c(c(paste0("Inj_","R1_",timepoints),paste0("Inj_","R2_",timepoints),paste0("Inj_","R3_",timepoints)),c(paste0("Inc_","R1_",timepoints),paste0("Inc_","R2_",timepoints),paste0("Inc_","R3_",timepoints)),c(paste0("Untreated_",timepoints))))
cpm_samples = cpm_samples[,orderAll]
cpm_samples$names = countDataSets_quantSeq_data$Inj_R1_TP8$Name
write.table(cpm_samples,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/microsymosium2018/data/RPMData_customAnnotation.txt",sep="\t",quote = F,row.names = F)

#### cpms of T>C reads 

cpm_TCsamples = do.call(cbind.data.frame,lapply(countDataSets_quantSeq_data,function(x) x$ConversionRate))
orderAll = (c(c(paste0("Inj_","R1_",timepoints),paste0("Inj_","R2_",timepoints),paste0("Inj_","R3_",timepoints)),c(paste0("Inc_","R1_",timepoints),paste0("Inc_","R2_",timepoints),paste0("Inc_","R3_",timepoints)),c(paste0("Untreated_",timepoints))))
cpm_TCsamples = cpm_TCsamples[,orderAll]
cpm_TCsamples$names = countDataSets_quantSeq_data$Inj_R1_TP8$Name
write.table(cpm_TCsamples,"/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/microsymosium2018/data/TCConversionRate_customAnnotation.txt",sep="\t",quote = F,row.names = F)










