library(dplyr)
library(tidyr)
library(reshape)
library(ggplot2)

mutationRatesFiles = list.files(path = "/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11//utrRates/",pattern = ".csv")
mutationRatesFiles_path = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11//utrRates/",mutationRatesFiles)
mutationRates = lapply(mutationRatesFiles_path,function(x) read.table(x,stringsAsFactors = F,sep="\t",header = T))
names(mutationRates) = mutationRatesFiles
######## sample info 
zygoticGenes = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/microsymosium2018/data/zygoticGenes.txt",sep="\t",stringsAsFactors = F)

sampleInfo = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt",sep="\t",stringsAsFactors = F)
sampleInfo = sampleInfo[order(sampleInfo$V2),]
names(mutationRates) = sampleInfo$V3
mutationRates = lapply(mutationRates, function(x) x[x$Name %in% zygoticGenes$external_gene_name,])
Injection_R1 = mutationRates[grep("Inj_R1",names(mutationRates))]
Injection_R2 =  mutationRates[grep("Inj_R2",names(mutationRates))]
Injection_R3 =  mutationRates[grep("Inj_R3",names(mutationRates))]


Incubation_R1 = mutationRates[grep("Inc_R1",names(mutationRates))]
Incubation_R2 =  mutationRates[grep("Inc_R2",names(mutationRates))]
Incubation_R3 =  mutationRates[grep("Inc_R3",names(mutationRates))]

Untreated = mutationRates[grep("Untreated",names(mutationRates))]

getTCrates = function(dataFrame){
  
  
  
  T_C_all = vector("list",length(dataFrame))
  names(T_C_all) = names(dataFrame)
  
  
  
  
  for(i in 1:length(dataFrame))
  {
    
    
    curTab = dataFrame[[i]]
    #curTab =  curTab[which(curTab$ReadCount > 10),]
    plusTab = curTab %>% dplyr::filter(Strand == "+")
    minusTab = curTab %>% dplyr::filter(Strand == "-") %>%
      dplyr::select(Name, Chr, Start, End,Strand,ReadCount, A_A = T_T, G_G = C_C, C_C = G_G, T_T = A_A, A_C = T_G, A_G = T_C, A_T = T_A, C_A = G_T, C_G = G_C, C_T = G_A, G_A = C_T, G_C = C_G, G_T = C_A, T_A = A_T, T_C = A_G, T_G = A_C)

    
    plusTab = plusTab  %>%
      mutate(Asum = A_A + A_C + A_G + A_T) %>%
      mutate(Csum = C_A + C_C + C_G + C_T) %>% 
      mutate(Gsum = G_A + G_C + G_G + G_T) %>%
      mutate(Tsum = T_A + T_C + T_G + T_T) %>%
      mutate_each(funs(. / Asum),matches("A_")) %>%
      mutate_each(funs(. / Csum),matches("C_")) %>%
      mutate_each(funs(. / Gsum),matches("G_")) %>%
      mutate_each(funs(. / Tsum),matches("T_")) 
   plusTab = plusTab %>%  dplyr::select(Name, Chr, Start, End,Strand,ReadCount, A_A , G_G , C_C , T_T, A_C, A_G , A_T , C_A, C_G, C_T , G_A, G_C, G_T, T_A , T_C, T_G , Asum, Csum,Gsum,Tsum)

    
    minusTab = minusTab %>%
      mutate(Asum = A_A + A_C + A_G + A_T) %>%
      mutate(Csum = C_A + C_C + C_G + C_T) %>% 
      mutate(Gsum = G_A + G_C + G_G + G_T) %>%
      mutate(Tsum = T_A + T_C + T_G + T_T) %>%
      mutate_each(funs(. / Asum),matches("A_")) %>%
      mutate_each(funs(. / Csum),matches("C_")) %>%
      mutate_each(funs(. / Gsum),matches("G_")) %>%
      mutate_each(funs(. / Tsum),matches("T_")) 
    
  
        T_C_all[[i]] = rbind(plusTab, minusTab)
        T_C_all[[i]] = T_C_all[[i]] %>% dplyr::filter(ReadCount > 0)  %>% dplyr::mutate(FractionPerRead = T_C/ReadCount)
  }
  
  return(T_C_all)
}


#### plotting stuff
TCrates_injection_R1 = getTCrates(dataFrame = Injection_R1)
TCrates_injection_R2 = getTCrates(dataFrame = Injection_R2)
TCrates_injection_R3 = getTCrates(dataFrame = Injection_R3)

TCrates_incubation_R1 = getTCrates(dataFrame = Incubation_R1)
TCrates_incubation_R2 = getTCrates(dataFrame = Incubation_R2)
TCrates_incubation_R3 = getTCrates(dataFrame = Incubation_R3)

TCrates_Untreated = getTCrates(dataFrame = Untreated)

fractionPerRead_untreated = melt(lapply(TCrates_Untreated,function(x) mean(x$FractionPerRead,na.rm=T))) %>% dplyr::mutate(type = "Untreated")
fractionPerRead_Inj1 = melt(lapply(TCrates_injection_R1,function(x) mean(x$FractionPerRead,na.rm=T))) %>% dplyr::mutate(type = "Inj1")
fractionPerRead_Inj2 = melt(lapply(TCrates_injection_R2,function(x) mean(x$FractionPerRead,na.rm=T))) %>% dplyr::mutate(type = "Inj2")
fractionPerRead_Inj3 = melt(lapply(TCrates_injection_R3,function(x) mean(x$FractionPerRead,na.rm=T))) %>% dplyr::mutate(type = "Inj3")

allFractionPerRead = rbind(fractionPerRead_untreated,fractionPerRead_Inj1,fractionPerRead_Inj2,fractionPerRead_Inj3) 
allFractionPerRead = allFractionPerRead[order(allFractionPerRead$L1),]
allFractionPerRead$time = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)
ggplot(data = allFractionPerRead,aes(x=time,y=log10(value),group=type,col=type)) + geom_line(size=1)



