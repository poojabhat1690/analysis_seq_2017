##### to determine the number of TC reads to use.. I want calculate the fraction of reads with TC conversions in the background 
##### and compare this to the fraction of TC conversions in the treatment samples...

library(dplyr)
library(ggplot2)
library(scales)
theme_ameres <- function (type) {  ### plotting function
  
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


######### preparing the data tables...of the fraction of reads with 1,2,...n TCs... (run on cluster)

####### getting the fraction of reads separately for mitochondrial genes


folder = "/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/probabilityOfFindingTC/zygoticTranscripts_vectorized/"
treatmentSamples = list.files(path = folder,pattern = "freqAllMutations_allTranscripts_totalSNP*")
treatmentSamples_Inj = treatmentSamples[grep("Inj",treatmentSamples)] %>% paste0(folder,.)
treatmentSamples_Unt = treatmentSamples[grep("Unt",treatmentSamples)] %>% paste0(folder,.)

getTabulatedTcs = function(sampleName){ ###3 get fraction of reads with 1TC, 2 TCs... etc...
  sampleData = read.table(sampleName)
  minusTab = sampleData %>% dplyr::filter(strand == "-") %>%
    dplyr::select( A_C = T_G, A_G = T_C, A_T = T_A, C_A = G_T, C_G = G_C, C_T = G_A, G_A = C_T, G_C = C_G, G_T = C_A, T_A = A_T, T_C = A_G, T_G = A_C,id,strand)
  plusTab = sampleData %>% dplyr::filter(strand == "+") %>%
    dplyr::select( A_C, A_G, A_T, C_A , C_G , C_T , G_A , G_C , G_T , T_A , T_C, T_G ,id,strand)
  sampleData  = rbind(plusTab,minusTab)
  sampleData_mt = dplyr::filter(sampleData, grepl("chrM",id))
  sampleData = dplyr::filter(sampleData, !grepl("chrM",id))

  ### all excluding mitochrondria
  #fractionTC_samples = table(sampleData$T_C)/nrow(sampleData)
  sampleData = sampleData %>% dplyr::select(-matches('A_A|G_G|C_C|T_T'))
  sampleData_mt = sampleData_mt %>% dplyr::select(-matches('A_A|G_G|C_C|T_T'))

  fractionTC_samples = apply(sampleData[,1:12],2,function(x) table(x)/length(x))
  totalReads =   apply(sampleData[,1:12],2,function(x) table(x))
  fractionTC_samples_reads = vector("list",2)
  names(fractionTC_samples_reads) = c("fractionTC","numberOfreads")
  fractionTC_samples_reads[[1]] = fractionTC_samples
  fractionTC_samples_reads[[2]] = totalReads

  ##### only mitochondria
  #fractionTC_samplesMt = table(sampleData_mt$T_C)/nrow(sampleData_mt)
  fractionTC_samplesMt = apply(sampleData_mt[,1:12],2,function(x) table(x)/length(x))
  totalReads_mt =   apply(sampleData_mt[,1:12],2,function(x) table(x))

  #totalReads_mt = table(sampleData_mt$T_C)
  fractionTC_samples_reads_mt = vector("list",2)
  names(fractionTC_samples_reads_mt) = c("fractionTC","numberOfreads")
  fractionTC_samples_reads_mt[[1]] = fractionTC_samplesMt
  fractionTC_samples_reads_mt[[2]] = totalReads_mt

  totalFracs = list(fractionTC_samples_reads_mt,fractionTC_samples_reads)
  names(totalFracs) = c("mt","allOther")
  return(totalFracs)
}


sample_TCs = lapply(treatmentSamples_Unt,function(x) getTabulatedTcs(sampleName = x))
names(sample_TCs) = treatmentSamples_Unt
save(sample_TCs ,file = "//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/FractionOfTCreads_totalNumberOfreads_unInjectedSamples.Rdata")

sample_TCs_injection = lapply(treatmentSamples_Inj,function(x) getTabulatedTcs(sampleName = x))
names(sample_TCs_injection) = treatmentSamples_Inj
save(sample_TCs_injection ,file = "//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/FractionOfTCreads_totalNumberOfreads_InjectedSamples.Rdata")


####### plotting these datasets...

load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/FractionOfTCreads_totalNumberOfreads_InjectedSamples.Rdata")
names(sample_TCs_injection) = substr(names(sample_TCs_injection),start = 227,stop = 250)

load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/FractionOfTCreads_totalNumberOfreads_unInjectedSamples.Rdata")
names(sample_TCs) = substr(names(sample_TCs),start = 211,stop = 224)


      #### for each timepoint... 
library(data.table)
library(ggplot2)
library(scales)
# 
# pdf("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/errorRates_untreatedTCs/FractionOfReadsWithTCs_perTImepoint.pdf",height=4,width=5)
#     TPs = paste0("TP",c(1:9))
#   for(i in 1:length(TPs)){
#     load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/FractionOfTCreads_totalNumberOfreads_InjectedSamples.Rdata")
#     names(sample_TCs_injection) = substr(names(sample_TCs_injection),start = 211,stop = 224)
#     
#     load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/FractionOfTCreads_totalNumberOfreads_unInjectedSamples.Rdata")
#     names(sample_TCs) = substr(names(sample_TCs),start = 211,stop = 224)
#     
#     # p = lapply(sample_TCs_injection[grep(TPs[i],names(sample_TCs_injection))] ,function(x) as.data.frame(x[[1]])) %>% bind_rows(.id = 'Var1') %>% 
#     #   stats::setNames(c("sample", "Var1", "Freq")) %>%  mutate(Var1 = as.numeric(Var1)) %>%
#     #   ggpubr::ggline(data = .,x = 'Var1',y = 'Freq',group = 'sample',color = 'sample')  + ylab("FractionOfreads") + theme_ameres(type = "barplot") +
#     #     scale_color_brewer(palette = "Dark2") + xlab("Number of T_Cs per read") +  scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x ),labels = trans_format("log10", math_format(10^.x)))
#     # 
#     # print(p)
#     # 
#     #### combining the replicates and also plotting the untreated sample
#     sample_TCs_injection_total = sample_TCs_injection
#     sample_TCs_injection = lapply(sample_TCs_injection,function(x) x$allOther)
#     
#     
#     
#     numberOfTCs = data.frame(Var1 = unique(unlist(lapply(sample_TCs_injection[grep(TPs[i],names(sample_TCs_injection))],function(x) as.numeric(names(x[[2]]))))) )
#     numberOfReads = lapply(sample_TCs_injection[grep(TPs[i],names(sample_TCs_injection))],function(x) data.frame(x[[2]]))
#     numberOfReads$numberOfTCs = numberOfTCs
#   
#     sample_TCs_all = sample_TCs
#     sample_TCs =  lapply(sample_TCs,function(x) x$allOther)
#     numberOfTCs_uninjected = data.frame(Var1 = unique(unlist(lapply(sample_TCs[grep(TPs[i],names(sample_TCs))],function(x) as.numeric(names(x[[2]]))))) )
#     numberOfReads_uninjected = lapply(sample_TCs[grep(TPs[i],names(sample_TCs))],function(x) data.frame(x[[2]]))
#     numberOfReads$uninjected = numberOfReads_uninjected[[1]]
#     
#     q =  plyr::join_all(numberOfReads, by='Var1') %>% stats::setNames(c("nTC", "R1", "R2","R3","Uninjected")) %>%  replace(., is.na(.), 0) %>% 
#       mutate(nTotalReads = R1+R2+R3) %>% mutate(fractionTC = nTotalReads/sum(nTotalReads)) %>%  mutate(fractionTC_untreated = Uninjected/sum(Uninjected)) %>%
#       dplyr::select('nTC','fractionTC','fractionTC_untreated')  %>% tidyr::gather(nTC,value) %>% mutate(numTC = rep(1:(nrow(.)/2),2)-1) %>%
#       ggpubr::ggline(.,x = 'numTC',y='value',group = 'nTC',color = 'nTC')+ theme_ameres(type = "barplot") + ylab("Fraction of reads") + 
#       xlab("number of TCs") + scale_color_brewer(palette = "Set1") + ggtitle(TPs[i])+  scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x ),
#                                                                                                           labels = trans_format("log10", math_format(10^.x)))
#     print(q)
#     
#     
# 
#     sample_TCs_injection = lapply(sample_TCs_injection_total,function(x) x$mt)
#     numberOfTCs = data.frame(Var1 = unique(unlist(lapply(sample_TCs_injection[grep(TPs[i],names(sample_TCs_injection))],function(x) as.numeric(names(x[[2]]))))) )
#     numberOfReads = lapply(sample_TCs_injection[grep(TPs[i],names(sample_TCs_injection))],function(x) data.frame(x[[2]]))
#     numberOfReads$numberOfTCs = numberOfTCs
#     
#     sample_TCs =  lapply(sample_TCs_all,function(x) x$mt)
#     numberOfTCs_uninjected = data.frame(Var1 = unique(unlist(lapply(sample_TCs[grep(TPs[i],names(sample_TCs))],function(x) as.numeric(names(x[[2]]))))) )
#     numberOfReads_uninjected = lapply(sample_TCs[grep(TPs[i],names(sample_TCs))],function(x) data.frame(x[[2]]))
#     numberOfReads$uninjected = numberOfReads_uninjected[[1]]
#     
#     r =  plyr::join_all(numberOfReads, by='Var1') %>% stats::setNames(c("nTC", "R1", "R2","R3","Uninjected")) %>%  replace(., is.na(.), 0) %>% 
#       mutate(nTotalReads = R1+R2+R3) %>% mutate(fractionTC = nTotalReads/sum(nTotalReads)) %>%  mutate(fractionTC_untreated = Uninjected/sum(Uninjected)) %>%
#       dplyr::select('nTC','fractionTC','fractionTC_untreated')  %>% tidyr::gather(nTC,value) %>% mutate(numTC = rep(1:(nrow(.)/2),2)-1) %>%
#       ggpubr::ggline(.,x = 'numTC',y='value',group = 'nTC',color = 'nTC')+ theme_ameres(type = "barplot") + ylab("Fraction of reads") + 
#       xlab("number of TCs") + scale_color_brewer(palette = "Set1") + ggtitle(paste0(TPs[i],"-mt"))+  scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x ),
#                                                                                                           labels = trans_format("log10", math_format(10^.x)))
#     print(r)
#     
#     
#     
#     
#     
#   }
#     
# dev.off()
#      
# 
# 


# 
# numberOfTCs = data.frame(Var1 = unique(unlist(lapply(sample_TCs_injection[grep(TPs[i],names(sample_TCs_injection))],function(x) as.numeric(names(x[[2]]))))) )
# numberOfReads = lapply(sample_TCs_injection[grep(TPs[i],names(sample_TCs_injection))],function(x) data.frame(x[[2]]))
# numberOfReads$numberOfTCs = numberOfTCs
# numberOfTCs_uninjected = data.frame(Var1 = unique(unlist(lapply(sample_TCs[grep(TPs[i],names(sample_TCs))],function(x) as.numeric(names(x[[2]]))))) )
# numberOfReads_uninjected = lapply(sample_TCs[grep(TPs[i],names(sample_TCs))],function(x) data.frame(x[[2]]))
# numberOfReads$uninjected = numberOfReads_uninjected[[1]]
# 





#########################
TPs = paste0("TP",c(1:9))
library(RColorBrewer)
colsUse = c(brewer.pal(n = 8,"Dark2"),brewer.pal(n = 4,"Set1"))

pdf("//Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/errorRates_untreatedTCs/otherMutations_compareTC.pdf",height = 5,width = 5)
for(i in 1:length(TPs)){
  
  load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/FractionOfTCreads_totalNumberOfreads_InjectedSamples.Rdata")
  names(sample_TCs_injection) = substr(names(sample_TCs_injection),start = 227,stop = 240)
  
  load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/FractionOfTCreads_totalNumberOfreads_unInjectedSamples.Rdata")
  names(sample_TCs) = substr(names(sample_TCs),start = 230,stop = 243)
  
  sample_TCs_injection_total = sample_TCs_injection
  sample_TCs_injection = lapply(sample_TCs_injection,function(x) x$allOther)
  
  sample_TCs_TP = sample_TCs_injection[grep(TPs[i],names(sample_TCs_injection))] 
  sample_TCs_TP = lapply(sample_TCs_TP,function(x) x$numberOfreads)
  
 p =  lapply(sample_TCs_TP,function(x) lapply(x,function(y) data.frame(y))) %>% lapply(. %>% plyr::ldply(., rbind)) %>% 
  lapply(. %>% mutate(id = paste0(.id,"_",x))) %>% purrr::reduce(., full_join, by = 'id') %>% replace(., is.na(.), 0) %>% 
  tidyr::separate(data = .,col = id,sep="_",into=c('firstNt','lastNt','numberNt')) %>% 
      dplyr::mutate(totalReads = Freq.x + Freq.y + Freq) %>% dplyr::select(-matches('x.x|x.y|.id.x|.id.y')) %>%
          dplyr::mutate(id_final = paste(firstNt,lastNt,sep="_")) %>% dplyr::group_by(id_final) %>% dplyr::mutate(fractionReads = totalReads/sum(totalReads)) %>%
            dplyr::mutate(numberNt = as.numeric(numberNt)) 
 
 ##### 
 p_TC = p %>% filter(id_final=='T_C') %>% mutate(type="T_C")
 p_other = p %>% filter(id_final !='T_C') %>% mutate(type = 'other')
 p_total = rbind(p_TC,p_other)
 
 q = ggplot(data = p_total %>% filter(numberNt<5) ,aes(x=type,y=fractionReads,fill=type,group=type)) + geom_violin()  +   scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x ),labels = trans_format("log10", math_format(10^.x))) + theme_ameres(type = 'barplot') + ggtitle("Distribution of conversions-allGenes")
 print(q)
 p =  ggplot() + geom_point(data = p %>% filter(id_final == "T_C"),aes(x=numberNt,y=fractionReads),size=2.5)  + geom_line(data = p %>% filter(id_final == "T_C"),aes(x=numberNt,y=fractionReads),col='darkgreen',size=1.5) +  geom_boxplot(data = p %>% filter(id_final != "T_C"),aes(x=numberNt,y=fractionReads,group=numberNt),fill='red') +  
   geom_line(data = p %>% filter(id_final != "T_C"),aes(x=numberNt,y=fractionReads,group=id_final),col='grey') +   scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x ),labels = trans_format("log10", math_format(10^.x))) + theme_ameres(type = "barplot") +
   scale_color_manual(values = colsUse) + ggtitle(paste0("allGenes-",TPs[i]))  + xlab("Number of conversions") + ylab("Fraction of reads") + scale_x_continuous(breaks = c(0,1,2,3,4,5, 6,7,8,9,10)) + xlim(c(0,5))
 
 print(p)
 
 #%>% ggpubr::ggline(.,x = 'numberNt',y='fractionReads',group='id_final',col='id_final' ) +  
   
  sample_TCs_TP_mt = lapply(sample_TCs_injection_total,function(x) x$mt)
  
  
  sample_TCs_TP_mt = sample_TCs_TP_mt[grep(TPs[i],names(sample_TCs_TP_mt))] 
  
  sample_TCs_TP_mt = lapply(sample_TCs_TP_mt,function(x) x$numberOfreads)
  
  
  r = lapply(sample_TCs_TP_mt,function(x) lapply(x,function(y) data.frame(y))) %>% lapply(. %>% plyr::ldply(., rbind)) %>% 
    lapply(. %>% mutate(id = paste0(.id,"_",x))) %>% purrr::reduce(., full_join, by = 'id') %>% replace(., is.na(.), 0) %>% 
    tidyr::separate(data = .,col = id,sep="_",into=c('firstNt','lastNt','numberNt')) %>% 
    dplyr::mutate(totalReads = Freq.x + Freq.y + Freq) %>% dplyr::select(-matches('x.x|x.y|.id.x|.id.y')) %>%
    dplyr::mutate(id_final = paste(firstNt,lastNt,sep="_")) %>% dplyr::group_by(id_final) %>% dplyr::mutate(fractionReads = totalReads/sum(totalReads)) %>%
    dplyr::mutate(numberNt = as.numeric(numberNt))
  r_TC = r %>% filter(id_final=='T_C') %>% mutate(type="T_C")
  r_other = r %>% filter(id_final !='T_C') %>% mutate(type = 'other')
  r_total = rbind(r_TC,r_other)
  s = ggplot(data = r_total %>% filter(numberNt<5) ,aes(x=type,y=fractionReads,fill=type,group=type)) + geom_violin()  +   scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x ),labels = trans_format("log10", math_format(10^.x))) + ggtitle("mt-distributionFractionReads") + theme_ameres(type = 'barplot')
  print(s)
  
    r = ggplot() + geom_point(data = r %>% filter(id_final == "T_C"),aes(x=numberNt,y=fractionReads),size=2.5)  + geom_line(data = r %>% filter(id_final == "T_C"),aes(x=numberNt,y=fractionReads),col='darkgreen',size=1.5) +  geom_boxplot(data = r %>% filter(id_final != "T_C"),aes(x=numberNt,y=fractionReads,group=numberNt),fill='red') +  
    geom_line(data = r %>% filter(id_final != "T_C"),aes(x=numberNt,y=fractionReads,group=id_final),col='grey') +   scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x ),labels = trans_format("log10", math_format(10^.x))) + theme_ameres(type = "barplot") +
    scale_color_manual(values = colsUse) + ggtitle(paste0("mt-",TPs[i]))  + xlab("Number of conversions") + ylab("Fraction of reads")  + scale_x_continuous(breaks = c(0,1,2,3,4,5, 6,7,8,9,10)) + xlim(c(0,5))
  print(r)
  
  

  #### doing the same for the background data...
  
  sample_TCs_TP_bg = lapply(sample_TCs,function(x) x$allOther)
  sample_TCs_TP_bg = sample_TCs_TP_bg[grep(TPs[i],names(sample_TCs_TP_bg))] 
  sample_TCs_TP_bg = sample_TCs_TP_bg[[1]]$numberOfreads
  
  bg_all = lapply(sample_TCs_TP_bg,function(x) data.frame(x)) %>% plyr::ldply(., rbind)%>%
    dplyr::group_by(.data = .,.id) %>% dplyr::mutate(fractionReads = Freq/sum(Freq))
  
  bg_TC = bg_all %>% dplyr::filter(.id=="T_C") 
  
  bg_all = as.data.frame(bg_all)
  bg_all$x = as.numeric(bg_all$x)
  bg_all = ggplot() + geom_point(data =as.data.frame(bg_all) %>% filter(.id == "T_C"),aes(x=x,y=fractionReads),size=2.5)  + geom_line(data = as.data.frame(bg_all) %>% filter(.id == "T_C"),aes(x=x,y=fractionReads,group=.id),col='darkgreen',size=1.5) +  geom_boxplot(data = as.data.frame(bg_all) %>% filter(.id != "T_C"),aes(x=x,y=fractionReads,group=x),fill='red') +  
    geom_line(data = as.data.frame(bg_all) %>% filter(.id != "T_C"),aes(x=x,y=fractionReads,group=.id),col='grey') +   scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x ),labels = trans_format("log10", math_format(10^.x))) + theme_ameres(type = "barplot") +
    scale_color_manual(values = colsUse) + ggtitle(paste0("allGenes_bg-",TPs[i]))  + xlab("Number of conversions") + ylab("Fraction of reads")  + scale_x_continuous( limits = c(0,5),breaks = c(0,1,2,3,4,5, 6,7,8,9,10)) 
  print(bg_all)
  
  
  
  ##### plotting only the TC fractions from the background and treatment... 
  
 p_TC =  p_TC %>% dplyr::select('fractionReads','type','numberNt') %>% dplyr::mutate(category="treated") %>% dplyr::ungroup() %>% select(-id_final) %>% dplyr::ungroup() 
  bg_TC  = bg_TC %>% dplyr::ungroup() %>% plyr::rename(c('.id'='type','x'='numberNt')) %>% dplyr::mutate(category = 'untreated') %>% dplyr::select('fractionReads','type','numberNt','category')
untreated_treatedTC = rbind(as.data.frame(p_TC),as.data.frame(bg_TC))
untreated_treatedTC$numberNt = as.numeric(untreated_treatedTC$numberNt )
tcplots = ggplot(untreated_treatedTC,aes(x=numberNt,y=fractionReads,group=category,col=category)) + geom_line()+   scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x ),labels = trans_format("log10", math_format(10^.x)))+ scale_x_continuous( limits = c(0,5))
tcplots = tcplots +  ggtitle(paste0("TCfractions all genes-",TPs[i]))+ xlab("Number of conversions") + ylab("Fraction of reads") + theme_ameres(type = "barplot")  + scale_color_brewer(palette = 'Dark2')
print(tcplots)
  # ggpubr::ggline( .,x = 'x',y='fractionReads',group = '.id',col='.id') + 
  #   scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x ),labels = trans_format("log10", math_format(10^.x)))+
  #   scale_color_manual(values = colsUse)+ ggtitle(paste0("BG_allGenes-",TPs[i])) + theme_ameres(type = "barplot")  
  # print(bg_all) 
  # 
  
  sample_TCs_TP_bg_mt = lapply(sample_TCs,function(x) x$mt)
  sample_TCs_TP_bg_mt = sample_TCs_TP_bg_mt[grep(TPs[i],names(sample_TCs_TP_bg_mt))] 
  sample_TCs_TP_bg_mt = sample_TCs_TP_bg_mt[[1]]$numberOfreads
  
  bg_mt = lapply(sample_TCs_TP_bg_mt,function(x) data.frame(x)) %>% plyr::ldply(., rbind)%>%
    dplyr::group_by(.data = .,.id) %>% gplyr::mutate(fractionReads = Freq/sum(Freq)) 
  bg_mt = as.data.frame(bg_mt)
  bg_mt$x = as.numeric(bg_mt$x)
  bg_mt = ggplot() + geom_point(data =as.data.frame(bg_mt) %>% filter(.id == "T_C"),aes(x=x,y=fractionReads),size=2.5)  + geom_line(data = as.data.frame(bg_mt) %>% filter(.id == "T_C"),aes(x=x,y=fractionReads,group=.id),col='darkgreen',size=1.5) +  geom_boxplot(data = as.data.frame(bg_mt) %>% filter(.id != "T_C"),aes(x=x,y=fractionReads,group=x),fill='red') +  
    geom_line(data = as.data.frame(bg_mt) %>% filter(.id != "T_C"),aes(x=x,y=fractionReads,group=.id),col='grey') +   scale_y_continuous(trans='log10',breaks = trans_breaks("log10", function(x) 10^x ),labels = trans_format("log10", math_format(10^.x))) + theme_ameres(type = "barplot") +
    scale_color_manual(values = colsUse) + ggtitle(paste0("mt_bg-",TPs[i]))  + xlab("Number of conversions") + ylab("Fraction of reads")  + scale_x_continuous( limits = c(0,5),breaks = c(0,1,2,3,4,5, 6,7,8,9,10)) 
  
  
  print(bg_mt)
  
  
}

dev.off()


########### i would now like to calculate the over representation of TC conversions over the other conversions... 



load("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/errorRates_untreatedSamples_zebrafish/FractionOfTCreads_totalNumberOfreads_InjectedSamples.Rdata")
names(sample_TCs_injection) = substr(names(sample_TCs_injection),start = 227,stop = 250)

###################

allOtherInjection = lapply(sample_TCs_injection,function(x) lapply(x[[2]],function(y) lapply(y, function(z) z )))
allOtherInjection = lapply(allOtherInjection,function(x) x[[2]])

allOtherInjection_greaterThanEqual2  = lapply(allOtherInjection,function(x) lapply(x,function(y) sum(y[3:5],na.rm = T)))



timePoints = paste0("TP",c(1:9))

for(i in 1:length(TPs)){
 
  timepointSpecific_treated = allOtherInjection_greaterThanEqual2[grep(timePoints[i],names(allOtherInjection_greaterThanEqual2))]
  TC_timePoint = rowSums(do.call(cbind,lapply(timepointSpecific_treated,function(x) unlist(x))))["T_C"]
  #Other_timePoint = round(sum(rowSums(do.call(cbind,lapply(timepointSpecific_treated,function(x) unlist(x))))[c("A_C","A_G","A_T","C_A","C_G","C_T","G_A","G_C","T_A","T_G")]),0)
  Other_timePoint = rowSums(do.call(cbind,lapply(timepointSpecific_treated,function(x) unlist(x))))["T_G"]
  
  allOtherInjection_timeSpecific = allOtherInjection[grep(timePoints[i],names(allOtherInjection))]
  readsContainingNoConversions = lapply(allOtherInjection_timeSpecific,function(x) lapply(x,function(y) y[1]))
  readsContainingNoTCConversions = sum(unlist(lapply(readsContainingNoConversions,function(x) unlist(x)["T_C.0"])))
  #readsContainingNoOTHERConversions = round(sum(unlist(lapply(readsContainingNoConversions,function(x) unlist(x)[c("A_C.0","A_G.0","A_T.0","C_A.0","C_G.0","C_T.0","G_A.0","G_C.0","G_T.0","T_A.0","T_G.0")]))),0)
  readsContainingNoOTHERConversions = round(sum(unlist(lapply(readsContainingNoConversions,function(x) unlist(x)[c("T_G.0")]))),0)
  
  contTable = matrix(data = c(TC_timePoint,Other_timePoint,readsContainingNoTCConversions,readsContainingNoOTHERConversions),byrow = T,nrow = 2)
  fisher.test(contTable,alternative = "t") 
  
}
  