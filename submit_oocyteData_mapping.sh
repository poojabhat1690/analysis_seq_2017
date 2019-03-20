#!/bin/bash

cd /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/rawData/oocyte_quantSeq_fromQuio/

for files in *.gz
		do
				
mkdir -p  /scratch/pooja/mapping_oocyteData/
mkdir -p /scratch/pooja/mapping_dr11_oocyteData/

sbatch  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/scripts/map_ii2.sh -i "$files" -o /scratch/pooja/mapping_oocyteData/ -d /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/rawData/oocyte_quantSeq_fromQuio/ -c /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/allAnnotations_withPriMirna.bed -p /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/countingWindows_transcriptionalOutput_priMirna.bed -s /scratch/pooja/mapping_dr11_oocyteData/

done




