#!/bin/bash

cd /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/rawData/combinedFiles

for files in *.gz
		do
				

			  sbatch /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/scripts/map_ii2.sh -i "$files" -o /scratch/pooja/mapping_dr10_06022018/ -d /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/rawData/combinedFiles/ -c /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/allAnnotations_withPriMirna.bed -p /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/countingWindows_transcriptionalOutput_priMirna.bed -s /scratch/pooja/mapping_dr11/

done




