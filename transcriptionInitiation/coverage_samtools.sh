#!/usr/bin/env bash                                                                                                                                                                                                            
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-24:00:00     # 2 minutes
#SBATCH --job-name=allTCpositions
####SBATCH --array=1-4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pooja.bhat@imba.oeaw.ac.at
#SBATCH  --qos=long

ml samtools/1.9-foss-2017a 
mkdir -p /scratch/pooja/coverage_allTC/

### take the negative SNPs. 
cd /scratch/pooja/tmp_processing_old/tmp_processing/finalBedFiles/

### just separating th bed file into positive and negative


for i in  *nk_mapped_filtered.bam.txt_sorted.bed
do
		cat "$i" |  awk '$6 == "-"' - | awk -F" " '{print $1":" $2"-" $3}' - > /scratch/pooja/coverage_allTC/"$i"_negative.bed
		cat "$i" | awk '$6 == "+"' - | awk -F" " '{print $1":" $2"-" $3}'  - > /scratch/pooja/coverage_allTC/"$i"_positive.bed
done

##### ]



mkdir -p /scratch/pooja/coverage_allTC/coverage_negative/

cd /scratch/pooja/coverage_allTC/ 

for i in *_negative.bed
	do
		while read p
							
			do
	samtools view  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/finalBamFile_R3.bam "$p" | awk '$2==16' - | wc -l >> /scratch/pooja/coverage_allTC/coverage_negative/"$i"_R2.txt
	samtools view  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/finalBamFile_R2.bam "$p" | awk '$2==16' - | wc -l >>/scratch/pooja/coverage_allTC/coverage_negative/"$i"_R3.txt
	samtools view  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/finalBamFile.bam "$p" | awk '$2==16' - | wc -l >>/scratch/pooja/coverage_allTC/coverage_negative/"$i"_R1.txt
																				
	done <  /scratch/pooja/coverage_allTC/"$i"

	paste  /scratch/pooja/coverage_allTC/"$i"_negative.bed /scratch/pooja/coverage_allTC/coverage_negative/"$i"_R2.txt > /scratch/pooja/coverage_allTC/coverage_negative/"$i"_R2_withCounts.txt
	paste  /scratch/pooja/coverage_allTC/"$i"_negative.bed /scratch/pooja/coverage_allTC/coverage_negative/"$i"_R2.txt > /scratch/pooja/coverage_allTC/coverage_negative/"$i"_R3_withCounts.txt

	paste  /scratch/pooja/coverage_allTC/"$i"_negative.bed /scratch/pooja/coverage_allTC/coverage_negative/"$i"_R1.txt > /scratch/pooja/coverage_allTC/coverage_negative/"$i"_R1_withCounts.txt
																								
awk -F'\t' 'print $1 $2 "-"' /scratch/pooja/coverage_allTC/coverage_negative/"$i"_R2_withCounts.txt > /scratch/pooja/coverage_allTC/coverage_negative/"$i"_R2_withCounts_final.txt

awk -F'\t' 'print $1 $2 "-"' /scratch/pooja/coverage_allTC/coverage_negative/"$i"_R3_withCounts.txt > /scratch/pooja/coverage_allTC/coverage_negative/"$i"_R3_withCounts_final.txt

awk -F'\t' 'print $1 $2 "-"' /scratch/pooja/coverage_allTC/coverage_negative/"$i"_R1_withCounts.txt > /scratch/pooja/coverage_allTC/coverage_negative/"$i"_R1_withCounts_final.txt

done



##### plus 



mkdir -p /scratch/pooja/coverage_allTC/coverage_positive/

cd /scratch/pooja/coverage_allTC/ 

for i in *_positive.bed
																								do
																									while read p
																										do
																									samtools view  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/finalBamFile_R3.bam "$p" | awk '$2!=16' - | wc -l >> /scratch/pooja/coverage_allTC/coverage_positive/"$i"_R2.txt
																									samtools view  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/finalBamFile_R2.bam "$p" | awk '$2!=16' - | wc -l >>/scratch/pooja/coverage_allTC/coverage_positive/"$i"_R3.txt
																									samtools view  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/finalBamFile.bam "$p" | awk '$2!=16' - | wc -l >>/scratch/pooja/coverage_allTC/coverage_positive/"$i"_R1.txt
																																
			done <  /scratch/pooja/coverage_allTC/"$i"

																																														paste  /scratch/pooja/coverage_allTC/"$i"_positive.bed /scratch/pooja/coverage_allTC/coverage_positive/"$i"_R2.txt > /scratch/pooja/coverage_allTC/coverage_positive/"$i"_R2_withCounts.txt
																						
																																														paste  /scratch/pooja/coverage_allTC/"$i"_positive.bed /scratch/pooja/coverage_allTC/coverage_positive/"$i"_R2.txt > /scratch/pooja/coverage_allTC/coverage_positive/"$i"_R3_withCounts.txt
																								paste  /scratch/pooja/coverage_allTC/"$i"_positive.bed /scratch/pooja/coverage_allTC/coverage_positive/"$i"_R1.txt > /scratch/pooja/coverage_allTC/coverage_positive/"$i"_R1_withCounts.txt
																																																	
																																																	awk -F'\t' 'print $1 $2 "+"' /scratch/pooja/coverage_allTC/coverage_positive/"$i"_R2_withCounts.txt > /scratch/pooja/coverage_allTC/coverage_positive/"$i"_R2_withCounts_final.txt
																																																		awk -F'\t' 'print $1 $2 "+"' /scratch/pooja/coverage_allTC/coverage_positive/"$i"_R3_withCounts.txt > /scratch/pooja/coverage_allTC/coverage_positive/"$i"_R3_withCounts_final.txt
																																																			awk -F'\t' 'print $1 $2 "+"' /scratch/pooja/coverage_allTC/coverage_positive/"$i"_R1_withCounts.txt > /scratch/pooja/coverage_allTC/coverage_positive/"$i"_R1_withCounts_final.txt

		done

