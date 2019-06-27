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

#### maing this overlap a bit faster

	#### from the bam file --> get reads that overlp with interesting counting findows
	#### from the SNP file, get the SNPs that overlap with the interesting couting windows
	#### get the reads that from the subsetted bam file that overlap the interesting SNPs.



ml bedtools/2.27.1-foss-2017a  
ml samtools/1.9-foss-2017a
######### overlap with counting windows first... 

cd /scratch/pooja/tmp_processing_old/tmp_processing/finalBedFiles/

mkdir -p /scratch/pooja/bedFiles_split/

### splitting the counting window file


awk '$5 == "+"' /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/countinWindows_expressedAccordingToTime.bed > /scratch/pooja/Cws_plus.bed
awk '$5 == "-"' /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/countinWindows_expressedAccordingToTime.bed > /scratch/pooja/Cws_minus.bed

##$###### subset bam file to get reads that overlap the counting windows of interest






#### also subset the SNP file to get SNPs that overlap the counting windows

#for i in *lamdunk_mapped_filtered.bam.txt_sorted.bed
#	do
			
#			awk '$6 == "-"' "$i"  > /scratch/pooja/bedFiles_split/"$i"_minus.bed
#			awk '$6 == "+"' "$i" > /scratch/pooja/bedFiles_split/"$i"_plus.bed
			
#			bedtools intersect -b /scratch/pooja/Cws_plus.bed  -a /scratch/pooja/bedFiles_split/"$i"_plus.bed > /scratch/pooja/bedFiles_split/"$i"_plus_overlapping.bed
#			bedtools intersect -b /scratch/pooja/Cws_minus.bed  -a /scratch/pooja/bedFiles_split/"$i"_minus.bed > /scratch/pooja/bedFiles_split/"$i"_minus_overlapping.bed	
			
#			awk -F' ' '{print $1":" $2"-" $3}' /scratch/pooja/bedFiles_split/"$i"_plus_overlapping.bed > /scratch/pooja/bedFiles_split/"$i"_plus_overlapping_positions.bed
#			awk -F' ' '{print $1":" $2"-" $3}' /scratch/pooja/bedFiles_split/"$i"_minus_overlapping.bed > /scratch/pooja/bedFiles_split/"$i"_minus_overlapping_positions.bed
#	done









##### step 3: overlap the identified positions subsets of bam files with the positions recorded



#### first for the minus strand 
cd  /scratch/pooja/bedFiles_split/
mkdir -p  /scratch/pooja/bedFiles_split/coverage/
OUTDIR=/scratch/pooja/bedFiles_split/coverage/
        
OUTDIR_FINAL="$OUTDIR"/finalBed/
mkdir -p "$OUTDIR_FINAL"
#for i in *_minus_overlapping_positions.bed
#	        do
#			while read p
						                                                        
#				do
				
			
#					samtools view  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/finalBamFile_R3.bam  "$p" | awk '$2==16' - | wc -l >> "$OUTDIR"/"$i"_R3.txt
#					samtools view  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/finalBamFile_R2.bam  "$p" | awk '$2==16' - | wc -l >> "$OUTDIR"/"$i"_R2.txt
#					samtools view  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/finalBamFile.bam  "$p" | awk '$2==16' - | wc -l >> "$OUTDIR"/"$i"_R1.txt

#			done <  /scratch/pooja/bedFiles_split/"$i"
														        
	

																                                                                                 
#																done

###### plus strand


for i in *_plus_overlapping_positions.bed
	     do
		 while read p
								                                                                                                        
			do
													                                
													                        
			samtools view  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/finalBamFile_R3.bam  "$p" | awk '$2!=16' - | wc -l >> "$OUTDIR"/"$i"_R3.txt
																		                                        samtools view  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/finalBamFile_R2.bam  "$p" | awk '$2!=16' - | wc -l >> "$OUTDIR"/"$i"_R2.txt
																							                   samtools view  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/finalBamFile.bam  "$p" | awk '$2!=16' - | wc -l >> "$OUTDIR"/"$i"_R1.txt
																											                                        
																												                        done <  /scratch/pooja/bedFiles_split/"$i"


																									done
