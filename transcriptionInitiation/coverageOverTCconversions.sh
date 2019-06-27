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

ml bedtools/2.27.1-foss-2017a  

cd /scratch/pooja/tmp_processing/finalBedFiles/

mkdir /scratch/pooja/tmp_processing/finalBedFiles/coverage/
COVERAGE=/scratch/pooja/tmp_processing/finalBedFiles/coverage/

############### coverage around all TC positions..  from combined bam files (all sample together per genotype).


#for j in *adapterTrimmed_slamdunk_mapped_filtered.bam.txt_sorted.bed
	#do
	#bedtools coverage -d -sorted -s -a /scratch/pooja/tmp_processing/finalBedFiles/"$j" -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/finalBamFile.bam > "$COVERAGE"/"$j"_R1.bed
	#bedtools coverage -d -sorted -s -a /scratch/pooja/tmp_processing/finalBedFiles/"$j" -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/finalBamFile_R2.bam > "$COVERAGE"/"$j"_R2.bed
	#bedtools coverage -d -sorted -s -a /scratch/pooja/tmp_processing/finalBedFiles/"$j" -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/finalBamFile_R3.bam > "$COVERAGE"/"$j"_R3.bed

	#done


cd "$COVERAGE"
mkdir -p "$COVERAGE"/filtered/
mkdir -p "$COVERAGE"/filtered/overlap_countingWindows/

for i in *.bed

do
	#awk '$8 < 10' "$i" > "$COVERAGE"/filtered/"$i"

	### intersectiong the low coverage TC conversions with counting windows
	bedtools intersect -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/allAnnotations_withPriMirna.bed -a "$COVERAGE"/filtered/"$i" -wa -wb -s > "$COVERAGE"/filtered/overlap_countingWindows/"$i"

done

