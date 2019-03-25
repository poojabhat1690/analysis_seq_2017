#!/usr/bin/env bash                                                                                                                                                                
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=100G
#SBATCH --time=0-24:00:00     # 2 minutes
#SBATCH --job-name=VarScan_mergedSamples
####SBATCH --array=1-4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pooja.bhat@imba.oeaw.ac.at
#SBATCH  --qos=long

### identifying SNPs using bcf tools
ml samtools/1.9-foss-2017a 
module load bcftools/1.6-foss-2017a  
ml bedtools/2.27.1-foss-2017a
#### i want to know all the conversions and the coverage over them

OUTDIR="/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/snpFiltering"

mkdir -p "$OUTDIR"

#bcftools mpileup -O v -I -R /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/allAnnotations_withPriMirna.bed -f /groups/ameres/bioinformatics/references/danio_rerio/dr11/danRer11.fa  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/filter/finalBamFile_R1.bam  > "$OUTDIR"/finalBamFile_R1.vcf

# bcftools mpileup -O v -I -R /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/allAnnotations_withPriMirna.bed  -f /groups/ameres/bioinformatics/references/danio_rerio/dr11/danRer11.fa  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/filter/finalBamFile_R2.bam  > "$OUTDIR"/finalBamFile_R2.vcf

 #bcftools mpileup -O v -I -R /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/allAnnotations_withPriMirna.bed -f /groups/ameres/bioinformatics/references/danio_rerio/dr11/danRer11.fa  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/filter/finalBamFile_R3.bam  > "$OUTDIR"/finalBamFile_R3.vcf



######## now i can get the all Ts and all As  that are covered by <10 reads in the total combined files... 

awk '$6 == "+"' /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/allAnnotations_withPriMirna.bed > /scratch/pooja/Cws_plus.bed
awk '$6 == "-"' /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/allAnnotations_withPriMirna.bed > /scratch/pooja/Cws_minus.bed



replicates=(R1 R2 R3)
 

#for i in "${replicates[@]}"
#do
#	cat  "$OUTDIR"/finalBamFile_"$i".vcf  | awk '$4 == "T" || $4=="A"' - | awk  '{split($8,a,";"); print $0,a[1]}'  | awk '{split($NF,a,"="); print $0,a[2]}' | awk  -v OFS="\t" '{print $1,$2,$4,$5,$12}' - | awk '$5<10' - | awk -v OFS="\t" '{split($4,a,","); print $0,a[1]}' | awk -v OFS="\t" '{print $0,$3 $6}'  |  awk -v s=1 '{print $1, $2,$2+s,$7,$5, $3}' -  | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' | awk '$6 == "A"' - | bedtools intersect -a "stdin" -b /scratch/pooja/Cws_minus.bed | awk -v OFS="\t" '{print $1, $2, $3, $4,$5, "-"}' - > /scratch/pooja/"$i"_minusOverlapping.bed


#	  cat  "$OUTDIR"/finalBamFile_"$i".vcf  | awk '$4 == "T" || $4=="A"' - | awk  '{split($8,a,";"); print $0,a[1]}'  | awk '{split($NF,a,"="); print $0,a[2]}' | awk  -v OFS="\t" '{print $1,$2,$4,$5,$12}' - | awk '$5<10' - | awk -v OFS="\t" '{split($4,a,","); print $0,a[1]}' | awk -v OFS="\t" '{print $0,$3 $6}'  |  awk -v s=1 '{print $1, $2,$2+s,$7,$5, $3}' -  | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' | awk '$6 == "T"' - | bedtools intersect -a "stdin" -b /scratch/pooja/Cws_plus.bed | awk -v OFS="\t" '{print $1, $2, $3, $4,$5, "+"}' - > /scratch/pooja/"$i"_plusOverlapping.bed


#	  cat /scratch/pooja/"$i"_minusOverlapping.bed /scratch/pooja/"$i"_plusOverlapping.bed > "$OUTDIR"/finalBamFile_"$i"_position_lessthan10.bed
#done





 grep R1 /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt | cut -f2 > /scratch/pooja/R1samples.txt
 grep R2 /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt | cut -f2 > /scratch/pooja/R2samples.txt
 grep R3 /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt | cut -f2 > /scratch/pooja/R3samples.txt
 grep Unt /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt | cut -f2 > /scratch/pooja/Untreatedsamples.txt

DIR="/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/filter/"
OP="/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/filter_revised/"
mkdir -p "$OP" 

for i in "${replicates[@]}"
do

	while read p; do

		bedtools intersect -a "$DIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/snpFiltering//finalBamFile_"$i"_position_lessthan10.bed  -v >  "$OP"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam
	       samtools index  "$OP"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OP"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam.bai
	 
       done <  /scratch/pooja/"$i"samples.txt	
done


##### the untreated sample comes from R1, so just removing R1 from the untreated samples. 


  while read p; do

	                  bedtools intersect -a "$DIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/snpFiltering//finalBamFile_R1_position_lessthan10.bed  -v >  "$OP"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam
			                 samtools index  "$OP"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OP"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam.bai
					          
done <  /scratch/pooja/Untreatedsamples.txt 
