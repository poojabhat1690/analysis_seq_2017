#!/usr/bin/env bash                                                                                                                                                                  #SBATCH --nodes=1                  #SBATCH --ntasks=1
##SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=100G
#SBATCH --time=0-8:00:00     # 2 minutes
#SBATCH --job-name=bamReadsOverlappingSNPs
####SBATCH --array=1-4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pooja.bhat@imba.oeaw.ac.at

set -x



### i want to find the coverage across all mis match positions identified.. 

	#### find reads that overlap with counting windows.
	#### find all positions with mismatches in each sample. 
	#### pool all the mismatches and takr the unique ones. 
	#### find the coverage across the positions in each sample. 
	#### check if the sum is > 100 or not. 


ml   samtools/1.4-foss-2017a

######### find reads that overlap with counting windows using samtools. ##############


OUTDIR=/scratch/pooja/coverageOverMismatches/
mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR"/finalMP/

cd /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/

for i in combined*.bam
	        
do     

	        ### getting a bam file that overlaps with the SNPs... 
		
#		samtools view -b  -L /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/allAnnotations_withPriMirna.bed  /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/"$i" > "$OUTDIR"/"$i"_relavant.bam

#		samtools view	"$OUTDIR"/"$i"_relavant.bam  | awk -v OFS='\t' '{print $3,$4, $24}' - > "$OUTDIR"/"$i"_relavant_readsWithMPtags.txt

#		 awk  '$3!=""'   "$OUTDIR"/"$i"_relavant_readsWithMPtags.txt   | awk '{split($3,a,"MP:Z:"); print a[3] a[2] a[1]}' -  |  awk -F"," '$1=$1' OFS="\t" -  > "$OUTDIR"/"$i"_relavant_readsMPtagsSeparated.txt

#		samtools view "$OUTDIR"/"$i"_relavant.bam  | grep MP: - | awk -v OFS='\t' '{print $3,$4, $24}' - > "$OUTDIR"/"$i"_relavantCols.txt 

#		 awk '{split($3,a,"MP:Z:"); print a[3] a[2] a[1]}' "$OUTDIR"/"$i"_relavantCols.txt  |  awk -F"," '$1=$1' OFS="\t" - > "$OUTDIR"/"$i"_MPtags_split.txt

		 MAXCOL=`awk '{if (NF > max) {max = NF; line=$0}} END{print line}' "$OUTDIR"/"$i"_MPtags_split.txt  | awk '{print NF}'` 


		for ((j=1;j<="$MAXCOL";j+=1))

		 do

#		 	echo "$j"

#			awk -v f="$j" -F '\t' '{print $f}' "$OUTDIR"/"$i"_MPtags_split.txt | awk '{split($1,a,":"); print a[3]}'   >"$OUTDIR"/"$i"_MPtags_split_"$j".txt


			###### i also want to get the identity of these SNPs
#			awk -v f="$j" -F '\t' '{print $f}' "$OUTDIR"/"$i"_MPtags_split.txt | awk '{split($1,a,":"); print a[1]}'   >"$OUTDIR"/"$i"_identityOfSNP_split_"$j".txt

#			paste "$OUTDIR"/"$i"_relavantCols.txt "$OUTDIR"/"$i"_MPtags_split_"$j".txt | awk '{ print $2+$4; }' - | paste "$OUTDIR"/"$i"_relavantCols.txt - | awk '$2!=$4' - | awk '{print $1, $4}' > "$OUTDIR"/finalMP/SNP_"$i"_"$j".txt

	paste "$OUTDIR"/"$i"_relavantCols.txt "$OUTDIR"/"$i"_MPtags_split_"$j".txt | awk '{ print $2+$4; }' - | paste "$OUTDIR"/"$i"_relavantCols.txt - | awk '$2!=$4' - | awk '{print $1, $4}' - | paste  - "$OUTDIR"/"$i"_identityOfSNP_split_"$j".txt | awk -v OFS='\t' '{print $1, $2, $3}' - | awk '$3!=""' -  > "$OUTDIR"/finalMP/SNPidentity_"$i"_"$j".txt

		done

	


#	cat "$OUTDIR"/finalMP/SNP_"$i"_* | awk -F"\t" '!seen[$1, $2]++' - > "$OUTDIR"/finalMP/"$i"_allSNPs_unique.txt

#	awk '{$3=$2+1}1' "$OUTDIR"/finalMP/"$i"_allSNPs_unique.txt | awk  'OFS="\t" {print $1,$2,$3}' - > "$OUTDIR"/finalMP/"$i"_allSNPs_unique.bed

#	 sort -k1,1 -k2,2n "$OUTDIR"/finalMP/"$i"_allSNPs_unique.bed  >"$OUTDIR"/finalMP/"$i"_allSNPs_unique_sorted.bed

	###I also want to add the information of identity of the SNP..  

	cat "$OUTDIR"/finalMP/SNPidentity_"$i"_*.txt | awk -F"\t" '!seen[$1, $2]++' - > "$OUTDIR"/finalMP/SNPidentity_combined_"$i".txt    ### getting the unique lines... 
	awk '{$4=$2+1}1' "$OUTDIR"/finalMP/SNPidentity_combined_"$i".txt | awk  'OFS="\t" {print $1,$2,$4, $3}' - > "$OUTDIR"/finalMP/"$i"_allSNPs_unique.bed


done




##### now i have all the SNP positions and i need to know what os the coverage across these in the merged bam files .... 



		#### so first converting the merged bam files to bed files... 

ml bedtools/2.27.1-foss-2017a  

cd /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/


#for i in finalBam*.bam

#do

#	bedtools bamtobed -i "$i" > "$OUTDIR"/finalMP/"$i".bed

#	sort -k1,1 -k2,2n "$OUTDIR"/finalMP/"$i".bed > "$OUTDIR"/finalMP/"$i"_sorted.bed
	






#done


#also convertin indicidual bam files to bed files 

#mkdir -p "$OUTDIR"/finalMP_bamFiles/

#cd /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/

#for i in  combinedFile_*.bam

#	do
#		bedtools bamtobed -i "$i" > "$OUTDIR"/finalMP_bamFiles/"$i".bed
#	done

####### now used bed tools coverage to calculate reads spanning the identified SNPS.. 


mkdir "$OUTDIR"/finalMP/coverage/
mkdir /scratch/pooja/R1_coverage/



ml bedtools/2.27.1-foss-2017a  



	cd "$OUTDIR"/finalMP/

	#for i in  combin*allSNPs_unique_sorted.bed
	#do 
	#	name=a
	#	sed -i "s/$/\t$name/"  "$OUTDIR"/finalMP/"$i"
	#	score=0
	#	sed -i "s/$/\t$score/"  "$OUTDIR"/finalMP/"$i"
	#	strand=+
	#	sed -i "s/$/\t$strand/"  "$OUTDIR"/finalMP/"$i"
	#	echo  "$i"
		
	#bedtools coverage  -a "$i" -b "$OUTDIR"/finalMP/finalBamFile.bam_sorted.bed  > "$OUTDIR"/finalMP/coverage/R1/"$i"_coverage.bed

	#done




#mkdir "$OUTDIR"/finalMP/coverage/R2/


 #       for i in  combin*_allSNPs_unique_sorted.bed     


#		 do

#				bedtools coverage -d -a "$i" -b finalBamFile_R2.bam.bed  > "$OUTDIR"/finalMP/coverage/R2/"$i"_coverage.bed
								                
#		done



#mkdir "$OUTDIR"/finalMP/coverage/R3/


 #       for i in  combin*_allSNPs_unique_sorted.bed     


#		                 do

#					bedtools coverage -d -a "$i" -b finalBamFile_R3.bam.bed  > "$OUTDIR"/finalMP/coverage/R3/"$i"_coverage.bed
									                                                                                 
#				done

