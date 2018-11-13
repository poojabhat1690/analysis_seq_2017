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
		

		

#		 awk  '$3!=""'   "$OUTDIR"/"$i"_relavant_readsWithMPtags.txt   | awk '{split($3,a,"MP:Z:"); print a[3] a[2] a[1]}' -  |  awk -F"," '$1=$1' OFS="\t" -  > "$OUTDIR"/"$i"_relavant_readsMPtagsSeparated.txt

		#samtools view "$OUTDIR"/"$i"_relavant.bam  | grep MP: - | awk -v OFS='\t' '{print $3,$4, $24, $2}' - > "$OUTDIR"/"$i"_relavantCols.txt 

	#awk '$4==16' "$OUTDIR"/"$i"_relavantCols.txt > "$OUTDIR"/"$i"_relavantCols_minus.txt ### splitting into plus and minus strand reads and splitting the MP tags separately
	#awk '$4!=16' "$OUTDIR"/"$i"_relavantCols.txt > "$OUTDIR"/"$i"_relavantCols_plus.txt
	 
	#awk '{split($3,a,"MP:Z:"); print a[3] a[2] a[1] }' "$OUTDIR"/"$i"_relavantCols_minus.txt  |  awk -F"," '$1=$1' OFS="\t" - > "$OUTDIR"/"$i"_MPtags_split_minus.txt 
	#awk '{split($3,a,"MP:Z:"); print a[3] a[2] a[1]}' "$OUTDIR"/"$i"_relavantCols_plus.txt  |  awk -F"," '$1=$1' OFS="\t" - > "$OUTDIR"/"$i"_MPtags_split_plus.txt             

	#awk '{print $1, $2}' "$OUTDIR"/"$i"_relavantCols_minus.txt  > "$OUTDIR"/"$i"_metadata_minus.txt ######keep metadata
	#awk '{print $1, $2}' "$OUTDIR"/"$i"_relavantCols_plus.txt  > "$OUTDIR"/"$i"_metadata_plus.txt




		 MAXCOLMINUS=`awk '{if (NF > max) {max = NF; line=$0}} END{print line}' "$OUTDIR"/"$i"_MPtags_split_minus.txt  | awk '{print NF}'`
		 MAXCOLPLUS=`awk '{if (NF > max) {max = NF; line=$0}} END{print line}' "$OUTDIR"/"$i"_MPtags_split_plus.txt  | awk '{print NF}'` 


		for ((j=1;j<="$MAXCOLMINUS";j+=1))
		do
			awk -v f="$j" -F '\t' '{print $f}' "$OUTDIR"/"$i"_MPtags_split_minus.txt | awk '{split($1,a,":"); print a[3]}' -    > "$OUTDIR"/"$i"_MPtags_split_minus"$j".txt #### get the reference position from MP tag
			awk -v f="$j" -F '\t' '{print $f}' "$OUTDIR"/"$i"_MPtags_split_minus.txt | awk '{split($1,a,":"); print a[1]}' -     > "$OUTDIR"/"$i"_MPtags_identity_minus"$j".txt #### get SNP identity from MP tag
			paste  "$OUTDIR"/"$i"_metadata_minus.txt "$OUTDIR"/"$i"_MPtags_split_minus"$j".txt "$OUTDIR"/"$i"_MPtags_identity_minus"$j".txt | awk '$3 >=0' -  > "$OUTDIR"/"$i"_MPtag_SNPidentity_minus"$j".txt #### adding the tag and identity together.
			  ####also adding tge metadata
		done


		#for ((j=1;j<="$MAXCOLPLUS";j+=1))
		#do
		 # 	awk -v f="$j" -F '\t' '{print $f}' "$OUTDIR"/"$i"_MPtags_split_plus.txt | awk '{split($1,a,":"); print a[3]}'  -   > "$OUTDIR"/"$i"_MPtags_split_plus"$j".txt
		#	awk -v f="$j" -F '\t' '{print $f}' "$OUTDIR"/"$i"_MPtags_split_plus.txt | awk '{split($1,a,":"); print a[1]}'  -   > "$OUTDIR"/"$i"_MPtags_identity_plus"$j".txt
		#	paste  "$OUTDIR"/"$i"_metadata_plus.txt "$OUTDIR"/"$i"_MPtags_split_plus"$j".txt "$OUTDIR"/"$i"_MPtags_identity_plus"$j".txt | awk '$3 >=0' -  > "$OUTDIR"/"$i"_MPtag_SNPidentity_plus"$j".txt #### adding the tag and identity together.
		#done


		### combine all minus and all plus together. 
		### create bed format file, which also contains the information about strand and the identity of the SNP
		### remove duplicated entries from the bed file. 
		### combine plus and minus files.

		minus=-
		plus=+
		cat "$OUTDIR"/"$i"_MPtag_SNPidentity_minus* | awk -v OFS='\t' {'print $1, $2 + $3, $2 + $3 + 1 , $3, $4'} - | awk -F"\t" '!seen[$1, $2, $3]++' -   > "$OUTDIR"/"$i"_combined_minus.txt 
		cat "$OUTDIR"/"$i"_MPtag_SNPidentity_plus* | awk -v OFS='\t'  {'print $1, $2 + $3, $2 + $3 + 1 , $3, $4'} - |  awk -F"\t" '!seen[$1, $2, $3]++' -  > "$OUTDIR"/"$i"_combined_plus.txt
		sed -i "s/$/\t$minus/" "$OUTDIR"/"$i"_combined_minus.txt  
		sed -i "s/$/\t$plus/" "$OUTDIR"/"$i"_combined_plus.txt  
		
		cat "$OUTDIR"/"$i"_combined_minus.txt  "$OUTDIR"/"$i"_combined_plus.txt |  awk -F"\t" '!seen[$1, $2, $3]++' - > "$OUTDIR"/"$i"_combined_total.txt

	 	sort -k1,1 -k2,2n "$OUTDIR"/"$i"_combined_total.txt  > "$OUTDIR"/"$i"_combined_total_sorted.bed
		
	

done




##### now i have all the SNP positions and i need to know what os the coverage across these in the merged bam files .... 



		#### so first converting the merged bam files to bed files... 

ml bedtools/2.27.1-foss-2017a  

cd /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/


######################### converting merged bam files to bed files and sorting by chromosome and position. 


#for i in finalBam*.bam

#do

#	bedtools bamtobed -i "$i" > "$OUTDIR"/finalMP/"$i".bed

#	sort -k1,1 -k2,2n "$OUTDIR"/finalMP/"$i".bed > "$OUTDIR"/finalMP/"$i"_sorted.bed
	






#done


































		
















								                













									                                                                                 



