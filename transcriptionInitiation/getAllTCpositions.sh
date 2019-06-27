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

### MP TAG format : nucleotide:positionInRead:poitionInReference

ml samtools/1.9-foss-2017a 

cd /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/filter

OUTDIR="/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/filter/"


####### for the negative strand

for j in *_filtered.bam

do
	  
	cd "$OUTDIR"
	samtools view -h "$j" | grep -v TC:i:0 | sed '/^@/ d' | cut -f 2,3,4,6,24,23,22 | awk '$1==16' -  > /scratch/pooja/bam_subset.txt  ### reads on the negative strand
	awk '{print substr($7,6)}' /scratch/pooja/bam_subset.txt | paste /scratch/pooja/bam_subset.txt - | cut -f 8 - | sed 's/\,/\ \t/g' - > /scratch/pooja/bam_subset_MPseparated.txt  ### getting the MP taga, each tag separated by a '\t'
	paste /scratch/pooja/bam_subset.txt  /scratch/pooja/bam_subset_MPseparated.txt > /scratch/pooja/MPtags_separated.txt ##### just pasting to the metadata columns.

	### getting the number of columns : 
	nCols=`awk '{print NF}'  /scratch/pooja/bam_subset_MPseparated.txt  | sort -nk1,1 -  | tail -n1` 
							
	### just modify the file /scratch/pooja/MPtags_separated.txt, from the 7th column to 7+ $nCols
	mkdir -p /scratch/pooja/tmp_processing/
	touch  /scratch/pooja/tmp_processing/combined.txt

	otherCols=7
	endCols=`expr $nCols + $otherCols`
	#awk '{for(i=8;i<=NF;++i)printf $i""FS ; print "\t"}' /scratch/pooja/MPtags_separated.txt > /scratch/pooja/tmp_processing/mpTagFile.txt
	
	##### creating the reference 'Start' file.
	cut -f3 /scratch/pooja/bam_subset.txt > /scratch/pooja/tmp_processing/positions.txt ### getting the read start position
	cut -f2 /scratch/pooja/bam_subset.txt > /scratch/pooja/tmp_processing/chromosome.txt #### getting the chromosome name
	paste -d '_' /scratch/pooja/tmp_processing/chromosome.txt /scratch/pooja/tmp_processing/positions.txt > /scratch/pooja/tmp_processing/chromosome_position.txt  ### creating an identifier with chromosome name and the position of the read.
															
	for ((i = 1; i <= "$nCols"; i++)); do
																	
		cut -f "$i"  /scratch/pooja/bam_subset_MPseparated.txt | sed '/^2:/!s/.*/\t/' - | awk '{split($0,a,":"); print a[3]}' - | paste /scratch/pooja/tmp_processing/positions.txt - | awk '{for(i=1;i<=NF;i++) t+=$i;print t;t=0}' - | paste -d '_' /scratch/pooja/tmp_processing/chromosome.txt -  | paste /scratch/pooja/tmp_processing/chromosome_position.txt - | awk '$1!=$2 {print $2}' - > /scratch/pooja/tmp_processing/file"$i".txt ##### already getting the position of the SNP
		#### addinf the chromosome information.
		#### filtering out if the chromosome_position == chromomdosome_position+value
																								
	done

	cd /scratch/pooja/tmp_processing/
	### just pasting all the files together:

	STR=""
		for ((i=1; i <="$nCols"; i++)); do
			STR=$STR"file"$i."txt"" "
		done

	cat $STR > /scratch/pooja/tmp_processing/combinedFile.txt

						#### if there is only one value per column, there is a duplication of the values, i will remove these events later. 
	sort /scratch/pooja/tmp_processing/combinedFile.txt | uniq - | sed 's/_/\t/' - > TCsnpFile_minus_"$j".txt

done





######################## for the positive strand



cd "$OUTDIR"

for j in *_filtered.bam

do
	         
	cd "$OUTDIR"
	samtools view -h "$j" | grep -v TC:i:0 | sed '/^@/ d' | cut -f 2,3,4,6,24,23,22 | awk '$1!=16' -  > /scratch/pooja/bam_subset.txt 
	awk '{print substr($7,6)}' /scratch/pooja/bam_subset.txt | paste /scratch/pooja/bam_subset.txt - | cut -f 8 - | sed 's/\,/\ \t/g' - > /scratch/pooja/bam_subset_MPseparated.txt 
	paste /scratch/pooja/bam_subset.txt  /scratch/pooja/bam_subset_MPseparated.txt > /scratch/pooja/MPtags_separated.txt

	### getting the number of columns : 
	nCols=`awk '{print NF}'  /scratch/pooja/bam_subset_MPseparated.txt  | sort -nk1,1 -  | tail -n1` 
							                                                        
	### just modify the file /scratch/pooja/MPtags_separated.txt, from the 7th column to 7+ $nCols
	mkdir -p /scratch/pooja/tmp_processing/
	touch  /scratch/pooja/tmp_processing/combined.txt

	otherCols=7
	endCols=`expr $nCols + $otherCols`
	#awk '{for(i=8;i<=NF;++i)printf $i""FS ; print "\t"}' /scratch/pooja/MPtags_separated.txt > /scratch/pooja/tmp_processing/mpTagFile.txt
													        
	##### creating the reference 'Start' file.
	cut -f3 /scratch/pooja/bam_subset.txt > /scratch/pooja/tmp_processing/positions.txt
	cut -f2 /scratch/pooja/bam_subset.txt > /scratch/pooja/tmp_processing/chromosome.txt
	paste -d '_' /scratch/pooja/tmp_processing/chromosome.txt /scratch/pooja/tmp_processing/positions.txt > /scratch/pooja/tmp_processing/chromosome_position.txt
																	                                                                                                                        
	for ((i = 1; i <= "$nCols"; i++)); do
	cut -f "$i"  /scratch/pooja/bam_subset_MPseparated.txt | sed '/^16:/!s/.*/\t/' - | awk '{split($0,a,":"); print a[3]}' - | paste /scratch/pooja/tmp_processing/positions.txt - | awk '{for(i=1;i<=NF;i++) t+=$i;print t;t=0}' - | paste -d '_' /scratch/pooja/tmp_processing/chromosome.txt -  | paste /scratch/pooja/tmp_processing/chromosome_position.txt - | awk '$1!=$2 {print $2}' - > /scratch/pooja/tmp_processing/file"$i".txt ##### already getting the position of the SNP
#### addinf the chromosome information.
#### filtering out if the chromosome_position == chromomdosome_position+value
																									                                                                                                                                                                                                
done

	cd /scratch/pooja/tmp_processing/
																											        ### just pasting all the files together:

	STR=""
	for ((i=1; i <="$nCols"; i++)); do
		STR=$STR"file"$i."txt"" "
	done

	cat $STR > /scratch/pooja/tmp_processing/combinedFile.txt
																																					                                                #### if there is only one value per column, there is a duplication of the values, i will remove these events later.
															
	sort /scratch/pooja/tmp_processing/combinedFile.txt | uniq - | sed 's/_/\t/' - > TCsnpFile_plus_"$j".txt


done


######### combining the plus and minus together...

 #cd /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter

#for j in *_filtered.bam

#do
#	cat /scratch/pooja/tmp_processing/TCsnpFile_plus_"$j".txt /scratch/pooja/tmp_processing/TCsnpFile_minus_"$j".txt > /scratch/pooja/tmp_processing/TCsnpFile_total_"$j".txt
#done 





cd /scratch/pooja/tmp_processing/

mkdir /scratch/pooja/tmp_processing/finalBedFiles/


for i in TCsnpFile_plus_combinedFile**.txt
	 do

		awk -v OFS='\t' '{print $1,$2,$2+1,$1$2,0,"+"}' "$i" > /scratch/pooja/tmp_processing/finalBedFiles/"$i".bed
	done

for i in TCsnpFile_minus_combinedFile**.txt
	do
	
		awk -v OFS='\t' '{print $1,$2,$2+1,$1$2,0,"-"}' "$i" > /scratch/pooja/tmp_processing/finalBedFiles/"$i".bed
	done


#### putting the bed files together and sorting them.
	
cd "$OUTDIR"
for j in *_filtered.bam
	do
		cat /scratch/pooja/tmp_processing/finalBedFiles/TCsnpFile_plus_"$j".txt.bed	/scratch/pooja/tmp_processing/finalBedFiles/TCsnpFile_minus_"$j".txt.bed > /scratch/pooja/tmp_processing/finalBedFiles/TCsnpFile_total_"$j".txt.bed

		sort -k1,1 -k2,2n /scratch/pooja/tmp_processing/finalBedFiles/TCsnpFile_total_"$j".txt.bed > /scratch/pooja/tmp_processing/finalBedFiles/TCsnpFile_total_"$j".txt_sorted.bed
	done




