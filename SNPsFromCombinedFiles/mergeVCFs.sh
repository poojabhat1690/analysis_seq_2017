#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=100G
#SBATCH --time=0-5:00:00     # 2 minutes
#SBATCH --job-name=mergeVCFFiles
####SBATCH --array=1-4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pooja.bhat@imba.oeaw.ac.at

### this script merges VCF files from diff

INDIR=/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/filter/

### adding back the genomtype specific SNPs to individual files... 

 grep R1 /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt | cut -f2 > /scratch/pooja/R1samples.txt
  grep R2 /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt | cut -f2 > /scratch/pooja/R2samples.txt
   grep R3 /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt | cut -f2 > /scratch/pooja/R3samples.txt
   grep Unt /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/sampleInfo/barcodes_description.txt | cut -f2 > /scratch/pooja/Untreatedsamples.txt

SNPDIR=/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/snp

cd "$SNPDIR"

OUTDIR=/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/snp_combined/

mkdir -p "$OUTDIR"


while read p; do

	 cat combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf "$INDIR"/varSample_R1_minCov.vcf > "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf
	echo "$p" "finished"  
	                

	awk -F"\t" '!seen[$1, $2, $3, $4]++' "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf > "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf_tmp


	mv "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf_tmp "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf





done <  /scratch/pooja/R1samples.txt



### R2



while read p; do

         cat combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf "$INDIR"/varSample_R2_minCov.vcf > "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf
	         echo "$p" "finished"  
	awk -F"\t" '!seen[$1, $2, $3, $4]++' "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf > "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf_tmp
	
	mv "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf_tmp "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf

done <  /scratch/pooja/R2samples.txt



### R3




while read p; do

	         cat combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf "$INDIR"/varSample_R3_minCov.vcf > "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf
		                  echo "$p" "finished"         
				 
				  awk -F"\t" '!seen[$1, $2, $3, $4]++' "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf > "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf_tmp
				          
				          mv "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf_tmp "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf

				  done <  /scratch/pooja/R3samples.txt





#### untreated 



while read p; do

	         cat combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf "$INDIR"/varSample_R1_minCov.vcf > "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf
		  
		 echo "$p" "finished"  
			                         

		awk -F"\t" '!seen[$1, $2, $3, $4]++' "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf > "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf_tmp


		mv "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf_tmp "$OUTDIR"/combinedFile_"$p".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_snp.vcf


done <  /scratch/pooja/Untreatedsamples.txt





	




