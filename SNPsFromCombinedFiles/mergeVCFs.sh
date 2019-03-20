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

INDIR=/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/

############ just adding these together 

cat "$INDIR"/varSample_R1_minCov.vcf "$INDIR"/varSample_R2_minCov.vcf "$INDIR"/varSample_R3_minCov.vcf > "$INDIR"/allConversions.vcf

sed '/#/d' "$INDIR"/allConversions.vcf  > "$INDIR"/allConversions_cleaned.vcf 


##### removing duplicates... 

awk -F"\t" '!seen[$1, $2, $3, $4]++' "$INDIR"/allConversions_cleaned.vcf > "$INDIR"/allConversions_uniqueSNPs.vcf

#### now i want to add the 'total' vcf files to the file specific vcf files...just loop over the file... 

SNPDIR=/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/snp/

cd "$SNPDIR"

OUTDIR=/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/snp_combined/

mkdir -p "$OUTDIR"

for i in combinedFile**.vcf

	do

		cat "$i" "$INDIR"/allConversions_uniqueSNPs.vcf > "$OUTDIR"/"$i"	
		echo "$i" done	
		
		###### now I will again remove duplicates.. 


		awk -F"\t" '!seen[$1, $2, $3, $4]++' "$OUTDIR"/"$i" > "$OUTDIR"/"$i"_tmp


		mv "$OUTDIR"/"$i"_tmp "$OUTDIR"/"$i"
	
	done



