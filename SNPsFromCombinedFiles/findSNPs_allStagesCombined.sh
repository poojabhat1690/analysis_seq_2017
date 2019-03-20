#!/usr/bin/env bash                                                                                                                                                                
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=100G
#SBATCH --time=0-5:00:00     # 2 minutes
#SBATCH --job-name=VarScan_mergedSamples
####SBATCH --array=1-4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pooja.bhat@imba.oeaw.ac.at




cd /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/
OUTDIR=/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11/filter/

module load samtools/1.4-foss-2017a


samtools merge "$OUTDIR"/finalBamFile.bam  "$OUTDIR"/combinedFile_CGGTTA.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_TTAACT.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_ATGAAC.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_CCTAAG.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_AATCCG.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_GGCTGC.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_TACCTT.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_TCTTAA.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_GTCAGG.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_ATACTG.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_TATGTC.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_GAGTCC.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_GGAGGT.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_CACACT.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_CCGCAA.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_TTTATG.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_AACGCC.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_CAAGCA.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam

samtools index "$OUTDIR"/finalBamFile.bam  

samtools mpileup -B -A -l /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/countingWindows_transcriptionalOutput_priMirna.bed -f /groups/ameres/bioinformatics/references/danio_rerio/dr11/danRer11.fa  "$OUTDIR"/finalBamFile.bam >  "$OUTDIR"/finalBam_R1_pileup.mpileup
cd /groups/ameres/bioinformatics/tools/slamdunk_new-2/slamdunk/slamdunk/contrib/


java -jar VarScan.v2.4.1.jar  mpileup2snp "$OUTDIR"/finalBam_R1_pileup.mpileup  --min-coverage 10 --min-var-freq 0.20 --strand-filter 0 --variants 1 --output-vcf 1 > "$OUTDIR"/varSample_R1_minCov.vcf


### replicate 2

samtools merge  "$OUTDIR"/finalBamFile_R2.bam  "$OUTDIR"/combinedFile_GCTCGA.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_GCGAAT.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_TGGATT.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_ACCTAC.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_CGAAGG.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_AGATAG.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_TTGGTA.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_GTTACC.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_CGCAAC.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_TGGCGA.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_ACCGTG.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_CAACAG.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_GATTGT.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_CTCTCG.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_TGACAC.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_AAGACA.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_ACAGAT.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam  "$OUTDIR"/combinedFile_TAGGCT.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam

#samtools index "$OUTDIR"/finalBamFile_R2.bam

samtools mpileup -B -A -l /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/countingWindows_transcriptionalOutput_priMirna.bed -f /groups/ameres/bioinformatics/references/danio_rerio/dr11/danRer11.fa  "$OUTDIR"/finalBamFile_R2.bam >  "$OUTDIR"/finalBam_R2_pileup.mpileup
cd /groups/ameres/bioinformatics/tools/slamdunk_new-2/slamdunk/slamdunk/contrib/

java -jar VarScan.v2.4.1.jar  mpileup2snp "$OUTDIR"/finalBam_R2_pileup.mpileup  --min-coverage 10 --min-var-freq 0.20 --strand-filter 0 --variants 1 --output-vcf 1 > "$OUTDIR"/varSample_R2_minCov.vcf

##replica 3
samtools merge "$OUTDIR"/finalBamFile_R3.bam "$OUTDIR"/combinedFile_CTCCAT.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam "$OUTDIR"/combinedFile_GCATGG.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam "$OUTDIR"/combinedFile_AATAGC.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam "$OUTDIR"/combinedFile_GTGCCA.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam "$OUTDIR"/combinedFile_TCGAGG.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam "$OUTDIR"/combinedFile_CACTAA.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam "$OUTDIR"/combinedFile_GGTATA.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam "$OUTDIR"/combinedFile_CGCCTG.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam "$OUTDIR"/combinedFile_AATGAA.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam "$OUTDIR"/combinedFile_ACAACG.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam "$OUTDIR"/combinedFile_ATATCC.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam "$OUTDIR"/combinedFile_AGTACT.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam "$OUTDIR"/combinedFile_ATAAGA.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam "$OUTDIR"/combinedFile_GGTGAG.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam "$OUTDIR"/combinedFile_TTCCGC.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam "$OUTDIR"/combinedFile_GAAGTG.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam "$OUTDIR"/combinedFile_CAATGC.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam "$OUTDIR"/combinedFile_ACGTCT.fastq.gz_adapterTrimmed_slamdunk_mapped_filtered.bam


samtools index "$OUTDIR"/finalBamFile_R3.bam

samtools mpileup -B -A -l /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/countingWindows_transcriptionalOutput_priMirna.bed -f /groups/ameres/bioinformatics/references/danio_rerio/dr11/danRer11.fa  "$OUTDIR"/finalBamFile_R3.bam >  "$OUTDIR"/finalBam_R3_pileup.mpileup
cd /groups/ameres/bioinformatics/tools/slamdunk_new-2/slamdunk/slamdunk/contrib/


java -jar VarScan.v2.4.1.jar  mpileup2snp "$OUTDIR"/finalBam_R3_pileup.mpileup  --min-coverage 10 --min-var-freq 0.20 --strand-filter 0 --variants 1 --output-vcf 1 > "$OUTDIR"/varSample_R3_minCov.vcf



####### i want to now convert the VCF files to bed files... 

ml  bedops/2.4.30

vcf2bed < "$OUTDIR"/varSample_R1_minCov.vcf > "$OUTDIR"/varSample_R1_minCov.bed
vcf2bed < "$OUTDIR"/varSample_R2_minCov.vcf > "$OUTDIR"/varSample_R2_minCov.bed
vcf2bed < "$OUTDIR"/varSample_R3_minCov.vcf > "$OUTDIR"/varSample_R3_minCov.bed


cat "$OUTDIR"/varSample_R2_minCov.bed "$OUTDIR"/varSample_R1_minCov.bed "$OUTDIR"/varSample_R3_minCov.bed > "$OUTDIR"/varSample_all_minCov.bed
cut -f1-7 "$OUTDIR"/varSample_all_minCov10.bed  | sort -k1,1 -k2,2n - | uniq - > "$OUTDIR"/allSNPS_combinedReplicated.bed


