#!/usr/bin/env bash                                                                                                                                 #SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=7      
#SBATCH --mem-per-cpu=20G   
#SBATCH --time=0-24:00:00     # 2 minutes
#SBATCH --job-name=runSLAMdunk  
####SBATCH --array=1-4          
#SBATCH --mail-type=ALL         
#SBATCH --mail-user=pooja.bhat@imba.oeaw.ac.at
#SBATCH  --qos=long   


module load cutadapt/1.16-foss-2017a-python-2.7.13
module load  python/2.7.13-foss-2017a
module load python/2.7.13-foss-2017a
module load cutadapt/1.9.1-foss-2017a-python-2.7.13
module load fastx-toolkit/0.0.14-foss-2017a
module load fastqc/0.11.5-java-1.8.0_121
module load anaconda2/5.1.0
module load  bedtools/2.25.0-foss-2017a
module load r/3.4.1-foss-2017a-x11-20170314
ml singularity/2.5.2 


cd /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/filter_revised/
 
                  
		
 		   for i in combined**.bam
				  do
		                    
	                  
			 cd          
			   
sbatch	singularity exec ./slamdunk-v0.3.4.simg  slamdunk count -o /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/count_revised/ -s /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/snp_combined/ -r /groups/ameres/bioinformatics/references/danio_rerio/dr11/danRer11.fa -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/allAnnotations_withPriMirna.bed -l 88 -q 27 /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/filter_revised/"$i"

###  also simultaneously i want the count of reads with 2TC   

sbatch	singularity exec ./slamdunk-v0.3.4.simg slamdunk count -o /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/count_2tc_revised/ -s /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/snp_combined/ -r /groups/ameres/bioinformatics/references/danio_rerio/dr11/danRer11.fa -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/allAnnotations_withPriMirna.bed -l 88 -q 27 -c 2 /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/filter_revised/"$i"

		

#sbatch  singularity exec ./slamdunk-v0.3.4.simg  slamdunk count -o /scratch/pooja/count_withAllSNPs_noGenotypePreference/ -s /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/snp_combined_withoutGenotype/ -r /groups/ameres/bioinformatics/references/danio_rerio/dr11/danRer11.fa -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/allAnnotations_withPriMirna.bed -l 88 -q 27 /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_0.3.4/filter_revised/"$i"

done   
	
						                                              
