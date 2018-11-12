#!/usr/bin/env bash                                                                                                                                                                
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-24:00:00     # 2 minutes
#SBATCH --job-name=runSLAMdunk
####SBATCH --array=1-4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pooja.bhat@imba.oeaw.ac.at
#SBATCH  --qos=long

        
while getopts 'i: o: d: s:' OPTION; do
          case "$OPTION" in   
                     i)        
                            ivalue="$OPTARG"
                            echo "the input file  $OPTARG"
                            ;;
                    o) 
                    		ovalue_adapter="$OPTARG"
                    		echo "the output directory for adapter trimmed files is $OPTARG"
                    		;;       
                    d)
                    		indir="$OPTARG"
                    		echo "the input directory is $OPTARG"
                    		;;                
                    
                    s) 
                    		sdir="$OPTARG"
                    		echo "the SLAMdunk optput is in $OPTARG"
				;;
                                                                         
                                                           ?)
                                                                	echo "script usage: $(basename $0) [-a adapter]" >&2
                                                                	exit 1
                                                                	;;
                    esac
                                        done
                shift "$(($OPTIND -1))"


mkdir -p "$sdir"

module load cutadapt/1.16-foss-2017a-python-2.7.13
module load  python/2.7.13-foss-2017a
module load python/2.7.13-foss-2017a
module load cutadapt/1.9.1-foss-2017a-python-2.7.13
module load fastx-toolkit/0.0.14-foss-2017a
module load fastqc/0.11.5-java-1.8.0_121
module load anaconda2/5.1.0
module load  bedtools/2.25.0-foss-2017a
module load r/3.4.1-foss-2017a-x11-20170314




cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -o "$ovalue_adapter"/"$ivalue"_adapterTrimmed.fastq -m 18 --trim-n "$indir"/"$ivalue"



               export PYTHONNOUSERSITE=1
                 source activate slamdunk0.3.3

        
       

			slamdunk all -r /groups/ameres/bioinformatics/references/danio_rerio/dr11/danRer11.fa -o "$sdir" -b /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/allAnnotations_withPriMirna.bed  -fb /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/allStagesCombined_new/countingWindows_transcriptionalOutput_priMirna.bed -a 5 -5 12 -n 100 -mv 0.2 -t 15 -mq 0 -mi 0.95 -m -rl 88 "$ovalue_adapter"/"$ivalue"_adapterTrimmed.fastq

       
       
                 source deactivate
        
