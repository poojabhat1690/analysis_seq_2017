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



while getopts 'i: o: d: c: p: s:' OPTION; do
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

		    c)
			    CWs="$OPTARG"
			    echo "the counting windows used are in file $OPTARG"
			    ;;

		    p)
			   MMAP="$OPTARG"
			  echo "the windows used for multimapping are $OPTARG"
			    ;;
                    
                    s) 
                   		sdir="$OPTARG"
                    		echo "the SLAMdunk optput is in $OPTARG"
				;;
                                                                         
                     ?)
                               echo "script usage: $(basename $0) " >&2
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
module load python/2.7.13-foss-2017a
module load r/3.5.1-foss-2017a-x11-20170314-bare
module load java/1.8.0_121
module load singularity/2.5.2 


if [ -e "$ovalue_adapter"/"$ivalue"_adapterTrimmed.fastq ]; then
	echo "adapter trimmed file exists."
else
	echo "performing adapter clipping"

	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -o "$ovalue_adapter"/"$ivalue"_adapterTrimmed.fastq -m 18 --trim-n "$indir"/"$ivalue"
fi




               #export PYTHONNOUSERSITE=1
               
	       #source activate slamdunk0.3.4

        	#slamdunk --version 
	       	#slamdunk all -r /groups/ameres/bioinformatics/references/danio_rerio/dr11/danRer11.fa -o "$sdir" -b "$CWs"  -fb "$MMAP" -a 5 -5 12 -n 100 -mv 0.2 -t 1 -mq 0 -mi 0.95 -m -rl 88 "$ovalue_adapter"/"$ivalue"_adapterTrimmed.fastq

  
         	#source deactivate
        



		#### using the singularity module:

	cd
		singularity exec ./slamdunk-v0.3.4.simg slamdunk all -r /groups/ameres/bioinformatics/references/danio_rerio/dr11/danRer11.fa -o "$sdir" -b "$CWs"  -fb "$MMAP" -a 5 -5 12 -n 100 -mv 0.2 -t 1 -mq 0 -mi 0.95 -m -rl 88 "$ovalue_adapter"/"$ivalue"_adapterTrimmed.fastq

