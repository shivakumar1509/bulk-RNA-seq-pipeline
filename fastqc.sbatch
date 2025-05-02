#!/bin/bash
#SBATCH -J shivrv  # job name
#SBATCH -o log_shivrv.o%j    # output and error file name (%j expands to jobID)
#SBATCH -e log_shivrv.e%j    # output and error file name (%j expands to jobID)
#SBATCH -n 1               # Run one process
#SBATCH --cpus-per-task=28 # with a full 28 cores available
#SBATCH -p defq            # queue (partition) -- defq, ipowerq, eduq, gpuq.
#SBATCH -t 0-8:00:00      # run time (d-hh:mm:ss)
#SBATCH --account=mjurynec
#SBATCH --partition=lonepeak
ulimit -v unlimited
ulimit -s unlimited
ulimit -u 10000

date
date
date

# Execute the program:
fastqc -o ~/experiments/Mick_mice_RNA_seq/fastqc_out/ -f fastq ~/experiments/Mick_mice_RNA_seq/rawdata/*.gz

date
date
date

echo "DONE---------------------"