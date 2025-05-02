Cd ..

#!/bin/bash
#SBATCH -J shivrv  # job name
#SBATCH -o log_shivrv.o%j    # output and error file name (%j expands to jobID)
#SBATCH -n 1               # Run one process
#SBATCH --cpus-per-task=28 # with a full 28 cores available
#SBATCH -p defq            # queue (partition) -- defq, ipowerq, eduq, gpuq.
#SBATCH -t 0-8:00:00      # run time (d-hh:mm:ss)
ulimit -v unlimited
ulimit -s unlimited
ulimit -u 10000


# Execute the program:
hisat2-build /data/shiv/mick_mice/mice_ref_genome/GRCm39_r105.all.fa /data/shiv/mick_mice/mice_index_output/GRCm39_r105.all
