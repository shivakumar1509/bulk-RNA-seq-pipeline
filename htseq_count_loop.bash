#!/bin/bash
#SBATCH -J SRaya  # job name
#SBATCH -o log_srayavara.o%j    # output and error file name (%j expands to jobID)
#SBATCH -n 1               # Run one process
#SBATCH --cpus-per-task=28 # with a full 28 cores available
#SBATCH -p defq            # queue (partition) -- defq, ipowerq, eduq, gpuq.
#SBATCH -t 0-72:00:00      # run time (d-hh:mm:ss)
ulimit -v unlimited
ulimit -s unlimited
ulimit -u 10000


# Execute the program:
conda activate py2
/data/shiv/mick_mice/samfiles/
for file in ./*.sam; do echo ${file}; htseq-count ${file} --type=gene --idattr=gene_id --additional-attr=gene_name --stranded=no //data/shiv/mick_mice/mice_ref_genome/Mus_musculus.GRCm39.105.gtf > ${file}.HTSEQ.txt; done

cd /data/shiv/samfiles
/data/shiv/samfiles/*.sam; do echo ${file}; htseq-count ${file} --type=gene --idattr=gene_id --additional-attr=gene_name --stranded=no /data/shiv/human_ref_genome/Homo_sapiens.GRCh38.104.gtf > ${file}.HTSEQ.txt; done


conda activate py2
/data/shiv/mick_mice/samfiles/
for file in ./*.sam; do echo ${file}; htseq-count ${file} --type=gene --idattr=gene_id --additional-attr=gene_name --stranded=no //data/shiv/mick_mice/mice_ref_genome/Mus_musculus.GRCm39.105.gtf > ${file}.HTSEQ.txt; done