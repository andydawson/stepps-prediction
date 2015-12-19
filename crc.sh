#!/bin/sh
# options
#$ -N {short}_{run}
#$ -m be
#$ -M andria.dawson@gmail.com
#$ -pe smp {threads}
#$ -V
#$ -R y
#$ -q *@@bio

cd $HOME/Documents/projects/stepps-prediction/runs/{model}/{run}
cp $HOME/Documents/projects/stepps-prediction/cpp/{exe} .

fsync $SGE_STDOUT_PATH &
fsync output.csv &

export OMP_NUM_THREADS={threads}
./{exe} sample num_warmup={num_warmup} \
        num_samples={num_samples} \
        save_warmup={save_warmup} \
        data file=input.dump \
        output file=output.csv
