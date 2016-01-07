#!/bin/sh
#SBATCH --job-name={short}
#SBATCH -c {threads}
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andria.dawson@gmail.com
#SBATCH --exclude=scf-sm20
#SBATCH --workdir=/accounts/gen/vis/adawson/Documents/projects/stepps-prediction/runs/{model}/{run}
#SBATCH --error="slurm-%j.err"
#SBATCH --output="slurm-%j.out"

cp $HOME/Documents/projects/stepps-prediction/cpp/{cpp} .
cp $HOME/Documents/projects/stepps-prediction/cpp/{exe} .

export OMP_NUM_THREADS={threads}
srun ./{exe} sample num_warmup={num_warmup} \
        num_samples={num_samples} \
        save_warmup={save_warmup} \
        data file=input.dump \
        output file=output.csv