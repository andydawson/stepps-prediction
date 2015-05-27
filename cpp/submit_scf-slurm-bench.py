import subprocess

runs = [ ('pred_bench',
          './pred_base_262.exe \
          sample num_warmup=5 num_samples=5 save_warmup=1\
          data file=../r/dump/12taxa_699cells_120knots_0to2000ypb_G_umw_3by_v0.3.dump \
          output file=12taxa_699cells_120knots_0to2000ypb_G_umw_3by.dump \
          random seed=42')
]

qsub = """\
#!/bin/sh
#SBATCH --job-name={name}
#SBATCH -c {threads}
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andria.dawson@gmail.com
#SBATCH --exclude=scf-sm20

cd $HOME/Documents/projects/stepps-prediction/cpp
export OMP_NUM_THREADS={threads}
srun {command}
"""
#SBATCH --nodelist=scf-sm21,scf-sm22,scf-sm23

dry_run = False

for name, command in runs:
#    sub = qsub.format(queue="low.q", walltime="672:00:00", command=run, threads=1)
    sub = qsub.format(command=command, threads=12, name=name)
    with open(name + ".sh", 'w') as f:
        f.write(sub)
    print "submitting:", name
    if not dry_run:
        subprocess.check_call(['sbatch', name + '.sh'])
    else:
        print sub
