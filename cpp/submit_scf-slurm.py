import subprocess

runs = [ ('pred_G_3by_ALL',
          './pred_od_mpp_full_nug_mu0.exe \
          sample num_warmup=75 num_samples=1000 save_warmup=1\
          data file=../r/dump/12taxa_699cells_120knots_0to2000ypb_G_umw_3by_v0.3.dump \
          output file=../output/12taxa_699cells_120knots_0to2000ypb_G_umw_3by.csv\
          random seed=42'),
         ('pred_PL_1by_ALL',
          './pred_od_mpp_full_nug_mu0_262_test.exe \
          sample num_warmup=75 num_samples=1000 save_warmup=1\
          data file= ../r/dump/12taxa_6341cells_120knots_0to2000ypb_PL_umw_1by_v0.3.dump\           
          output file=../output/12taxa_6341cells_120knots_0to2000ypb_PL_umw_3by.csv\
          random seed=42')
]

qsub = """\
#!/bin/sh
#SBATCH --job-name={name}
#SBATCH -c {threads}
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andria.dawson@gmail.com

cd $HOME/Documents/projects/stepps-prediction/cpp
export OMP_NUM_THREADS={threads}
srun {command}
"""

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
