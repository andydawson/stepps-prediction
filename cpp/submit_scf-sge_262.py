
import subprocess

runs = [ ('pred_1by_ALL',
          './pred_od_mpp_full_nug_mu0_262.exe \
          sample num_warmup=175 num_samples=1000 \
          data file=../r/dump/12taxa_4199cells_79knots_0to2000ypb_PL_umwW_1by_v0.3.dump \
          output file=../output/12taxa_4199cells_79knots_0to2000ypb_PL_umwW_1by_v0.3.csv \
          random seed=42')
]

         # data file=../r/dump/12taxa_6341cells_120knots_0to2000ypb_PL_umw_1by_v0.3.dump \
         #  output file=../output/12taxa_6341cells_120knots_0to2000ypb_PL_umw_1by_v0.3.csv \

qsub = """
#!/bin/sh
# options
#$ -N {name}
#$ -m be
#$ -M andria.dawson@gmail.com
#$ -pe smp {threads}
#$ -q low.q
#$ -l h_rt=672:00:00
#$ -V
#$ -l hostname=scf-sm12
#$ -R y

cd $HOME/Documents/projects/stepps-prediction/cpp
export OMP_NUM_THREADS={threads}
{command}
"""

dry_run = False

for name, command in runs:
#    sub = qsub.format(queue="low.q", walltime="672:00:00", command=run, threads=1)
    sub = qsub.format(command=command, threads=12, name=name)
    with open(name + ".qsub", 'w') as f:
        f.write(sub)
    print "submitting:", name
    if not dry_run:
        subprocess.check_call(['qsub', name + '.qsub'])
    else:
        print sub
