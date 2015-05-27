
import subprocess

runs = [ ('pred_bench',
          './pred_base_262.exe \
          sample num_warmup=5 num_samples=5 save_warmup=1 \
          data file=../r/dump/12taxa_699cells_120knots_0to2000ypb_G_umw_3by_v0.3.dump \
          output file=12taxa_699cells_120knots_0to2000ypb_G_umw_3by.csv \
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
