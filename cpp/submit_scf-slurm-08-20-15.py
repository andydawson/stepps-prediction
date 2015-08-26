import subprocess

runs = [ ('pred_PL',
          './pred_base_262.exe \
          sample num_warmup=100 num_samples=1000 save_warmup=1\
          data file=../r/dump/12taxa_699cells_120knots_50to2050ypb_PL_umw_3by_v0.3_bacon.dump \
          output file=../output/12taxa_699cells_120knots_50to2050ypb_PL_umw_3by_v3_bacon.dump \
          random seed=42'),
         ('pred_G',
          './pred_base_262.exe \
          sample num_warmup=100 num_samples=1000 save_warmup=1\
          data file=../r/dump/12taxa_699cells_120knots_50to2050ypb_G_umw_3by_v0.3_bacon.dump \
          output file=../output/12taxa_699cells_120knots_50to2050ypb_G_umw_3by_v2_bacon.dump \
          random seed=42'),
         ('pred_kw_kgamma_PL',
          './pred_kw_kgamma_262.exe \
          sample num_warmup=100 num_samples=1000 save_warmup=1\
          data file=../r/dump/12taxa_699cells_120knots_50to2050ypb_KW_KGAMMA_PL_umw_3by_v0.3_bacon.dump \
          output file=../output/12taxa_699cells_120knots_50to2050ypb_KW_KGAMMA_PL_umw_3by_v3_bacon.dump \
          random seed=42'),
         ('pred_kw_kgamma_G',
          './pred_kw_kgamma_262.exe \
          sample num_warmup=100 num_samples=1000 save_warmup=1\
          data file=../r/dump/12taxa_699cells_120knots_50to2050ypb_KW_KGAMMA_G_umw_3by_v0.3_bacon.dump \
          output file=../output/12taxa_699cells_120knots_50to2050ypb_KW_KGAMMA_G_umw_3by_v2_bacon.dump \
          random seed=42')
]

 # ('pred_G_3by_ALL',
 #          './pred_od_mpp_full_nug_mu0.exe \
 #          sample num_warmup=75 num_samples=1000 save_warmup=1\
 #          data file=../r/dump/12taxa_699cells_120knots_0to2000ypb_G_umw_3by_v0.3.dump \
 #          output file=../output/12taxa_699cells_120knots_0to2000ypb_G_umw_3by.csv\
 #          random seed=42'),

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