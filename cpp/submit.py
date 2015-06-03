import subprocess

dry_run = True

system = 'crc'

runs = [ { 'name':       'pred_kw_kgamma_PL',
           'exe':        'pred_kw_kgamma_262.exe',
           'data':       '12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_PL_umw_3by_v0.3.dump',
           'output':     '12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_PL_umw_3by_v3.dump',
           'num_warmup':  75,
           'num_samples': 1000,
           'save_warmup': 1,
           'threads':     12,
           },
         ]

qsub = {

    'crc': """\
#!/bin/sh
# options
#$ -N {name}
#$ -m be
#$ -M andria.dawson@gmail.com
#$ -pe smp {threads}
#$ -V
#$ -R y
#$ -q *@@bio

cd $HOME/Documents/projects/stepps-prediction/cpp

fsync $SGE_STDOUT_PATH &
fsync ../output/{output} &

export OMP_NUM_THREADS={threads}
./{exe} sample num_warmup={num_warmup} \
        num_samples={num_samples} \
        save_warmup={save_warmup} \
        data file=../r/dump/{data} \
        output file=../output/{output} \
        random seed=42 
""",

    'scf-slurm': """\
#!/bin/sh
#SBATCH --job-name={name}
#SBATCH -c {threads}
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andria.dawson@gmail.com

cd $HOME/Documents/projects/stepps-prediction/cpp
export OMP_NUM_THREADS={threads}
srun ./{exe} sample num_warmup={num_warmup} \
        num_samples={num_samples} \
        save_warmup={save_warmup} \
        data file=../r/dump/{data} \
        output file=../output/{output} \
        random seed=42 
""",

    'scf-sge': """\
""",

}

submit = {
    'crc': 'qsub',
    'scf-slurm': 'srun',
    'scf-sge': 'qsub',
    }

for run in runs:
    sub = qsub[system].format(**run)
    with open(run['name'] + ".sh", 'w') as f:
        f.write(sub)
    print "submitting:", run['name']
    if not dry_run:
        subprocess.check_call([submit[system], run['name'] + '.sh'])
    else:
        print sub
