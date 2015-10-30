
import itertools
from os.path import join as pjoin

#
# config
#

system  = 'crc'#'scf-slurm'
models  = [
   # { 'model':      '120knots_50to2050ybp_G_umw_3by_v0.4_bacon',
   #   'short':      'pred_g',
   #   'exe':        'pred_base_262.exe',
   #   'num_warmup':  75,
   #   'num_samples': 1000,
   #   'save_warmup': 1,
   #   'threads':     12,
   #},
    { 'model':      '120knots_50to2050ybp_KW_KGAMMA_PL_umw_3by_v0.4_bacon',
      'short':      'pred_kw_kgamma',
      'exe':        'pred_kw_kgamma_262.exe',
      'num_warmup':  150,
      'num_samples': 1000,
      'save_warmup': 1,
      'threads':     12,
  }#,
]

submit = {
    'crc': 'qsub',
    'scf-slurm': 'srun',
    }

#
#
#

with open(system + '.sh', 'r') as f:
    sub = f.read()

with open('submit-all.sh', 'w') as all:
    for model, run in itertools.product(models, range(10)):
        fname = pjoin(model['model'], 'run' + str(run+1), 'submit.sh')
        with open(fname, 'w') as f:
            f.write(sub.format(run='run'+str(run+1), **model))
        print "wrote:", fname
        all.write(submit[system] + " " + fname + '\n')
        # if not dry_run:
        #     subprocess.check_call([submit[system], run['name'] + '.sh'])
        # else:
        #     print sub

print 'wrote: submit-all.sh'
print ''
print 'run "submit-all.sh" on {system} to queue all runs'.format(system=system)
