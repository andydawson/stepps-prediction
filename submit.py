
import itertools
from os.path import join as pjoin

#
# config
#

#system  = 'scf-slurm' 
system  = 'crc' #'scf
models  = [
   { 'model':      '120knots_150to20150ybp_PL_umw_3by_v2.1_mean_ar',
     'short':      'pred_kw_kgamma_arfv_nb',
     'exe':        'pred_kw_kgamma_arfv_nb_262.exe',
     'cpp':        'pred_kw_kgamma_arfv_nb_262.cpp',   
     'num_warmup':  200,
     'num_samples': 2000,
     'save_warmup': 1,
     'threads':     12,
   }# ,
   # { 'model':      '120knots_50to2050ybp_KW_KGAMMA_PL_umw_3by_v0.5_bacon_mean_lambda',
   #   'short':      'pred_kw_kgamma_lambda_pl',
   #   'exe':        'pred_kw_kgamma_lambda_262.exe',
   #   'cpp':        'pred_kw_kgamma_lambda_262.cpp',   
   #   'num_warmup':  250,
   #   'num_samples': 2000,
   #   'save_warmup': 1,
   #   'threads':     12,
   # }#,
  #   { 'model':      '120knots_50to2050ybp_KW_KGAMMA_PL_umw_3by_v0.4_bacon_mean',
  #     'short':      'pred_kw_kgamma_pl',
  #     'exe':        'pred_kw_kgamma_262.exe',
  #     'cpp':        'pred_kw_kgamma_262.cpp',
  #     'num_warmup':  200,
  #     'num_samples': 1000,
  #     'save_warmup': 1,
  #     'threads':     12,
  # }#,
]

submit = {
    'crc': 'qsub',
    'scf-slurm': 'sbatch',
    }

#
#
#

with open(system + '.sh', 'r') as f:
    sub = f.read()

with open('submit-all.sh', 'w') as all:
    for model, run in itertools.product(models, range(1)):
        fname = pjoin('runs', model['model'], 'run' + str(run+1), 'submit.sh')
        with open(fname, 'w') as f:
            f.write(sub.format(run='run'+str(run+1), **model))
        print "wrote:", fname
        all.write(submit[system] + " " + fname + '\n')

print 'wrote: submit-all.sh'
print ''
print 'run "submit-all.sh" on {system} to queue all runs'.format(system=system)
