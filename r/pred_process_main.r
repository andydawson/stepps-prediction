# PL = list(suff_dat = '12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_PL_umw_3by_v0.3',
#           suff_fit = '12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_PL_umw_3by_v0.3')
# 
# G  = list(suff_dat = '12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_G_umw_3by_v0.3', 
#           suff_fit = '12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_G_umw_3by_v0.3')

PL = list(suff_dat  = '12taxa_699cells_120knots_0to2000ypb_PL_umw_3by_v0.3',
          suff_fit  = '12taxa_699cells_120knots_0to2000ypb_PL_umw_3by',
          suff_figs = 'PL')

G  = list(suff_dat  = '12taxa_699cells_120knots_0to2000ypb_G_umw_3by_v0.3', 
          suff_fit  = '12taxa_699cells_120knots_0to2000ypb_G_umw_3by', 
          suff_figs = 'G')

# PL = list(suff_dat = '12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_PL_umw_3by_v0.3',
#           suff_fit = '12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_PL_umw_3by_v2_tmp',
#           suff_figs = 'KW_KGAMMA_PL')
# 
# G  = list(suff_dat = '12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_G_umw_3by_v0.3', 
#           suff_fit = '12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_G_umw_3by_v2_tmp',
#           suff_figs = 'KW_KGAMMA_G')

runs = list(PL, G)

for (run in runs){
  suff_dat  = run$suff_dat
  suff_fit  = run$suff_fit
  suff_figs = run$suff_figs
  source('r/pred_process_full_test.r')
  source('r/pred_plot.r')
  
#   post = load_stan_output(suff_fit)
  
}