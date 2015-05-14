PL = list(suff_dat = '12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_PL_umw_3by_v0.3',
          suff_fit = '12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_PL_umw_3by_v0.3')

G  = list(suff_dat = '12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_G_umw_3by_v0.3', 
          suff_fit = '12taxa_699cells_120knots_0to2000ypb_KW_KGAMMA_G_umw_3by_v0.3')

runs = list(PL, G)

for (run in runs){
  suff_dat = run$suff_dat
  suff_fit = run$suff_fit
  source('r/pred_process_full_test.r')
}