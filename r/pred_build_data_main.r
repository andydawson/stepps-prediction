grids = c(1, 3, 5)
sides = c('', 'E', 'W')

run_g_all = list(suff_fit  = 'cal_g_ALL_v0.3', 
                 suff_dat = '12taxa_mid_comp_ALL_v0.2',
                 kernel    = 'gaussian', 
                 one_psi   = TRUE, 
                 one_gamma = TRUE, 
                 EPs       = FALSE)
run_pl_all = list(suff_fit  = 'cal_pl_ALL_v0.3', 
                  suff_dat = '12taxa_mid_comp_ALL_v0.2',
                  kernel    = 'pl', 
                  one_a     = TRUE,
                  one_b     = TRUE,
                  one_gamma = TRUE, 
                  EPs       = FALSE)
run_g_Kpsi_Kgamma = list(suff_fit  = 'cal_g_Kpsi_Kgamma_EPs_v0.3', 
            suff_dat = '12taxa_mid_comp_v0.1',
            kernel    = 'gaussian', 
            one_psi   = FALSE, 
            one_gamma = FALSE, 
            EPs       = TRUE)
run = list(suff_fit  = 'cal_g_Kpsi_EPs_v0.3', 
            suff_dat = '12taxa_mid_comp_v0.1',
            kernel    = 'gaussian', 
            one_psi   = FALSE, 
            one_gamma = TRUE, 
            EPs       = TRUE)

runs = list(run_g_all, run_pl_all)
grids = c(3)


for (run in runs){
  for (res in grids){
    for (side in sides){
      source('r/pred_build_data.r')
    }
  }
}