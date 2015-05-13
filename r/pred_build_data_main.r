grids = c(1, 3, 5)
sides = c('', 'E', 'W')

run_g = list(suff_fit  = 'cal_g_ALL_v0.3', 
             suff_dat = '12taxa_mid_comp_ALL_v0.2',
             kernel    = 'gaussian', 
             one_psi   = TRUE, 
             one_gamma = TRUE, 
             EPs       = FALSE)
run_pl = list(suff_fit  = 'cal_pl_ALL_v0.3', 
              suff_dat = '12taxa_mid_comp_ALL_v0.2',
              kernel    = 'pl', 
              one_a     = TRUE,
              one_b     = TRUE,
              one_gamma = TRUE, 
              EPs       = FALSE)
run_g_Kpsi_Kgamma = list(suff_fit  = 'cal_g_Kpsi_Kgamma_EPs_ALL_v0.3', 
                         suff_dat = '12taxa_mid_comp_v0.1',
                         kernel    = 'gaussian', 
                         one_psi   = FALSE, 
                         one_gamma = FALSE, 
                         EPs       = TRUE)
run_pl_Ka_Kgamma = list(suff_fit  = 'cal_pl_Ka_Kgamma_EPs_ALL_v0.3', 
                        suff_dat  = '12taxa_mid_comp_ALL_v0.2',
                        kernel    = 'pl', 
                        one_a     = FALSE, 
                        one_b     = TRUE,
                        one_gamma = FALSE, 
                        EPs       = TRUE)
# 
# run = list(suff_fit  = 'cal_g_Kpsi_EPs_v0.3', 
#             suff_dat = '12taxa_mid_comp_v0.1',
#             kernel    = 'gaussian', 
#             one_psi   = FALSE, 
#             one_gamma = TRUE, 
#             EPs       = TRUE)

runs = list(run_g, run_pl, run_g_Kpsi_Kgamma, run_pl_Ka_Kgamma)
runs = list(run_pl, run_pl_Ka_Kgamma)
grids = c(3)

for (run in runs){
  for (res in grids){
    for (side in sides){
      source('r/pred_build_data.r')
    }
  }
}
