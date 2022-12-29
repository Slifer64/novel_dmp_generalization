# Simulations

Containts scripts that run the simulations and plot the corresponding figures. In each case, the figures may not be concatenated as in the paper, or some additional plots may also appear.

The code was tested in Ubuntu 20.04 using MatlabR18b.

## Simulations
Run the following matlab scripts for each Figure:

- **Figure 1a**: `fig1a_sim_same_demo_start_goal.m`
- **Figure 1b**: `fig1b_sim_same_exec_start_end.m`
- **Figure 1c**: `fig1c_sim_mirroring.m`
- **Figure 2** : `fig2_sim_compare_scalings.m`
- **Figure 3** : `fig3_sim_variable_target.m`  
- **Figure 4** : `fig4_sim_couplings.m` \
  (set line 31 `adapt_to_ext_signals=false` to simulate the simple DMP)
- **Figure 6** : `fig6_sim_ineq_constr.m` \
  Requires that QSQP (matlab interface) is installed (see https://osqp.org/docs/get_started/matlab.html for installation)
