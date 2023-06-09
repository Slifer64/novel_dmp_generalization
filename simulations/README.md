# Simulations

Containts scripts that run the simulations and plots the corresponding figures. In each case, the figures may not be concatenated as in the paper, or some additional plots may also appear.

The code was tested in Ubuntu 20.04 using MatlabR18b.

## Simulations and results/plots
Run the following matlab scripts for each Figure:

- **Figure 1a**: 
```
sim_same_demo_start_goal
```
- **Figure 1b**: 
```
sim_same_exec_start_end.m
```
- **Figure 1c**: 
```
sim_mirroring.m
```
- **Figure 2** : For each of the 4 cases run: 
```
compare_scalings(target_id)
```
where `target_id` is 1, 2, 3 or 4.
- **Figures 3-4** :
```
viapoint_sim('DMP_pp') % for the proposed
viapoint_sim('VMP') % for VMP
```
- **Figure 6** : 
```
sim_ineq_constr.m
```
Requires that QSQP (matlab interface) is installed (see https://osqp.org/docs/get_started/matlab.html for installation)
