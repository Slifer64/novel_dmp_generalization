clc;
% close all;
clear;

%% =============  includes...  =============
addpath('utils/')
import_gmp_lib();

%% =============  Load demo data  =============                           
load('data/demo_data1.mat', 'Timed', 'Pd_data');

n_dof = size(Pd_data, 1);
% calc vel/accel numerically for plotting them later
dPd_data = [diff(Pd_data, 1, 2)./diff(Timed) zeros(n_dof,1)];
ddPd_data = [diff(dPd_data, 1, 2)./diff(Timed) zeros(n_dof,1)];

%% ============= Train DMP  =============
gmp = trainDMP(Timed, Pd_data);


%% =============  DMP simulation  =============
disp('Simulation...');
t_start = tic;

%% Initial/Final values
y0d = Pd_data(:,1);   % Initial demo position
gd = Pd_data(:,end); % Target demo position
y0 = y0d; % set initial position for execution (for simplicity lets leave it the same as the demo)
g = -gd;  % set target position for execution
Tf = Timed(end); % set the time duration of the executed motion

dt = 0.005; % time step for numerical integration

%% Simulate DMP
get_target_fun = @(t) g;
[Time, P_data, dP_data, ddP_data] = simulateModel(DMP_pp(gmp), dt, Tf, y0, 'get_target_fun',get_target_fun);

dmp_classic = DMP_classic(n_dof, 25);
classic_mse = dmp_classic.train(Timed, Pd_data, dPd_data, ddPd_data);
[Time2, P2_data, dP2_data, ddP2_data] = dmp_classic.generate_trajectory(y0, g, Tf, dt);

toc(t_start)

%% Plot results
plot_DMP_comparison_1DoF();


%% ============================================================
%% ============================================================
