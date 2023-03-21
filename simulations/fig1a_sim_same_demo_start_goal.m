clc;
close all;
clear;

%% =============  includes...  =============
addpath('utils/')
import_gmp_lib();


%% =============  Load demo data  =============                           
load('data/demo_same_start_goal.mat', 'Timed', 'Pd_data');

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
g = gd - 0.1;  % set target position for execution
Tf = Timed(end); % set the time duration of the executed motion

dt = 0.005; % time step for numerical integration

%% Simulate DMP
get_target_fun = @(t) g;
[Time, P_data, dP_data, ddP_data] = simulateModel(DMP_pp(gmp), dt, Tf, y0, 'get_target_fun',get_target_fun);
[Time2, P2_data, dP2_data, ddP2_data] = simulateModel(DMP_classic(gmp), dt, Tf, y0, 'get_target_fun',get_target_fun);

toc(t_start)

%% Plot results

% Scale DMP trajectory for visualization purposes
s = 400;
dP2_data = dP2_data/s;
ddP2_data = ddP2_data/s;
P2_data = P2_data/s;

plot_DMP_comparison_1DoF();

dmp_pl.DisplayName = ['DMP/' num2str(s)];

%% ============================================================
%% ============================================================
