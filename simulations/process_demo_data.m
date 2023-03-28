clc;
close all;
clear;

load('utils/data1.mat', 'sd_data', 'Pd_data');
s1_data = sd_data;
P1_data = Pd_data;

save('data/vp_demo2.mat', 'sd_data', 'Pd_data');
return


load('utils/data2.mat', 'sd_data', 'Pd_data');
s2_data = sd_data;
P2_data = Pd_data;

x_data = interp1(s2_data, P2_data(1,:), s1_data);

sd_data = s1_data;
Pd_data = [x_data; P1_data];

Tf = 4.0;
n_dof = size(Pd_data, 1);
Timed = sd_data * Tf;
dTime = diff(Timed);
dPd_data = [diff(Pd_data, 1, 2) ./ dTime, zeros(n_dof, 1)];
ddPd_data = [diff(dPd_data, 1, 2) ./ dTime, zeros(n_dof, 1)];

vel_norm = vecnorm(dPd_data, 2, 1);
accel_norm = vecnorm(ddPd_data, 2, 1);

figure;
plot3(Pd_data(1,:), Pd_data(2,:), Pd_data(3,:), 'LineWidth',2, 'color','blue');
ylabel('2D path');
axis equal;
grid on;

figure;
subplot(2,1,1);
plot(Timed, vel_norm, 'LineWidth',2, 'color','magenta');
ylabel('vel norm');
subplot(2,1,2);
plot(Timed, accel_norm, 'LineWidth',2, 'color','red');
ylabel('accel norm');

dist = norm(Pd_data(:,end) - Pd_data(:,1))
ampl = max(Pd_data(3,:)) - min(Pd_data(3,:))

save('data/comp_scalings_demo2.mat', 'Timed', 'Pd_data', 'dPd_data', 'ddPd_data');
