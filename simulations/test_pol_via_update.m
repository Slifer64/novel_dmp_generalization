clc;
close all;
clear;


% s_via = [0    0.4407    0.5165    0.7458    0.7458    0.7458    1.0000];
% 
% y_via = [0    0.4000    0.6000    0.9000    0.9000    0.9000    0.9000;
%          0    0.4000    0.4000    0.1800    0.1200    0.0600         0];
%      
s_via = [0    0.4407    0.5165    0.7458    1.0000];

y_via = [0    0.4000    0.6000    0.9000    0.9000;
         0    0.4000    0.4000    0.1800         0];
     
n_dof = size(y_via, 1);
mp = FifthOrderMP(n_dof);
 
 

fig = figure;
ax = axes(); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
mp_plt = plot(nan, nan, 'LineStyle','-', 'color',[0 0 1], 'LineWidth',2);
scatter(ax, y_via(1,:), y_via(2,:), 'SizeData',150, 'LineWidth',2, 'Marker','*', 'MarkerEdgeColor',[1 0 1]);
ax.XLim = [-0.05 1];
ax.YLim = [-0.05 0.45];


dt = 0.005;
Tf = 5;
t = 0;
s = 0;
s_dot = 1/Tf;
s_ddot = 0;

O_ndof = zeros(n_dof, 1);

nan_vec = nan(n_dof, 1);

y0 = y_via(:,1);
g = y_via(:, end);

mp.update(0.0, y0, O_ndof, O_ndof, 1.0, g, O_ndof, O_ndof, s_dot, s_ddot);

Time = [];
P_data = [];
dP_data = [];
ddP_data = [];

while true
    
    %% find the two active via-points
    i1 = 1;
    i2 = length(s_via);
    for j=2:length(s_via)
       if (s_via(j) > s)
           i1 = j-1;
           i2 = j;
           break;
       end
    end

    s1 = s_via(i1);
    s2 = s_via(i2);
    
%     yv1 = y_via(:, i1);
    s1 = s;
    if (abs(s1 - s2) < 5e-3), s1 = s2 - 5e-3; end
    yv1 = mp.getRefPos(s1);
    yv1_dot = mp.getRefVel(s1, s_dot); %nan_vec;
    yv1_ddot = nan*mp.getRefAccel(s1, s_dot, s_ddot); %nan_vec;
    yv2 = y_via(:, i2);
    yv2_dot = nan*mp.getRefVel(s2, s_dot); %nan_vec;
    yv2_ddot = nan*mp.getRefAccel(s2, s_dot, s_ddot); %nan_vec;
    
    mp.update(s1, yv1, yv1_dot, yv1_ddot, s2, yv2, yv2_dot, yv2_ddot, s_dot, s_ddot);
    
%     s1 = s;
%     if (abs(s1 - 1.0) < 5e-3), s1 = 1.0 - 5e-3; end
%     yv1 = mp.getRefPos(s1);
%     yv1_dot = mp.getRefVel(s1, s_dot);
%     yv1_ddot = mp.getRefAccel(s1, s_dot, s_ddot);
%     mp.update(s1, yv1, yv1_dot, yv1_ddot, 1.0, g, O_ndof, O_ndof, s_dot, s_ddot);
    
%     yv1_hat = mp.getRefPos(s1);
%     yv1_hat_dot = mp.getRefVel(s1, s_dot);
%     yv1_hat_ddot = mp.getRefAccel(s1, s_dot, s_ddot);
    
%     s
%     err1 = norm(yv1  -yv1_hat)
%     dot_err1 = norm(yv1_dot  -yv1_hat_dot)
%     ddot_err1 = norm(yv1_ddot  -yv1_hat_ddot)
%     pause
    
    p = mp.getRefPos(s);
    p_dot = mp.getRefVel(s, s_dot);
    p_ddot = mp.getRefAccel(s, s_dot, s_ddot);
    
    %% data logging
    Time = [Time t];
    P_data = [P_data p];
    dP_data = [dP_data p_dot];
    ddP_data = [ddP_data p_ddot];
    
    %% vizualization
    mp_plt.XData = [mp_plt.XData p(1)];
    mp_plt.YData = [mp_plt.YData p(2)];
    drawnow();
    pause(0.001);
    
    %% stopping criteria
    if s >= 1 && norm(p - g) < 1e-3
        break;
    end
    
    if t > 1.4*Tf
        warning('Time limit exceeded...');
        break;
    end
    
    
    %% numerical integration
    t = t + dt;
    s = s + s_dot*dt;
    
end

% Trajectories
ax = {};
figure;
ax{1} = subplot(2,1,1); hold on;  grid on;
plot(Time, vecnorm(dP_data, 2, 1), 'LineWidth',2.0, 'Color', 'green');
ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',15);
axis tight;
box on;
hold off;

ax{2} = subplot(2,1,2); hold on;  grid on;
plot(Time, vecnorm(ddP_data, 2, 1), 'LineWidth',2.0, 'Color', 'red');
ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',15);
axis tight;
box on;
hold off;

for j=1:length(ax)
    ax{j}.XLim = ax{j}.XLim + 0.01*(ax{j}.XLim(2) - ax{j}.XLim(1)) * [0 1];
    ax{j}.YLim = ax{j}.YLim + 0.02*(ax{j}.YLim(2) - ax{j}.YLim(1)) * [-1 1];
end