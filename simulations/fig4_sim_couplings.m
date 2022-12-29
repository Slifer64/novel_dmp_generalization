clc;
close all;
clear;

%% =============  includes...  =============
addpath('utils/')
import_gmp_lib();

%% ========= Load demo data ===========
load('data/couplings_demo_data.mat');

Yd_data = Pd_data;

Timed = sd_data * 8;
y0d = Yd_data(:, 1);
ygd = Yd_data(:, end);

n_dof = size(Yd_data, 1);


%% ========= Define obstacles ===========
obstacle = {};
obstacle = [obstacle {EllipsoidObstacle('center',[-0.02; 0.33], 'angle',-20, 'lambda_x',0.11, 'lambda_y',0.07, 'd0',1.0, 'color',0.3*[1 1 1])}];
obstacle = [obstacle {EllipsoidObstacle('center',[-0.3; 0.53], 'angle',75, 'lambda_x',0.15, 'lambda_y',0.06, 'd0',1.0, 'color',0.3*[1 1 1])}];

obstacle = [obstacle {PlaneObstacle('point',[0.16; 0.73], 'normal_angle',-120, 'color', 0.3*[1 1 1])}];
obstacle = [obstacle {PlaneObstacle('point',[0.21; 0.64], 'normal_angle',180, 'color', 0.3*[1 1 1])}];

k_rep = 0.1; % repulsive force gain

adapt_to_ext_signals = true;

%% ========= Train model ===========
import_gmp_lib();
gmp = GMP(n_dof, 25, 1.5);
t_start = tic;
sd_data = Timed/Timed(end);
offline_train_mse = gmp.train('LS', sd_data, Yd_data);
offline_train_mse
toc(t_start)


%% ========== Simulation ============
%% set initial values
y0 = y0d;
yg = ygd + [0.07; 0.15];
Tf = Timed(end);

y = y0;
y_dot = zeros(n_dof,1);
y_ddot = zeros(n_dof,1);
O_ndof = zeros(n_dof, 1);

t = 0.0;
dt = 0.001;

t_end = Tf;
tau = t_end;

% phase variable, from 0 to 1
s = 0.0;
s_dot = 1/tau;
s_ddot = 0; % since s_dot is constant here

f_rep = zeros(n_dof,1);
iters = 0;

% data to log
Time = [];
Y_data = [];
dY_data = [];
ddY_data = [];

model = DMP_pp(gmp);
model.init(s, y0, yg, Tf);

model.setAdaptToRobot(adapt_to_ext_signals);
model.r1(1) = 1e-4;
model.r1(2) = 1e-1;
model.r1(3) = 1e-1;

% generate nominal trajectory
s_data = linspace(0,1, 100);
Y0_data = zeros(n_dof, length(s_data));
for j=1:length(s_data), Y0_data(:,j) = model.getRefPos(s_data(j)); end

% Initialize online plot
plot_every = 50;
plot_count = 0;
fig = figure; fig.Position(3:4) = [582 769]; 

ax = axes('Position',[-2 -2 0.01 0.01]); hold(ax, 'on');
y_nom_pl = plot(Y0_data(1,:), Y0_data(2,:), 'LineStyle',':', 'LineWidth',2, 'Color',[0 0.8 0], 'DisplayName','nominal');
y0_pl = plot(y0d(1), y0d(2), 'LineStyle','None', 'LineWidth',3, 'MarkerSize',12, 'Marker','o', 'Color','green', 'DisplayName','$y_0$');
g_pl = plot(yg(1), yg(2), 'LineStyle','None', 'LineWidth',3, 'MarkerSize',12, 'Marker','x', 'Color','red', 'DisplayName','$g$');
lg2 = legend({}, 'interpreter','latex', 'fontsize',14, 'Orientation','horizontal', 'Box','off');

ax = axes('Parent',fig, 'FontSize',13); hold(ax, 'on');
ax.Box = 'on';
y_nom_pl = copyobj(y_nom_pl, ax); y_nom_pl.HandleVisibility='off';
y_pl = plot(nan, nan, 'LineStyle','-', 'LineWidth',2, 'Color','blue', 'DisplayName','DMP$^{++}$');
y_rev_pl = plot(nan, nan, 'LineStyle','-.', 'LineWidth',2, 'Color','cyan', 'DisplayName','rev-DMP$^{++}$');
% y_final_pl = plot(nan, nan, 'LineStyle','--', 'LineWidth',2, 'Color','magenta', 'DisplayName','final DMP$^{++}$');
y0_pl = copyobj(y0_pl, ax); y0_pl.HandleVisibility='off';
g_pl = copyobj(g_pl, ax); g_pl.HandleVisibility='off';
axis(ax, 'equal');
xlabel('X [m]', 'interpreter','latex', 'fontsize',16);
ylabel('Y [m]', 'interpreter','latex', 'fontsize',16);
ax.XLim = [-0.3770 0.2275];
ax.YLim = [-0.0415 0.8401];
ax.XLim = ax.XLim + 0.05*(ax.XLim(2) - ax.XLim(1))*[-1 1];
ax.YLim = ax.YLim + 0.05*(ax.YLim(2) - ax.YLim(1))*[-1 1];
for i=1:length(obstacle), obstacle{i}.plot(ax); end
lg = legend({}, 'interpreter','latex', 'fontsize',16, 'Orientation','horizontal', 'Box','off');
lg.Position(1) = ax.Position(1) + (ax.Position(3)-lg.Position(3))/2;
lg.Position(2) = ax.Position(2) + ax.Position(4) + 0.005;
lg2.Position(1) = lg.Position(1);
lg2.Position(2) = lg.Position(2) + 0.035;

fig = figure; fig.Position(3:4)=[582 404];
ax = axes('Parent',fig); 
hold(ax, 'on');
grid(ax, 'on');
ax.Box = 'on';
ax.FontSize = 12;
frep_pl = plot(nan, nan, 'LineStyle','-', 'LineWidth',2, 'Color','red', 'DisplayName','$||f_{rep}||$');
ax.XLim = [0 2*Tf];
ax.YLim = [-0.01 20];
plot([Tf Tf], ax.YLim, 'LineWidth',2, 'LineStyle',':', 'Color',0.7*[1 1 1], 'HandleVisibility','off');
xlabel('time [$s$]', 'interpreter','latex', 'fontsize',16);
ylabel('[$N$]', 'interpreter','latex', 'fontsize',17);
legend({}, 'interpreter','latex', 'fontsize',19, 'Box','off', 'Position',[0.6956 0.8218 0.1882 0.0949]);

frw_txt = annotation(fig, 'textbox', 'String','Forward', 'Color',y_pl.Color, 'FontSize',18, 'LineStyle','none', 'Position',[0.2204 0.9312 0.1937 0.0577]);
frw_txt.Position(1) = ax.Position(1) + 0.25*ax.Position(3) - 0.5*frw_txt.Position(3);
frw_txt.Position(2) = ax.Position(2)+ ax.Position(4) + 0.02;

rev_txt = annotation(fig, 'textbox', 'String','Reverse', 'Color',y_rev_pl.Color, 'FontSize',18, 'LineStyle','none', 'Position',[0.6221 0.9312 0.1937 0.0601]);
rev_txt.Position(1) = ax.Position(1) + 0.75*ax.Position(3) - 0.5*rev_txt.Position(3);
rev_txt.Position(2) = ax.Position(2)+ ax.Position(4) + 0.02;


%% simulation loop
while (true)

    %% data logging
    Time = [Time t];
    Y_data = [Y_data y];
    dY_data = [dY_data y_dot];
    ddY_data = [ddY_data y_ddot];
    
    %% plot
    y_pl.XData = [y_pl.XData y(1)];
    y_pl.YData = [y_pl.YData y(2)];
    frep_pl.XData = [frep_pl.XData t];
    frep_pl.YData = [frep_pl.YData norm(f_rep)];
    
    plot_count = plot_count + 1;
    if (plot_count == plot_every)
        plot_count = 0;
        drawnow();
        pause(0.001);
    end

    %% Update DMP_pp
    model.update(yg, s, s_dot, y, y_dot);

    %% Get reference
    y_s = model.getRefPos(s);
    dy_s = model.getRefVel(s, s_dot);
    ddy_s = model.getRefAccel(s, s_dot, s_ddot);

    K = 300; % set the DMP stiffness
    D = 60; % set the DMP damping

    f_rep = zeros(n_dof, 1); % repulsive force
    for i=1:length(obstacle)
        f_rep = f_rep + obstacle{i}.repulsive_force(y, k_rep);
    end

    % GMP equation (equivalent to DMP) 
    y_ddot = ddy_s + D*(dy_s - y_dot) + K*(y_s - y) + f_rep;

    %% Stopping criteria
    if (t>=1.0*t_end && norm(y-yg)<1e-3 && norm(y_dot)<5e-3)
        break;
    end

    if (t >= 1.4*t_end)
        warning('Time limit exceeded!');
        break; 
    end

    %% Numerical integration
    iters = iters + 1;
    t = t + dt;
    s = s + s_dot*dt;
    s_dot = s_dot + s_ddot*dt;
    y = y + y_dot*dt;
    y_dot = y_dot + y_ddot*dt;

end

fprintf('Error: pos=%e , vel=%e, accel=%e \n', norm(y - yg), norm(y_dot), norm(y_ddot));


%% run reverse
t = 0;
s_dot = -s_dot;
model.setAdaptToRobot(adapt_to_ext_signals);
model.r1(1) = 1e-5;
enable_rep_force = 1;

while (true)
    
    %% plot
    y_rev_pl.XData = [y_rev_pl.XData y(1)];
    y_rev_pl.YData = [y_rev_pl.YData y(2)];
    frep_pl.XData = [frep_pl.XData t+Tf];
    frep_pl.YData = [frep_pl.YData norm(f_rep)];
    
    plot_count = plot_count + 1;
    if (plot_count == plot_every)
        plot_count = 0;
        drawnow();
        pause(0.001);
    end
    
    %% Update DMP_pp
    model.update(yg, s, s_dot, y, y_dot);
    
    %% Get reference
    y_s = model.getRefPos(s);
    dy_s = model.getRefVel(s, s_dot);
    ddy_s = model.getRefAccel(s, s_dot, s_ddot);

    K = 300; % set the DMP stiffness
    D = 60; % set the DMP damping
    
    f_rep = 0; % optionally add some external signal
    for i=1:length(obstacle)
        f_rep = f_rep + obstacle{i}.repulsive_force(y, k_rep);
    end
    
    % Track it using a 2nd order dynamical system. This is actually the DMP. 
    y_ddot = ddy_s + D*(dy_s - y_dot) + K*(y_s - y) + enable_rep_force*f_rep;
    
    %% Stopping criteria
    if (t>=1.0*t_end && norm(y-y0)<1e-3 && norm(y_dot)<5e-3)
        break;
    end
    
    %% Numerical integration
    t = t + dt;
    s = s + s_dot*dt;
    s_dot = s_dot + s_ddot*dt;
    y = y + y_dot*dt;
    y_dot = y_dot + y_ddot*dt;
    
end

fprintf('Error: pos=%e , vel=%e, accel=%e \n', norm(y - y0), norm(y_dot), norm(y_ddot));

    
s_data = linspace(0, 1, 100);
y_data = zeros(n_dof, length(s_data));
for j=1:length(s_data)
    y_data(:,j) = model.getRefPos(s_data(j));
end

% y_final_pl.XData = y_data(1,:);
% y_final_pl.YData = y_data(2,:);

%% ========= Plot results ===========

figure;
ax = axes();
hold(ax, 'on');

y_adapt_pl = plot(nan, nan, 'LineWidth',2, 'LineStyle','-.', 'Color','magenta', 'Parent',ax);

y_adapt_pl.XData = y_data(1,:);
y_adapt_pl.YData = y_data(2,:);

y_nom_pl = copyobj(y_nom_pl, ax);

axis(ax, 'tight');
axis(ax, 'equal');


%% ###################################







