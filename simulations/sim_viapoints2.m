clc;
close all;
clear;

%% =============  includes...  =============
addpath('utils/')
import_gmp_lib();

%% =============  Load train data  =============
load('data/vp_demo2.mat', 'sd_data', 'Pd_data');

n_dof = size(Pd_data,1);
y0d = Pd_data(:,1);
gd = Pd_data(:,end);


%% =============  Train DMP++  =============
gmp = GMP(n_dof, 25);
gmp_mse = gmp.train('LS', sd_data, Pd_data)
dmp_pp = DMP_pp(gmp);

vmp = VMP(gmp);

% model = dmp_pp;
model = vmp;

%% =============  Set init/target  =============
t = 0;
Tf = 6;
tau = Tf;
dt = 0.002;

g = gd;
y0 = y0d;

%% =============  Set up vizualization environment  =============
fig = figure;
ax = axes(); hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
% ax.XLim = [-0.05 1.2];
% ax.YLim = [-0.05 0.6];

plot(y0(1), y0(2), 'LineStyle','None', 'Marker','o', 'MarkerSize',14, 'LineWidth',3, 'Color','green', 'Parent',ax);
plot(g(1), g(2), 'LineStyle','None', 'Marker','x', 'MarkerSize',14, 'LineWidth',3, 'Color','red', 'Parent',ax);

dmp_path = plot(nan, nan, 'LineStyle','-', 'LineWidth',2, 'Color','blue', 'Parent',ax);
ref_path = plot(nan, nan, 'LineStyle',':', 'LineWidth',2, 'Color','cyan', 'Parent',ax);

plot_every = 20;
plot_count = 0;

%% =============  Run simulation  =============
y = y0; % position
y_dot = zeros(size(y)); % velocity
y_ddot = zeros(size(y)); % acceleration

Time = [];
Y_data = [];
dY_data = [];
ddY_data = [];

can_sys = CanonicalSystem(Tf);

model.rf = [1e-11, 1e-7, 1e-7];
model.rv = 1e-8;
model.r1 = [1e-7, 1e-7, 1e-6];

model.init(can_sys.s, y0, g, Tf);
model.K = 300; % DMP stiffness
model.D = 2*sqrt(model.K + 10); % DMP damping

plot_future_path(model, can_sys.s, ax, 'color',[0 0 1 0.2]);
% pause

via_points = [-0.75;
              0.08];
s_v = model.updateViapoints(can_sys.s, via_points, 'vp_1');
plot_future_path(model, can_sys.s, ax, 'color',[0 0 1 0.4]);
% pause

plot(via_points(1, :), via_points(2, :), 'LineStyle','None', 'Marker','*', 'MarkerSize',14, 'LineWidth',2, 'Color','magenta', 'Parent',ax);

% model.gmp.W = 0*model.gmp.W;

while (true)
    
    %% data logging
    Time = [Time t];
    Y_data = [Y_data y];
    dY_data = [dY_data y_dot];
    ddY_data = [ddY_data y_ddot];

    %% Get reference
    s_ddot = can_sys.getPhaseDDot();
    
    %% Update model
    model.update(g, can_sys.s, can_sys.s_dot, y, y_dot);

    y_ref = model.getRefPos(can_sys.s);
    dy_ref = model.getRefVel(can_sys.s, can_sys.s_dot);
    ddy_ref = model.getRefAccel(can_sys.s, can_sys.s_dot, s_ddot);
    
    y_ddot = ddy_ref;
    y_dot = dy_ref;
    y = y_ref;
    
    % vizualization
    dmp_path.XData = [dmp_path.XData y(1)];
    dmp_path.YData = [dmp_path.YData y(2)];
    
    ref_path.XData = [ref_path.XData y_ref(1)];
    ref_path.YData = [ref_path.YData y_ref(2)];
    
    plot_count = plot_count + 1;
    if (plot_count == plot_every)
        plot_count = 0;
        drawnow();
        pause(0.001);
    end
    
    %% Stopping criteria
    if (t>=(1.0*Tf+0.05) && norm(y-g)<1e-3 && norm(y_dot)<5e-3)
        break;
    end

    if (t >= 1.4*Tf)
        warning('Time limit exceeded!');
        break;
    end

    %% Numerical integration
    can_sys.integrate(t, t+dt);
    t = t + dt;
%     y = y + y_dot*dt;
%     y_dot = y_dot + y_ddot*dt;
    
end

fprintf('Error: pos=%e , vel=%e, accel=%e \n', norm(y - g), norm(y_dot), norm(y_ddot));


dY2_data = [diff(Y_data, 1, 2) ./ diff(Time) zeros(n_dof,1)];
dY2_norm = vecnorm(dY2_data, 2, 1);

% Trajectories
ax = {};
figure;
ax{1} = subplot(2,1,1); hold on;  grid on;
plot(Time, vecnorm(dY_data, 2, 1), 'LineWidth',2.0, 'Color', 'green');
plot(Time, dY2_norm, 'LineWidth',2.0, 'Color', 'magenta', 'LineStyle',':');
ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',15);
axis tight;
box on;
hold off;

ddY2_norm = [vecnorm( diff(dY2_data, 1, 2) ./ diff(Time), 2, 1) 0];
ddY3_norm = [vecnorm( diff(dY_data, 1, 2) ./ diff(Time), 2, 1) 0];


ax{2} = subplot(2,1,2); hold on;  grid on;
plot(Time, vecnorm(ddY_data, 2, 1), 'LineWidth',2.0, 'Color', 'red');
plot(Time, ddY2_norm, 'LineWidth',2.0, 'Color', 'magenta', 'LineStyle',':');
plot(Time, ddY3_norm, 'LineWidth',2.0, 'Color', 'cyan', 'LineStyle','--');
ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',15);
box on;
axis tight;
ax{2}.YLim = [0 min([max(ax{2}.YLim) 3.0])];

for j=1:length(ax)
    ax{j}.XLim = ax{j}.XLim + 0.01*(ax{j}.XLim(2) - ax{j}.XLim(1)) * [0 1];
    ax{j}.YLim = ax{j}.YLim + 0.02*(ax{j}.YLim(2) - ax{j}.YLim(1)) * [-1 1];
end


%% ========================================================


function h = plot_future_path(model, s, ax, varargin)
    
    n_points = 70;
    s_data = linspace(s, 1, n_points);
    P_data = zeros(2, n_points);
    for j=1:n_points
       P_data(:, j) = model.getRefPos(s_data(j)); 
    end
    
    h = plot(P_data(1,:), P_data(2,:), 'LineStyle',':', 'LineWidth',2, 'Color',[0 1 0 0.6], 'Parent',ax);
    set(h, varargin{:});
end
