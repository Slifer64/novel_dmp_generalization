clc;
% close all;
clear;

%% =============  includes...  =============
addpath('utils/')
import_gmp_lib();

%% =============  Load train data  =============
load('data/vp_demo.mat', 'sd_data', 'Pd_data');

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
dt = 0.005;

box_pos = [0.9; 0];
g = box_pos;
g_vp_offsets = [0.18 0.12 0.06];
g_viapoints = g + [zeros(size(g_vp_offsets)); g_vp_offsets];

obst_pos = [0.5; 0];
obst_h = 0.35;
obst_w = 0.2;
obst_vp_offsets = [-obst_w/2, obst_w/2; obst_h + 0.05, obst_h + 0.05];
obst_vp = obst_pos + obst_vp_offsets;
y0 = [0; 0];


%% =============  Set up vizualization environment  =============
fig = figure;
ax = axes(); hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
ax.XLim = [-0.05 1.2];
ax.YLim = [-0.05 0.6];
obst = draw_obstacle(obst_pos, obst_h, obst_w, obst_vp, ax);

box = draw_target_box(box_pos, g_viapoints, ax);

plot(y0(1), y0(2), 'LineStyle','None', 'Marker','o', 'MarkerSize',14, 'LineWidth',3, 'Color','green', 'Parent',ax);

dmp_path = plot(nan, nan, 'LineStyle','-', 'LineWidth',2, 'Color','blue', 'Parent',ax);

plot_every = 2;
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

target_changed = false;
target_change_on = 0*true;

obst_changed = false;
obst_change_on = 0*true;

plot_future_path(model, can_sys.s, ax, 'color',[0 0 1 0.2]);
% pause

model.updateViapoints(can_sys.s, obst_vp, 'obst_vp');
model.updateViapoints(can_sys.s, g_viapoints, 'target_vp');
% model.update(g, can_sys.s, can_sys.s_dot, y, y_dot);

plot_future_path(model, can_sys.s, ax, 'color',[0 0 1 0.4]);
% pause

% sv = [0 0.4576 0.5496 1];
% pv = zeros(n_dof, length(sv));
% for j=1:length(sv), pv(:,j) = vmp.getRefPos(sv(j)); end
% plot(pv(1,:), pv(2,:), 'LineStyle','None', 'Marker','o', 'MarkerSize',14, 'LineWidth',3, 'Color','green', 'Parent',ax);
% pause

while (true)
    
    %% data logging
    Time = [Time t];
    Y_data = [Y_data y];
    dY_data = [dY_data y_dot];
    ddY_data = [ddY_data y_ddot];
        
    % simulate obstacle or target changes, e.g. tracked through vision
    if (obst_change_on && t > 0.3*Tf)
        obst_change_on = false; % turn off, to not enter again
        obst_changed = true;
        delete_object(obst);
        obst_pos = obst_pos + [0.1; 0.15];
        obst_vp = obst_pos + obst_vp_offsets;
    end
    
    if (target_change_on && t > 0.5*Tf)
        target_change_on = false; % turn off, to not enter again
        target_changed = true;
        box_pos = box_pos + [0.2; 0.15];
        g = box_pos;
        g_viapoints = g + [zeros(size(g_vp_offsets)); g_vp_offsets];
    end
    
    %% Update model
    model.update(g, can_sys.s, can_sys.s_dot, y, y_dot);
    if (target_changed)
        target_changed = false; % acknowledged, so disable
        model.updateViapoints(can_sys.s, g_viapoints, 'target_vp');
        % for vizualization:
        delete_object(box);
        box = draw_target_box(box_pos, g_viapoints, ax);
%         pause
        plot_future_path(model, can_sys.s, ax, 'color',[0 0 1 0.8]);
%         pause
    end
    if (obst_changed)
        obst_changed = false; % acknowledged, so disable
        model.updateViapoints(can_sys.s, obst_vp, 'obst_vp');
        obst = draw_obstacle(obst_pos, obst_h, obst_w, obst_vp, ax);
%         pause
        plot_future_path(model, can_sys.s, ax, 'color',[0 0 1 0.6]);
%         pause
    end

    %% Get reference
    s_ddot = can_sys.getPhaseDDot();

    external_signal = 0; % optionally add some external signal

    % DMP transformation system: 
    y_ddot = model.goal_attractor(y, y_dot, tau) + model.shape_attractor(can_sys.s, can_sys.s_dot, s_ddot, tau) + external_signal;
    
%     y = model.getRefPos(can_sys.s);
    
    % vizualization
    dmp_path.XData = [dmp_path.XData y(1)];
    dmp_path.YData = [dmp_path.YData y(2)];
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
    y = y + y_dot*dt;
    y_dot = y_dot + y_ddot*dt;
    
end

fprintf('Error: pos=%e , vel=%e, accel=%e \n', norm(y - g), norm(y_dot), norm(y_ddot));

% Trajectories
ax = {};
figure;
ax{1} = subplot(2,1,1); hold on;  grid on;
plot(Time, vecnorm(dY_data, 2, 1), 'LineWidth',2.0, 'Color', 'green');
ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',15);
axis tight;
hold off;

ax{2} = subplot(2,1,2); hold on;  grid on;
plot(Time, vecnorm(ddY_data, 2, 1), 'LineWidth',2.0, 'Color', 'red');
ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',15);
axis tight;
hold off;

for j=1:length(ax)
    ax{j}.XLim = ax{j}.XLim + 0.01*(ax{j}.XLim(2) - ax{j}.XLim(1)) * [0 1];
    ax{j}.YLim = ax{j}.YLim + 0.02*(ax{j}.YLim(2) - ax{j}.YLim(1)) * [-1 1];
end


%% ========================================================

function delete_object(obj_handles)
    
    for i=1:length(obj_handles)
       delete(obj_handles(i)); 
    end
    
end

function h = plot_future_path(model, s, ax, varargin)
    
    n_points = 40;
    s_data = linspace(s, 1, n_points);
    P_data = zeros(2, n_points);
    for j=1:n_points
       P_data(:, j) = model.getRefPos(s_data(j)); 
    end
    
    h = plot(P_data(1,:), P_data(2,:), 'LineStyle',':', 'LineWidth',2, 'Color',[0 1 0 0.6], 'Parent',ax);
    set(h, varargin{:});
end

function plt_handles = draw_target_box(c, g_viapoints, ax)
    
    g = c;
    c(2) = c(2) - 0.02;
    h = 0.2;
    w = 0.1;
    lw = 4;
    color = [0.85, 0.33, 0.1];
    h1 = plot([c(1)-w/2, c(1)+w/2], [c(2) c(2)], 'LineWidth',lw, 'Color',color, 'Parent',ax);
    h2 = plot([c(1)-w/2, c(1)-w/2], [c(2) c(2)+h], 'LineWidth',lw, 'Color',color, 'Parent',ax);
    h3 = plot([c(1)+w/2, c(1)+w/2], [c(2) c(2)+h], 'LineWidth',lw, 'Color',color, 'Parent',ax);
    
    h4 = plot(g(1), g(2), 'LineStyle','None', 'Marker','x', 'MarkerSize',14, 'LineWidth',3, 'Color','red', 'Parent',ax);
    
    h5 = plot(g_viapoints(1,:), g_viapoints(2,:), 'LineStyle','None', 'Marker','*', 'MarkerSize',14, 'LineWidth',2, 'Color','magenta', 'Parent',ax);
    
    plt_handles = [h1 h2 h3 h4 h5];
end


function obst_handles = draw_obstacle(center, height, width, obst_vp, ax)

    center_x = center(1);
    center_y = center(2);

    % Define the corner points of the rectangle
    x_vals = [center_x - width/2, center_x + width/2, center_x + width/2, center_x - width/2];
    y_vals = [center_y, center_y, center_y + height, center_y + height];

    % Plot the filled rectangle
    h1 = fill(x_vals, y_vals, 'r');
    set(h1, 'LineStyle','none', 'FaceColor',0.5*[1 1 1], 'FaceAlpha',0.6, 'Parent',ax);
    
    h2 = plot(obst_vp(1,:), obst_vp(2,:), 'LineStyle','None', 'Marker','*', 'MarkerSize',14, 'LineWidth',2, 'Color','magenta', 'Parent',ax);
    
    obst_handles = [h1 h2];
end


function h = draw_ellipsoid(center, lambda_x, lambda_y, theta, varargin)

    center_x = center(1);           % x-coordinate of center point
    center_y = center(2);           % y-coordinate of center point
    a = lambda_x;                  % length of semi-axis along x direction
    b = lambda_y;                  % length of semi-axis along y direction
    theta = theta * pi/180; % angle of rotation (in radians) of the ellipsoid

    % Define the number of points to use in the ellipse approximation
    num_points = 100;

    % Create arrays of x and y values for the ellipse
    phi = linspace(0, 2*pi, num_points);     % angle values
    x_vals = center_x + a*cos(phi)*cos(theta) - b*sin(phi)*sin(theta);   % x values
    y_vals = center_y + a*cos(phi)*sin(theta) + b*sin(phi)*cos(theta);   % y values

    % Plot the filled ellipsoid
    h = fill(x_vals, y_vals, 'r');
    set(h, 'LineStyle','none', varargin{:});

end
