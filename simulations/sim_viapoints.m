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

model = dmp_pp;
model_color = [0 0 1];

% model = vmp;
% model_color = [0.85 0.33 0.1];

%% =============  Set init/target  =============
t = 0;
Tf = 6;
tau = Tf;
dt = 0.004;

box_pos = [0.9; 0];
box_h = 0.25;
box_w = 0.1;
box_offset = [0.25; 0.1];
t_box_change = 0.5*Tf;
target_changed = false;
target_change_on = 1*true;
g = box_pos;
g_vp_offsets = [0.24 0.12];
g_viapoints = g + [zeros(size(g_vp_offsets)); g_vp_offsets];

obst_pos = [0.5; 0];
obst_h = 0.35;
obst_w = 0.2;
obst_offset = [0.15; 0.2];
t_obst_change = 0.3*Tf;
obst_changed = false;
obst_change_on = 1*true;
obst_vp_offsets = [ [-obst_w/2, obst_h + 0.05];
%                     [        0,      obst_h + 0.07];
                    [obst_w/2,  obst_h + 0.05]
                  ]';
obst_vp = obst_pos + obst_vp_offsets;
y0 = [0; 0];


%% =============  Set up vizualization environment  =============
fig = figure;
ax = axes(); hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal'); box(ax, 'on');
ax.XLim = [-0.05 1.25];
ax.YLim = [-0.05 0.7];
xlabel('X [$m$]', 'interpreter','latex', 'fontsize',15, 'Parent',ax);
ylabel('Y [$m$]', 'interpreter','latex', 'fontsize',15, 'Parent',ax);
obst = draw_obstacle(obst_pos, obst_h, obst_w, obst_vp, ax);

target_box = draw_target_box(box_pos, box_h, box_w, g_viapoints, ax);

plot(y0(1), y0(2), 'LineStyle','None', 'Marker','o', 'MarkerSize',14, 'LineWidth',3, 'Color','green', 'Parent',ax);

dmp_path = plot(nan, nan, 'LineStyle','-', 'LineWidth',2, 'Color',model_color, 'Parent',ax);
ref_path = plot(nan, nan, 'LineStyle',':', 'LineWidth',2, 'Color','cyan', 'Parent',ax);

plot_every = 4;
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

model.updateViapoints(can_sys.s, obst_vp, 'obst_vp');
model.updateViapoints(can_sys.s, g_viapoints, 'target_vp');

plot_future_path(model, can_sys.s, ax, 'color',[0.5 0.5 1]);
pause

% sv = [0 0.4576 0.5496 1];
% pv = zeros(n_dof, length(sv));
% for j=1:length(sv), pv(:,j) = vmp.getRefPos(sv(j)); end
% plot(pv(1,:), pv(2,:), 'LineStyle','None', 'Marker','o', 'MarkerSize',14, 'LineWidth',3, 'Color','green', 'Parent',ax);
% pause

% model.gmp.W = 0*model.gmp.W;

while (true)
    
    %% data logging
    Time = [Time t];
    Y_data = [Y_data y];
    dY_data = [dY_data y_dot];
    ddY_data = [ddY_data y_ddot];
        
    % simulate obstacle or target changes, e.g. tracked through vision
    if (obst_change_on && t > t_obst_change)
        obst_change_on = false; % turn off, to not enter again
        obst_changed = true;
        delete_object(obst);
        obst_pos = obst_pos + obst_offset;
        obst_vp = obst_pos + obst_vp_offsets;
    end
    
    if (target_change_on && t > t_box_change)
        target_change_on = false; % turn off, to not enter again
        target_changed = true;
        box_pos = box_pos + box_offset;
        g = box_pos;
        g_viapoints = g + [zeros(size(g_vp_offsets)); g_vp_offsets];
    end
    
    %% Update model
    model.update(g, can_sys.s, can_sys.s_dot, y, y_dot);
    if (target_changed)
        target_changed = false; % acknowledged, so disable
        model.updateViapoints(can_sys.s, g_viapoints, 'target_vp');
        % for vizualization:
        delete_object(target_box);
        target_box = draw_target_box(box_pos, box_h, box_w, g_viapoints, ax);
%         pause
        plot_future_path(model, can_sys.s, ax, 'color',[0 0 1 0.8]);
        pause
    end
    if (obst_changed)
        obst_changed = false; % acknowledged, so disable
        model.updateViapoints(can_sys.s, obst_vp, 'obst_vp');
        obst = draw_obstacle(obst_pos, obst_h, obst_w, obst_vp, ax);
%         pause
        plot_future_path(model, can_sys.s, ax, 'color',[0 0 1 0.6]);
        pause
    end

    %% Get reference
    s_ddot = can_sys.getPhaseDDot();

    external_signal = 0; % optionally add some external signal

    % DMP transformation system: 
    y_ddot = model.goal_attractor(y, y_dot, tau) + model.shape_attractor(can_sys.s, can_sys.s_dot, s_ddot, tau) + external_signal;
    
    y_ref = model.getRefPos(can_sys.s);
    dy_ref = model.getRefVel(can_sys.s, can_sys.s_dot);
    ddy_ref = model.getRefAccel(can_sys.s, can_sys.s_dot, s_ddot);
    
    if isequal(class(model),'VMP')
        y_ddot = ddy_ref;
        y_dot = dy_ref;
        y = y_ref;
    end
    
    % vizualization
    dmp_path.XData = [dmp_path.XData y(1)];
    dmp_path.YData = [dmp_path.YData y(2)];
    
%     ref_path.XData = [ref_path.XData y_ref(1)];
%     ref_path.YData = [ref_path.YData y_ref(2)];
    
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

s_data = Time / Time(end);
Yd_data = zeros(n_dof, length(s_data));
dYd_data = zeros(size(Yd_data));
ddYd_data = zeros(size(Yd_data));
for j=1:length(s_data)
    Yd_data(:, j) = gmp.getYd(s_data(j));
    dYd_data(:, j) = gmp.getYdDot(s_data(j), 1/Tf);
    ddYd_data(:, j) = gmp.getYdDDot(s_data(j), 1/Tf, 0);
end
pd_plt = plot(Yd_data(1,:), Yd_data(2,:), 'color',[0.0 0.9 0.0], 'LineStyle','-.', 'LineWidth',2, 'Parent',ax);
scatter(ax, Yd_data(1,end), Yd_data(2,end), 'Marker','x', 'SizeData',170, 'LineWidth',2.5, 'MarkerEdgeColor',[1 0.5 0.5]);

fprintf('Error: pos=%e , vel=%e, accel=%e \n', norm(y - g), norm(y_dot), norm(y_ddot));

dY2_data = dY_data;
ddY2_data = ddY_data;

dY_data = [diff(Y_data, 1, 2) ./ diff(Time) zeros(n_dof,1)];
ddY_data = [diff(dY_data, 1, 2) ./ diff(Time) zeros(n_dof,1)];

% Trajectories
ax = {};
figure;
ax{1} = subplot(2,1,1); hold on;  grid on;
plot(Time, vecnorm(dY_data, 2, 1), 'LineWidth',2.0, 'Color', dmp_path.Color);
plot(Time, vecnorm(dYd_data, 2, 1), 'LineWidth',2.0, 'Color', pd_plt.Color, 'LineStyle',pd_plt.LineStyle);
% plot(Time, vecnorm(dY2_data, 2, 1), 'LineWidth',2.0, 'Color', 'magenta', 'LineStyle','--');
ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',16);
axis tight;
box on;
hold off;

ax{2} = subplot(2,1,2); hold on;  grid on;
plot(Time, vecnorm(ddY_data, 2, 1), 'LineWidth',2.0, 'Color', dmp_path.Color);
plot(Time, vecnorm(ddYd_data, 2, 1), 'LineWidth',2.0, 'Color', pd_plt.Color, 'LineStyle',pd_plt.LineStyle);
% plot(Time, vecnorm(ddY2_data, 2, 1), 'LineWidth',2.0, 'Color', 'magenta', 'LineStyle','--');
ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',16);
xlabel('time [$s$]', 'interpreter','latex', 'fontsize',16);
axis tight;
box on;
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

function plt_handles = draw_target_box(c, h, w, g_viapoints, ax)
    
    g = c;
    c(2) = c(2) - 0.02;
%     h = 0.2;
%     w = 0.1;
    lw = 4;
    color = 0.6*[1 1 1];
    h1 = plot([c(1)-w/2, c(1)+w/2], [c(2) c(2)], 'LineWidth',lw, 'Color',color, 'Parent',ax);
    h2 = plot([c(1)-w/2, c(1)-w/2], [c(2) c(2)+h], 'LineWidth',lw, 'Color',color, 'Parent',ax);
    h3 = plot([c(1)+w/2, c(1)+w/2], [c(2) c(2)+h], 'LineWidth',lw, 'Color',color, 'Parent',ax);
    
    h4 = scatter(ax, g(1), g(2), 'Marker','x', 'SizeData',170, 'LineWidth',2.5, 'MarkerEdgeColor','red');
    
    h5 = scatter(ax, g_viapoints(1,:), g_viapoints(2,:), 'Marker','*', 'SizeData',170, 'LineWidth',1.5, 'MarkerEdgeColor','magenta');
    
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
    
    h2 = scatter(ax, obst_vp(1,:), obst_vp(2,:), 'Marker','*', 'SizeData',170, 'LineWidth',1.5, 'MarkerEdgeColor','magenta');
    
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
