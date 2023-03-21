clc;
close all;
clear;

%% =============  includes...  =============
addpath('utils/')
import_gmp_lib();


%% =============  Create demo  =============   
Timed = 0:0.005:3;
Pd0 = [5; 5; 5]; %[0; 0; 0];
Pgd = [5.5; 5.5; 5.5]; %[0.8; 0.9; 1.1];

[Pd_data, dPd_data, ddPd_data] = get5thOrderTraj(Pd0, Pgd, Timed);

n_dof = size(Pd_data, 1);

%% ============= Train DMP  =============
gmp = trainDMP(Timed, Pd_data);


%% =============  DMP simulation  =============
disp('Simulation...');
t_start = tic;

%% Initial/Final values
y0d = Pd_data(:,1);   % Initial demo position
gd = Pd_data(:,end); % Target demo position
y0 = y0d + [0.25; -0.4; 0.1]; % set initial position for execution (for simplicity lets leave it the same as the demo)
g = gd + [-0.7; 0.5; -0.2];  % set target position for execution
Tf = Timed(end); % set the time duration of the executed motion

y0 = [0; 7.0711; 5];
g = [0; 7.7782; 5.5];

% g(1) = y0(1);
g(1) = 0;
y0(1) = 0;


dt = 0.005; % time step for numerical integration

%% Simulate DMP
get_target_fun = @(t) g;
[Time, P_data, dP_data, ddP_data] = simulateModel(gmp, y0, get_target_fun, Tf, dt);
[Time2, P2_data, dP2_data, ddP2_data] = simulateDMP(gmp, y0, get_target_fun, Tf, dt);

toc(t_start)

%% Plot results
ax_font = 13;
x_font = 16;
y_font = 16;
legend_font = 17;

% Plot 3D path
figure; hold on;
plot3(P_data(1,:), P_data(2,:), P_data(3,:), 'LineWidth',2, 'LineStyle','-', 'Color','blue', 'DisplayName','DMP$^{++}$');
plot3(P2_data(1,:), P2_data(2,:), P2_data(3,:), 'LineWidth',2, 'LineStyle',':', 'Color','magenta', 'DisplayName','DMP');
% plot3(Pd_data(1,:), Pd_data(2,:), Pd_data(3,:), 'LineWidth',2, 'LineStyle','--', 'Color','green', 'DisplayName','demo');
% plot3(Pd0(1), Pd0(2), Pd0(3), 'LineStyle','None', 'Marker','o', 'Color',[0.5 1 0.5], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
% plot3(Pgd(1), Pgd(2), Pgd(3), 'LineStyle','None', 'Marker','x', 'Color',[1 0.5 0.5], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
plot3(y0(1), y0(2), y0(3), 'LineStyle','None', 'Marker','o', 'Color',[0 0.85 0], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
plot3(g(1), g(2), g(3), 'LineStyle','None', 'Marker','x', 'Color',[0.85 0 0], 'LineWidth',3, 'MarkerSize',10, 'HandleVisibility','off');
xlabel('X [$m$]', 'interpreter','latex', 'fontsize',15);
ylabel('Y [$m$]', 'interpreter','latex', 'fontsize',15);
zlabel('Z [$m$]', 'interpreter','latex', 'fontsize',15);
legend({}, 'interpreter','latex', 'fontsize',16, 'Box','off');
grid on;
axis equal;
% view(138, 15.7);
view(28.2, 88.5);


% Plot trajectories
for i=1:n_dof
    fig = figure;
    ax_vec = [];
    fig.Position(3:4) = [581 656];
    % Plot Position
    ax = subplot(3,1,1); hold(ax, 'on');
    ax_vec = [ax_vec ax];
    if exist('Pg_data','var') == 1
        plot(Time, Pg_data(i,:), 'LineWidth',1.5, 'LineStyle','-', 'Color',[1, 0.0, 0.0], 'HandleVisibility','off');
    end
    plot(Time, P_data(i,:), 'LineWidth',2.0 , 'Color','blue', 'DisplayName', 'DMP$^{++}$');
    dmp_pl = plot(Time2, P2_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta', 'DisplayName','DMP');
    plot(Timed, Pd_data(i,:), 'LineWidth',2.0, 'LineStyle','-.' , 'Color','green', 'DisplayName', 'demo');
    plot(Time(end), gd(i), 'LineWidth',3, 'Marker','x', 'MarkerSize',12, 'LineStyle','none', 'Color',[1 0.6 0.6], 'DisplayName','$g_d$');
    plot(Time(end), g(i), 'LineWidth',3, 'Marker','x', 'MarkerSize',12, 'LineStyle','none', 'Color','red', 'DisplayName','$g$');
    ax.FontSize = ax_font;
    ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',y_font);
    % title(['DoF $' num2str(i) '$'], 'interpreter','latex', 'fontsize',18);
    legend({}, 'interpreter','latex', 'fontsize',legend_font, 'Position',[0.2251 0.9396 0.5708 0.0425], 'Orientation','horizontal', 'Box','off');
    axis tight;
    xlim([Time(1) Time(end)]);
    hold off;
    
    % Plot Velocity
    ax = subplot(3,1,2); hold(ax, 'on');
    ax_vec = [ax_vec ax];
    plot([Time(1) Time(end)], [0 0], 'LineWidth',1.0, 'LineStyle',':', 'Color',0.4*[1 1 1]);
    plot(Time, dP_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Time2, dP2_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    plot(Timed, dPd_data(i,:), 'LineWidth',2.0, 'LineStyle','-.' , 'Color','green');
    ax.FontSize = ax_font;
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',y_font);
    axis tight;
    xlim([Time(1) Time(end)]);
    hold off;
    
    % Plot Acceleration
    ax = subplot(3,1,3); hold(ax, 'on');
    ax_vec = [ax_vec ax];
    plot([Time(1) Time(end)], [0 0], 'LineWidth',1.0, 'LineStyle',':', 'Color',0.4*[1 1 1]);
    plot(Time, ddP_data(i,:), 'LineWidth',2.0, 'Color','blue');
    plot(Time2, ddP2_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color','magenta');
    plot(Timed, ddPd_data(i,:), 'LineWidth',2.0, 'LineStyle','-.' , 'Color','green');
    ax.FontSize = ax_font;
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',y_font);
    xlabel('time [$s$]', 'interpreter','latex', 'fontsize',x_font);
    axis tight;
    xlim([Time(1) Time(end)]);
    hold off;
    
    linkaxes(ax_vec, 'x');
    for i=1:length(ax_vec)
       ax = ax_vec(i);
       ax.Box = 'on';
       ax.YLim = ax.YLim + 0.07*(ax.YLim(2)-ax.YLim(1))*[-1 1];
    end
    ax = ax_vec(i);
    ax.XLim(2) = ax.XLim(2) + 0.05;
end

%% ============================================================
%% ============================================================


function [y, y_dot, y_ddot] = get5thOrderTraj(y0, yf, Time)

    y0 = y0(:);
    yf = yf(:);
    Time = Time(:)';
    
    T = Time(end);
    t = Time/T;
    
    y = y0 + (yf - y0) * (10*t.^3 - 15*t.^4 + 6*t.^5 );

    if (nargout > 1)
        y_dot = (yf - y0) * (30*t.^2 - 60*t.^3 + 30*t.^4 ) / T;
    end

    if (nargout > 2)
        y_ddot = (yf - y0) * (60*t - 180*t.^2 + 120*t.^3 ) / T^2;
    end

end