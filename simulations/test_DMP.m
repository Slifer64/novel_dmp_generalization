clc;
close all;
clear;

%% =============  includes...  =============
addpath('utils/')
import_gmp_lib();

%% =============  Load train data  =============
load('data/comp_scalings_data.mat', 'Timed', 'Pd_data', 'dPd_data', 'ddPd_data');

Ts = Timed(2) - Timed(1);

%% =============  Create/Train GMP  =============
n_dof = size(Pd_data, 1); % number of DoFs
n_kernels = 25;

dmp_classic = DMP_classic(n_dof, n_kernels);
classic_mse = dmp_classic.train(Timed, Pd_data, dPd_data, ddPd_data)

dmp_rot = DMP_rot(n_dof, n_kernels);
rot_mse = dmp_rot.train(Timed, Pd_data, dPd_data, ddPd_data)

dmp_bio = DMP_bio(n_dof, n_kernels);
bio_mse = dmp_bio.train(Timed, Pd_data, dPd_data, ddPd_data)

dmp_bio_plus = DMP_bio_plus(n_dof, n_kernels);
bio_plus_mse = dmp_bio_plus.train(Timed, Pd_data, dPd_data, ddPd_data)

gmp = GMP(n_dof, n_kernels, 1.5);
gmp_mse = gmp.train('LS', Timed/Timed(end), Pd_data)
dmp_pp = DMP_pp(gmp);

y0d = Pd_data(:, 1);
gd = Pd_data(:, end);

y0 = y0d;
g = gd + [0.01; 0.05; 0.8]; %[-0.9; 0.8; -1.6];
g = 0.5*gd;
Tf = Timed(end);
dt = 0.005;

[Time_cl, P_data_cl, dP_data_cl, ddP_data_cl] = dmp_classic.generate_trajectory(y0, g, Tf, dt);
[Time_rot, P_data_rot, dP_data_rot, ddP_data_rot] = dmp_rot.generate_trajectory(y0, g, Tf, dt);

[Time_bio, P_data_bio, dP_data_bio, ddP_data_bio] = dmp_bio.generate_trajectory(y0, g, Tf, dt);
[Time_bio_plus, P_data_bio_plus, dP_data_bio_plus, ddP_data_bio_plus] = dmp_bio_plus.generate_trajectory(y0, g, Tf, dt);

[Time, P_data, dP_data, ddP_data] = dmp_pp.generate_trajectory(y0, g, Tf, dt);


%% Accumulate the results
dat = {};

dat{1} = struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, ...
                'Color','blue', 'LineStyle','-', 'DisplayName','DMP$^{++}$');

dat{2} = struct('Time',Time_cl, 'Pos',P_data_cl, 'Vel',dP_data_cl, 'Accel',ddP_data_cl, ...
                'Color','magenta', 'LineStyle',':', 'DisplayName','DMP');

dat{3} = struct('Time',Time_rot, 'Pos',P_data_rot, 'Vel',dP_data_rot, 'Accel',ddP_data_rot, ...
                'Color','cyan', 'LineStyle','-.', 'DisplayName','DMP-rot');

dat{4} = struct('Time',Time_bio, 'Pos',P_data_bio, 'Vel',dP_data_bio, 'Accel',ddP_data_bio, ...
                'Color',[0.85 0.4 0.1], 'LineStyle','-', 'DisplayName','DMP-bio');
            
dat{5} = struct('Time',Time_bio_plus, 'Pos',P_data_bio_plus, 'Vel',dP_data_bio_plus, 'Accel',ddP_data_bio_plus, ...
                'Color',[0.2 0.8 0.2], 'LineStyle','-', 'DisplayName','DMP-bio+');

demo = struct('Time',Timed, 'Pos',Pd_data, 'Vel',dPd_data, 'Accel',ddPd_data, ...
                'Color',0.5*[1 1 1], 'LineStyle',':', 'DisplayName','demo');
dat{6} = demo;

%% Plot results
% Trajectories
for i=1:3
    figure;
    subplot(3,1,1); hold on;
    plot(Tf, g(i), 'LineWidth',4, 'Marker','x', 'MarkerSize',10, 'LineStyle','none', 'Color','red', 'HandleVisibility','off');
    for k=1:length(dat)
        if (isempty(dat{k})), continue; end
        plot(dat{k}.Time, dat{k}.Pos(i,:), 'LineWidth',2.0, 'LineStyle',dat{k}.LineStyle, 'Color',dat{k}.Color, 'DisplayName',dat{k}.DisplayName);
    end
    ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',15);
    title(['DoF $' num2str(i) '$'], 'interpreter','latex', 'fontsize',18);
    legend({}, 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

    subplot(3,1,2); hold on;
    for k=1:length(dat)
        if (isempty(dat{k})), continue; end
        plot(dat{k}.Time, dat{k}.Vel(i,:), 'LineWidth',2.0, 'LineStyle',dat{k}.LineStyle, 'Color',dat{k}.Color);
    end
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

    subplot(3,1,3); hold on;
    for k=1:length(dat)
        if (isempty(dat{k})), continue; end
        plot(dat{k}.Time, dat{k}.Accel(i,:), 'LineWidth',2.0, 'LineStyle',dat{k}.LineStyle, 'Color',dat{k}.Color);
    end
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;
end

% 3D path
fig = figure;
fig.Position(3:4) = [686 616];

ax = subplot(2, 1, 1); hold(ax, 'on');
ax.FontSize = 12;
for k=1:length(dat)
    if (isempty(dat{k})), continue; end
    plot3(dat{k}.Pos(1,:), dat{k}.Pos(2,:), dat{k}.Pos(3,:), 'LineWidth',2.0, 'LineStyle',dat{k}.LineStyle, 'Color',dat{k}.Color, 'DisplayName',dat{k}.DisplayName);
end
demo_pl = plot3(demo.Pos(1,:), demo.Pos(2,:), demo.Pos(3,:), 'LineWidth',2.0, 'LineStyle',demo.LineStyle, 'Color',demo.Color, 'HandleVisibility','off');
p0_pl = plot3(y0(1), y0(2), y0(3), 'LineWidth',4, 'Marker','o', 'MarkerSize',10, 'LineStyle','none', 'Color','green', 'HandleVisibility','off');
pg_pl = plot3(g(1), g(2), g(3), 'LineWidth',4, 'Marker','x', 'MarkerSize',10, 'LineStyle','none', 'Color','red', 'HandleVisibility','off');
pgd_pl = plot3(gd(1), gd(2), gd(3), 'LineWidth',4, 'Marker','x', 'MarkerSize',10, 'LineStyle','none', 'Color',[1 0.6 0.6], 'HandleVisibility','off');
legend({}, 'interpreter','latex', 'fontsize',17, 'Position',[0.1052 0.9344 0.7555 0.0506], 'Orientation','horizontal', 'Box','off');
xlabel('$X$ [$m$]', 'interpreter','latex', 'fontsize',17);
ylabel('$Y$ [$m$]', 'interpreter','latex', 'fontsize',17);
zlabel('$Z$ [$m$]', 'interpreter','latex', 'fontsize',17);
% view(-8.7, 19.8);
view(7.7, 11.6);
axis tight;
grid on;
hold off;

ax2 = subplot(2, 1, 2); hold on;
demo_pl = copyobj(demo_pl, ax2); set(demo_pl, 'HandleVisibility','on', 'DisplayName',demo.DisplayName);
p0_pl = copyobj(p0_pl, ax2); set(p0_pl, 'HandleVisibility','on', 'DisplayName','$p_0$');
pg_pl = copyobj(pg_pl, ax2); set(pg_pl, 'HandleVisibility','on', 'DisplayName','$g$');
pgd_pl = copyobj(pgd_pl, ax2); set(pgd_pl, 'HandleVisibility','on', 'DisplayName','$g_d$');
legend({}, 'interpreter','latex', 'fontsize',17, 'Position',[0.1096 0.8870 0.4507 0.0506], 'Orientation','horizontal', 'Box','off');

ax.Position = [0.1300 0.1100 0.7750 0.8150];
ax2.Position = [-5 -5 0.001 0.001];

%% ============================================================
%% ============================================================


function [Time, Y_data, dY_data, ddY_data] = simulateBioDMP(dmp_bio, y0, g, T, dt)

    %% set initial values
    n_dofs = length(dmp_bio); % number of DoFs

    y = y0; % position
    dy = zeros(n_dofs,1); % velocity
    ddy = zeros(n_dofs,1); % acceleration

    z = zeros(n_dofs,1); % z state: scaled velocity
    dz = zeros(n_dofs,1); % 
    
    t = 0.0;

    t_end = T;
    tau = t_end;
    
    can_clock_ptr = dmp_bio{1}.can_clock_ptr;
    can_clock_ptr.setTau(tau);

    % phase variable
    s = 0.0;
    s_dot = can_clock_ptr.getPhaseDot(s);

    iters = 0;

    % data to log
    Time = [];
    Y_data = [];
    dY_data = [];
    ddY_data = [];

    for i=1:n_dofs
        dmp_bio{i}.setY0(y0(i));
    end

    count = 0;

    %% simulate
    while (true)

        %% data logging
        Time = [Time t];
        Y_data = [Y_data y];
        dY_data = [dY_data dy];
        ddY_data = [ddY_data ddy];

        %% DMP simulation
        s_dot = can_clock_ptr.getPhaseDot(s);
        for i=1:n_dofs
            dmp_bio{i}.update(s, y(i), z(i), g(i), 0, 0);
            dy(i) = dmp_bio{i}.getYdot();
            dz(i) = dmp_bio{i}.getZdot();
            ddy(i) = dmp_bio{i}.getYddot();
        end

        %% Stopping criteria
        if (t>=1.0*t_end) % && norm(y-g)<1e-3 && norm(dy)<5e-3)
            break;
        end

        count = count + 1;

        %% Numerical integration
        iters = iters + 1;
        t = t + dt;
        s = s + s_dot*dt;
        y = y + dy*dt;
        z = z + dz*dt;

    end


end

