clc;
close all;
clear;

%% =============  includes...  =============
addpath('utils/')
import_gmp_lib();
import_dmp_lib();

%% =============  Load train data  =============
load('data/comp_scalings_data.mat', 'Timed', 'Pd_data', 'dPd_data', 'ddPd_data');

Ts = Timed(2) - Timed(1);

%% =============  Create/Train GMP  =============
n_dof = size(Pd_data, 1); % number of DoFs

%% initialize and train GMP
gmp = GMP(n_dof, 25, 1.5);
t_start = tic;
offline_train_mse = gmp.train('LS', Timed/Timed(end), Pd_data);
offline_train_mse
toc(t_start)

%% classical DMP scaling
classic_traj_sc = TrajScale_Prop(n_dof);

%% Rot DMP scaling
rot_traj_sc = TrajScale_Rot_min();
% rot_traj_sc = TrajScale_Rot_wb();
% rot_traj_sc.setWorkBenchNormal([0; 0; 1]); % set also the workbench normal

%% Bio-DMP
can_clock_ptr = CanonicalClock();
t_start = tic;
for i=1:n_dof
    dmp_bio{i} = DMP_bio(25, 60, 5, can_clock_ptr, ExpGatingFunction()); %ExpGatingFunction, LinGatingFunction
    train_error(i) = dmp_bio{i}.train('LS', Timed, Pd_data(i,:), dPd_data(i,:), ddPd_data(i,:));
end
train_error
toc(t_start)

%% DMP simulation
disp('Simulation...');
t_start = tic;

%% Initial/Final values
P0d = Pd_data(:,1);   % Initial demo position
Pgd = Pd_data(:,end); % Target demo position
P0 = P0d; % set initial position for execution (for simplicity lets leave it the same as the demo)
Pg = Pgd + [-0.9; 0.8; -1.6];

Tf = 0.5*Timed(end); % set the time duration of the executed motion

dt = 0.005; % time step for numerical integration

sc = TrajScale_Rot_wb();
sc.setWorkBenchNormal([0 0 1]);


%% Execute the DMP
get_target_fun = @(t) Pg;
[Time_cl, P_data_cl, dP_data_cl, ddP_data_cl] = simulateModel(DMP_classic(gmp, classic_traj_sc), dt, Tf, P0, 'get_target_fun',get_target_fun);
[Time_rot, P_data_rot, dP_data_rot, ddP_data_rot] = simulateModel(DMP_classic(gmp, rot_traj_sc), dt, Tf, P0, 'get_target_fun',get_target_fun);
[Time_bio, P_data_bio, dP_data_bio, ddP_data_bio] = simulateBioDMP(dmp_bio, P0, Pg, Tf, dt);
[Time, P_data, dP_data, ddP_data] = simulateModel(DMP_pp(gmp), dt, Tf, P0, 'get_target_fun',get_target_fun);

[Time_rd, P_data_rd, dP_data_rd, ddP_data_rd] = simulateModel(DMP_classic(gmp, TrajScale_roto_dial(n_dof)), dt, Tf, P0, 'get_target_fun',get_target_fun);

toc(t_start)

%% Accumulate the results
dat{1} = struct('Time',Time, 'Pos',P_data, 'Vel',dP_data, 'Accel',ddP_data, ...
                'Color','blue', 'LineStyle','-', 'DisplayName','DMP$^{++}$');

dat{2} = struct('Time',Time_cl, 'Pos',P_data_cl, 'Vel',dP_data_cl, 'Accel',ddP_data_cl, ...
                'Color','magenta', 'LineStyle',':', 'DisplayName','DMP');

dat{3} = struct('Time',Time_rot, 'Pos',P_data_rot, 'Vel',dP_data_rot, 'Accel',ddP_data_rot, ...
                'Color','cyan', 'LineStyle','-.', 'DisplayName','DMP-rot');

dat{4} = struct('Time',Time_bio, 'Pos',P_data_bio, 'Vel',dP_data_bio, 'Accel',ddP_data_bio, ...
                'Color',[0.85 0.4 0.1], 'LineStyle','-', 'DisplayName','DMP-bio');
            
dat{5} = struct('Time',Time_rd, 'Pos',P_data_rd, 'Vel',dP_data_rd, 'Accel',ddP_data_rd, ...
                'Color',[0.2 0.8 0.2], 'LineStyle','-', 'DisplayName','DMP-rot-dial');
            
demo = struct('Time',Timed, 'Pos',Pd_data, 'Vel',dPd_data, 'Accel',ddPd_data, ...
                'Color',0.5*[1 1 1], 'LineStyle',':', 'DisplayName','demo');


%% Plot results

% Trajectories
for i=1:3
    figure;
    subplot(3,1,1); hold on;
    plot(Time(end), Pg(i), 'LineWidth',4, 'Marker','x', 'MarkerSize',10, 'LineStyle','none', 'Color','red', 'HandleVisibility','off');
    for k=1:length(dat)
        plot(dat{k}.Time, dat{k}.Pos(i,:), 'LineWidth',2.0, 'LineStyle',dat{k}.LineStyle, 'Color',dat{k}.Color, 'DisplayName',dat{k}.DisplayName);
    end
    ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',15);
    title(['DoF $' num2str(i) '$'], 'interpreter','latex', 'fontsize',18);
    legend({}, 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

    subplot(3,1,2); hold on;
    for k=1:length(dat)
        plot(dat{k}.Time, dat{k}.Vel(i,:), 'LineWidth',2.0, 'LineStyle',dat{k}.LineStyle, 'Color',dat{k}.Color);
    end
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',15);
    axis tight;
    hold off;

    subplot(3,1,3); hold on;
    for k=1:length(dat)
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
    plot3(dat{k}.Pos(1,:), dat{k}.Pos(2,:), dat{k}.Pos(3,:), 'LineWidth',2.0, 'LineStyle',dat{k}.LineStyle, 'Color',dat{k}.Color, 'DisplayName',dat{k}.DisplayName);
end
demo_pl = plot3(demo.Pos(1,:), demo.Pos(2,:), demo.Pos(3,:), 'LineWidth',2.0, 'LineStyle',demo.LineStyle, 'Color',demo.Color, 'HandleVisibility','off');
p0_pl = plot3(P0(1), P0(2), P0(3), 'LineWidth',4, 'Marker','o', 'MarkerSize',10, 'LineStyle','none', 'Color','green', 'HandleVisibility','off');
pg_pl = plot3(Pg(1), Pg(2), Pg(3), 'LineWidth',4, 'Marker','x', 'MarkerSize',10, 'LineStyle','none', 'Color','red', 'HandleVisibility','off');
pgd_pl = plot3(Pgd(1), Pgd(2), Pgd(3), 'LineWidth',4, 'Marker','x', 'MarkerSize',10, 'LineStyle','none', 'Color',[1 0.6 0.6], 'HandleVisibility','off');
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

