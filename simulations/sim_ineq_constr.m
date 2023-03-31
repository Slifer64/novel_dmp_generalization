clc;
close all;
clear;

addpath('utils/');
import_gmp_lib();

%% =============  Load train data  =============
load('data/demo_data3.mat', 'Timed', 'Pd_data', 'dPd_data', 'ddPd_data');

n_dofs = size(Pd_data, 1);
sd_data = Timed/Timed(end);

%% -------- load GMP model --------
gmp = GMP(n_dofs, 25, 1.5);
train_mse = gmp.train('LS', sd_data, Pd_data)

%% =============  DMP simulation  =============
disp('Simulation...');
t_start = tic;

%% Initial/Final values
y0d = Pd_data(:,1);   % Initial demo position
gd = Pd_data(:,end); % Target demo position
y0 = y0d; % set initial position for execution (for simplicity lets leave it the same as the demo)
g = gd + 0.3;  % set target position for execution
Tf = 1*Timed(end); % set the time duration of the executed motion

dt = 0.002; % time step for numerical integration


%% -------- via-points ---------
via_points = {struct('s',0.4, 'pos',0.5), struct('s',0.8, 'pos',0.5)}; %struct('s',1.6/Tf, 'pos',0.45)}; %, struct('s',3/Tf, 'pos',0.5)};

%% -------- Limits ---------
pos_lim = [-0.1 0.7];
vel_lim = repmat([ -0.5 , 0.5 ], n_dofs, 1);
accel_lim = repmat([ -0.8 , 0.8], n_dofs, 1);
slack_limits = [2e-2, 0.05, 0.2];

%% --------- Optimization objective ----------
opt_metric = struct('pos',0, 'vel',1);

%% Simulate
disp('Simulating DMP++ with ineq constraints...')
[Time, P_data, dP_data, ddP_data] = simulateModel(@DMP_pp, gmp, y0, g, Tf, dt, pos_lim, vel_lim, accel_lim, slack_limits, opt_metric, via_points);

disp('Simulating DMP with ineq constraints...')
[Time2, P2_data, dP2_data, ddP2_data] = simulateModel(@DMP_classic_wrapper, gmp, y0, g, Tf, dt, pos_lim, vel_lim, accel_lim, slack_limits, opt_metric, via_points);

%% Plot results

t_v = [];
p_v = [];
for k=1:length(via_points)
    t_v = [t_v via_points{k}.s*Tf];
    p_v = [p_v via_points{k}.pos];
end

ax_font = 13;
x_font = 16;
y_font = 16;
legend_font = 17;

Timed = Timed/Timed(end)*Tf;

for i=1:n_dofs
    fig = figure;
    ax_vec = [];
    fig.Position(3:4) = [647 832];
    % Plot Position
    ax = subplot(3,1,1); hold(ax, 'on');
    ax_vec = [ax_vec ax];
    pl = plot(Time, P_data(i,:), 'LineWidth',2.0 , 'Color','blue', 'DisplayName', 'DMP$^{++}$');
    pl2 = plot(Time2, P2_data(i,:), 'LineWidth',2.0, 'LineStyle',':', 'Color',[0.85 0.4 0.1], 'DisplayName','DMP');
    pl0 = plot(Timed, Pd_data(i,:), 'LineWidth',2.0, 'LineStyle','-.' , 'Color','green', 'DisplayName', 'demo');
    plot(Tf, gd(i), 'LineWidth',3, 'Marker','x', 'MarkerSize',12, 'LineStyle','none', 'Color',[1 0.6 0.6], 'DisplayName','$g_d$');
    plot(Tf, g(i), 'LineWidth',3, 'Marker','x', 'MarkerSize',12, 'LineStyle','none', 'Color','red', 'DisplayName','$g$');
    if (~isempty(t_v))
       plot(t_v, p_v(i,:), 'LineWidth',2, 'LineStyle','none', 'Marker','*', 'MarkerSize',14, 'Color',[1 0.2 0.2], 'DisplayName','via-points');
    end
    % plot limits
    plot_limits(ax, Time, pos_lim, slack_limits(1));
    ax.FontSize = ax_font;
    ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',y_font);
    legend({}, 'interpreter','latex', 'fontsize',legend_font, 'Position',[0.2251 0.9396 0.5708 0.0425], 'Orientation','horizontal', 'Box','off');
    axis tight;
    hold off;
    
    % Plot Velocity
    ax = subplot(3,1,2); hold(ax, 'on');
    ax_vec = [ax_vec ax];
    plot([Time(1) Time(end)], [0 0], 'LineWidth',1.0, 'LineStyle',pl2.LineStyle, 'Color',0.4*[1 1 1]);
    plot(Time, dP_data(i,:), 'LineWidth',pl.LineWidth, 'Color',pl.Color);
    plot(Time2, dP2_data(i,:), 'LineWidth',pl2.LineWidth, 'LineStyle',':', 'Color',pl2.Color);
%     plot(Timed, dPd_data(i,:), 'LineWidth',pl0.LineWidth, 'LineStyle','-.' , 'Color','green');
    % plot limits
    plot_limits(ax, Time, vel_lim, slack_limits(2));
    ax.FontSize = ax_font;
    ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',y_font);
    axis tight;
    hold off;
    
    % Plot Acceleration
    ax = subplot(3,1,3); hold(ax, 'on');
    ax_vec = [ax_vec ax];
    plot([Time(1) Time(end)], [0 0], 'LineWidth',1.0, 'LineStyle',':', 'Color',0.4*[1 1 1]);
    plot(Time, ddP_data(i,:), 'LineWidth',pl.LineWidth, 'Color',pl.Color);
    plot(Time2, ddP2_data(i,:), 'LineWidth',pl2.LineWidth, 'LineStyle',pl2.LineStyle, 'Color',pl2.Color);
%     plot(Timed, ddPd_data(i,:), 'LineWidth',pl0.LineWidth, 'LineStyle',pl0.LineStyle , 'Color',pl0.Color);
    % plot limits
    plot_limits(ax, Time, accel_lim, slack_limits(3));
    ax.FontSize = ax_font;
    ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',y_font);
    xlabel('time [$s$]', 'interpreter','latex', 'fontsize',x_font);
    axis tight;
    hold off;
    
    linkaxes(ax_vec, 'x');
    for i=1:length(ax_vec)
       ax = ax_vec(i);
       ax.Box = 'on';
       ax.YLim = ax.YLim + 0.07*(ax.YLim(2)-ax.YLim(1))*[-1 1];
    end
    ax = ax_vec(i);
    ax.XLim(2) = ax.XLim(2) + 0.02;
end



%% ===============================================
%% ===============================================

function plot_limits(ax, Time, soft_limit, slack_margin)
    
    plot([Time(1) Time(end)], soft_limit(1)*[1 1], 'LineWidth',1.2, 'LineStyle','--', 'Color','magenta', 'HandleVisibility','off', 'Parent',ax);
    plot([Time(1) Time(end)], soft_limit(2)*[1 1], 'LineWidth',1.2, 'LineStyle','--', 'Color','magenta', 'HandleVisibility','off', 'Parent',ax);
    plot([Time(1) Time(end)], (soft_limit(1)-slack_margin)*[1 1], 'LineWidth',1.2, 'LineStyle','--', 'Color',0.3*[1 1 1], 'HandleVisibility','off', 'Parent',ax);
    plot([Time(1) Time(end)], (soft_limit(2)+slack_margin)*[1 1], 'LineWidth',1.2, 'LineStyle','--', 'Color',0.3*[1 1 1], 'HandleVisibility','off', 'Parent',ax);

end

function [Time, Y_data, dY_data, ddY_data] = simulateModel(DMP_model, gmp, y0, yg, Tf, dt, pos_lim, vel_lim, accel_lim, slack_limits, opt_metric, via_points)

    %% set initial values
    n_dofs = length(y0);
    
    O_ndof = zeros(n_dofs,1);
    y = y0; % position
    y_dot = O_ndof; % velocity
    y_ddot = O_ndof; % acceleration

    t = 0.0;
    tau = Tf;

    % phase variable, from 0 to 1
    s = 0.0;

    % data to log
    Time = [];
    Y_data = [];
    dY_data = [];
    ddY_data = [];
    pos_slack_data = [];
    vel_slack_data = [];
    accel_slack_data = [];
    
    %% -------  DMP-Model  -------
    dmp_model = DMP_model(gmp);
    dmp_model.init(s, y0, yg, Tf);
    
    %% -------  Canonical System  -------
    can_sys = CanonicalSystem(tau, 30);
    
    %% -------  GMP-MPC  -------
    N_horizon = 10;
    pred_time_step = 0.1;
    N_kernels = 30;
    kernels_std_scaling = 1.5;
    
    final_state_err_tol = [1e-4; 1e-3; 1e-1];
    
    slack_gains = [1e5 100 1];
        
    %% --------  GMP - MPC  --------
    gmp_mpc = GMP_MPC(dmp_model, N_horizon, pred_time_step, N_kernels, kernels_std_scaling, slack_gains);
    
    gmp_mpc.settings.max_iter = 12000;
    gmp_mpc.settings.time_limit = 0; %2e-3;
    gmp_mpc.settings.abs_tol = 1e-3;
    gmp_mpc.settings.rel_tol = 1e-5;
    
    gmp_mpc.setObjCostGains(opt_metric.pos, opt_metric.vel);
    
    gmp_mpc.setPosLimits(pos_lim(:,1), pos_lim(:,2));
    gmp_mpc.setVelLimits(vel_lim(:,1), vel_lim(:,2));
    gmp_mpc.setAccelLimits(accel_lim(:,1), accel_lim(:,2));
    
    gmp_mpc.setPosSlackLimit(slack_limits(1));
    gmp_mpc.setVelSlackLimit(slack_limits(2));
    gmp_mpc.setAccelSlackLimit(slack_limits(3));

    gmp_mpc.setInitialState(y, y_dot, y_ddot, can_sys.s, can_sys.s_dot, 0);
    gmp_mpc.setFinalState(yg, O_ndof, O_ndof, 1, can_sys.s_dot, 0, final_state_err_tol);
    
    can_sys_fun = @(s, s_dot) [s_dot; can_sys.getPhaseDDot(s, s_dot)];
    gmp_mpc.setCanonicalSystemFunction(can_sys_fun);
    
    t_start = tic;
    
    opt_fail = false;
    on_opt_fail_exit = true;
    
    progress_cfg = struct('line_len',0, 'show_every',10, 'iter_count',9);
    
    %% simulate
    while (true)
        
        progress_cfg.iter_count = progress_cfg.iter_count + 1;
        if progress_cfg.iter_count == progress_cfg.show_every
            progress_cfg.iter_count = 0;
            fprintf(repmat('\b',1,progress_cfg.line_len))
            progress_cfg.line_len = fprintf('progress: %.2f %%\n', can_sys.s*100);
        end

        %% data logging
        Time = [Time t];
        Y_data = [Y_data y];
        dY_data = [dY_data y_dot];
        ddY_data = [ddY_data y_ddot];

        %% Update
        dmp_model.update(yg, can_sys.s, can_sys.s_dot, y, y_dot);
        gmp_mpc.setFinalState(yg, O_ndof, O_ndof, 1, can_sys.sd_dot, 0, final_state_err_tol);
        
        %% proccess via-points
        s_vp_ahead = 1 / tau; % process vp 1 sec before it is reached
        if (~isempty(via_points))
            
            rm_ind = [];
            for i=1:length(via_points)
                if (via_points{i}.s >= can_sys.s && via_points{i}.s <= can_sys.s + s_vp_ahead)
                   rm_ind = [rm_ind i];
                   dmp_model.updateViapoint(via_points{i}.s, via_points{i}.pos, false);
                   gmp_mpc.addViaPoint(via_points{i}.s, via_points{i}.pos, 1e-3);
                end
            end
            via_points(rm_ind) = [];
        end
        
        %% update final state
        
        
        %% optimization
        if (~opt_fail)

            % sovle
            sol = gmp_mpc.solve(can_sys.s, can_sys.s_dot);

            % check exit status
            if (sol.exit_flag)
                warning(sol.exit_msg);
                if (sol.exit_flag < 0)
                    opt_fail = true;
                    if (on_opt_fail_exit), break; end 
                end
            end
            
        end
        
        %% get optimal trajectory
        s = can_sys.s;
        s_dot = can_sys.s_dot;
        s_ddot = can_sys.getPhaseDDot(s, s_dot);
        
        y_ref = gmp_mpc.getYd(s); %sol.y;
        y_ref_dot = gmp_mpc.getYdDot(s, s_dot); %sol.y_dot;
        y_ref_ddot = gmp_mpc.getYdDDot(s, s_dot, s_ddot); %sol.y_ddot;
        
        y_ddot = 300*(y_ref - y) + 60*(y_ref_dot - y_dot) + y_ref_ddot;
        
        %% update initial state (optionally, if say a perturbation occurs)
        gmp_mpc.setInitialState(y, y_dot, y_ddot, s, s_dot, s_ddot);

        %% Log data
        Time = [Time t];
        Y_data = [Y_data y];
        dY_data = [dY_data y_dot];
        ddY_data = [ddY_data y_ddot];

        if (~opt_fail)
            pos_slack_data = [pos_slack_data sol.pos_slack];
            vel_slack_data = [vel_slack_data sol.vel_slack];
            accel_slack_data = [accel_slack_data sol.accel_slack];
        end
        
        %% Numerical integration
        can_sys.integrate(t, t+dt);
        t = t + dt;
        
        y = y + y_dot*dt;
        y_dot = y_dot + y_ddot*dt;

        %% Stopping criteria
        if (t>=1.0*Tf && norm(y-yg)<1e-3 && norm(y_dot)<5e-3)
            break;
        end
        
        if (t >= 1.4*Tf)
            warning('Time limit exceeded!');
            break; 
        end

    end

    fprintf('Error: pos=%e , vel=%e, accel=%e \n', norm(y - yg), norm(y_dot), norm(y_ddot));

end

