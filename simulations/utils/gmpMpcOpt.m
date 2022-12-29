function [Time, P_data, dP_data, ddP_data] = gmpMpcOpt(gmp0, Tf, y0, yg0, yg, t_g, pos_lim, vel_lim, accel_lim, slack_limits, opt_pos, opt_vel, vp_config)
      
    if (isempty(vp_config))
       vp = {};
       t_vp_ahead = nan;
    else
        vp = vp_config.via_points;
        t_vp_ahead = vp_config.t_ahead;
    end
    
    gmp = gmp0.deepCopy();
   
    n_dof = length(y0);

    O_ndof = zeros(n_dof,1);

    t = 0;
    dt = 0.002;
    tau = Tf;
    can_sys = CanonicalSystem(tau, 30);
    y = y0;
    y_dot = O_ndof;
    y_ddot = O_ndof;
    
    N_horizon = 10;
    pred_time_step = 0.1;
    N_kernels = 30;
    kernels_std_scaling = 1.5;
    
    final_state_err_tol = [1e-4; 1e-3; 1e-1];
    
    slack_gains = [1e5 100 1];

    time_limit = 0; %2e-3;
    max_iter = 12000;
    abs_tol = 1e-3;
    rel_tol = 1e-5;
        
    %% --------  GMP - MPC  --------
    gmp_mpc = GMP_MPC(gmp, N_horizon, pred_time_step, N_kernels, kernels_std_scaling, slack_gains);
    
    gmp_mpc.settings.max_iter = max_iter;
    gmp_mpc.settings.time_limit = time_limit;
    gmp_mpc.settings.abs_tol = abs_tol;
    gmp_mpc.settings.rel_tol = rel_tol;
    
    gmp_mpc.setObjCostGains(opt_pos, opt_vel);
    
    gmp_mpc.setPosLimits(pos_lim(:,1), pos_lim(:,2));
    gmp_mpc.setVelLimits(vel_lim(:,1), vel_lim(:,2));
    gmp_mpc.setAccelLimits(accel_lim(:,1), accel_lim(:,2));
    
    gmp_mpc.setPosSlackLimit(slack_limits(1));
    gmp_mpc.setVelSlackLimit(slack_limits(2));
    gmp_mpc.setAccelSlackLimit(slack_limits(3));

    gmp_mpc.setInitialState(y, y_dot, y_ddot, can_sys.s, can_sys.s_dot, 0);
    gmp_mpc.setFinalState(yg0, O_ndof, O_ndof, 1, can_sys.s_dot, 0, final_state_err_tol);
    
    can_sys_fun = @(s, s_dot) [s_dot; can_sys.getPhaseDDot(s, s_dot)];
    gmp_mpc.setCanonicalSystemFunction(can_sys_fun);
    
    gmp.setScaleMethod(TrajScale_Prop(n_dof));
    gmp.setY0(y0);
    gmp.setGoal(yg0);
    
    text_prog = ProgressText(40);
    text_prog.init();
    
    t_start = tic;
    
    s_vp_ahead = t_vp_ahead / tau;
    
    opt_fail = false;
    on_opt_fail_exit = true;
    
    Time = [];
    P_data = [];
    dP_data = [];
    ddP_data = [];
    pos_slack_data = [];
    vel_slack_data = [];
    accel_slack_data = [];
    
    %% -------  Simulation loop  --------
    while (true)
        
%         if (can_sys.s > 0.45), break; end
        
        %% proccess via-points
        if (~isempty(vp))
            
            rm_ind = [];
            for i=1:length(vp)
                if (vp{i}.s >= can_sys.s && vp{i}.s <= can_sys.s + s_vp_ahead)
                   rm_ind = [rm_ind i];
                   gmp_mpc.addViaPoint(vp{i}.s, vp{i}.pos, vp{i}.err_tol);
                end
            end
            vp(rm_ind) = [];
        end
        
        %% update final state
        if (t >= t_g)
            t_g = inf; % to avoid calling it again
            gmp.setGoal(yg);
        end
        gmp_mpc.setFinalState(yg, O_ndof, O_ndof, 1, can_sys.sd_dot, 0, final_state_err_tol);
            
        %% Stopping criteria
%         if (s >= 1 && norm(y_dot)<1e-2 && norm(y-yg)<1e-2 )
%             break;
%         elseif (s >= 1.3)
%             warning('Time limits exceeded...');
%             break;
%         end
        if (can_sys.s >= 1), break; end
        
        % display progress
        if (can_sys.s <= 1), text_prog.update(100*can_sys.s); end
        
        %% optimization
        if (~opt_fail)

            % sovle
            sol = gmp_mpc.solve(can_sys.s, can_sys.s_dot);

            % check exit status
            if (sol.exit_flag)
                warning(sol.exit_msg);
                text_prog.printInNewLine();
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
        
%         y_ddot = 300*(y_ref - y) + 80*(y_ref_dot - y_dot) + y_ref_ddot;
        y = y_ref;
        y_dot = y_ref_dot;
        y_ddot = y_ref_ddot;
        
        %% update initial state (optionally, if say a perturbation occurs)
        gmp_mpc.setInitialState(y, y_dot, y_ddot, s, s_dot, s_ddot);

        %% Log data
        Time = [Time t];
        P_data = [P_data y];
        dP_data = [dP_data y_dot];
        ddP_data = [ddP_data y_ddot];

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
    
    end

    if (can_sys.s >= 1), text_prog.update(100); end

    if (~isempty(Time))
        
        plotSlackVariables(Time, pos_slack_data, vel_slack_data, accel_slack_data, slack_limits);
        
        target_err = norm(P_data(:,end)-yg)
        vel_err = norm(dP_data(:,end))
        accel_err = norm(ddP_data(:,end))
        
        max_slack_violations = [max(abs(pos_slack_data(:))), max(abs(vel_slack_data(:))), max(abs(accel_slack_data(:)))]
    end
    
    fprintf('\n');
    fprintf('===> GMP-MPC optimization finished! Elaps time: %f ms\n',toc(t_start)*1000);
    
%     fprintf('H sp: %.2f (%d/%d)\n', gmp_mpc.H_sp, gmp_mpc.H_size(1), gmp_mpc.H_size(2));
%     fprintf('A sp: %.2f (%d/%d)\n', gmp_mpc.A_sp, gmp_mpc.A_size(1), gmp_mpc.A_size(2));


end


