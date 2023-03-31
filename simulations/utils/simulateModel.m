function [Time, Y_data, dY_data, ddY_data, Yg_data] = simulateModel(model, dt, Tf, y0, varargin)
    
    args = parse_args(Tf, varargin{:});

    %% set initial values
    n_dofs = length(y0);

    y = y0; % position
    y_dot = zeros(n_dofs,1); % velocity
    y_ddot = zeros(n_dofs,1); % acceleration

    t = 0.0;
    t_end = Tf;
    tau = Tf;

    % canonical system
    can_sys = CanonicalSystem(Tf);

    % data to log
    Time = [];
    s_data = [];
    Y_data = [];
    dY_data = [];
    ddY_data = [];
    Yg_data = [];
    Yg_filt_data = [];
    
    yg = args.get_target_fun(t);

    model.init(can_sys.s, y0, yg, Tf);
    model.setAdaptToRobot(false);
    % model.setRecursiveUpdate(true);
    
    model.K = 300; % DMP stiffness
    model.D = 2*sqrt(model.K + 10); % DMP damping
    
    %% simulate
    while (true)

        %% data logging
        Time = [Time t];
        Y_data = [Y_data y];
        dY_data = [dY_data y_dot];
        ddY_data = [ddY_data y_ddot];
        s_data = [s_data can_sys.s];

        yg_new = args.get_target_fun(t);
        if (isinf(args.goal_filt_coeff))
            yg = yg_new;
            yg_dot = zeros(n_dofs,1);
        else
            yg_dot = args.goal_filt_coeff*(yg_new - yg);
        end
        
        Yg_data = [Yg_data yg_new];
        Yg_filt_data = [Yg_filt_data yg];

        %% Update model
        model.update(yg, can_sys.s, can_sys.s_dot, y, y_dot);
        
        %% Get reference
        s_ddot = can_sys.getPhaseDDot();
        
        external_signal = 0; % optionally add some external signal

        % DMP transformation system: 
        y_ddot = model.goal_attractor(y, y_dot, tau) + model.shape_attractor(can_sys.s, can_sys.s_dot, s_ddot, tau) + external_signal;

        %% Stopping criteria
        if (t>=(1.0*t_end+0.05) && norm(y-yg_new)<1e-3 && norm(y_dot)<5e-3)
            break;
        end
        
        if (t >= 1.4*t_end)
            warning('Time limit exceeded!');
            break; 
        end

        %% Numerical integration
        can_sys.integrate(t, t+dt);
        t = t + dt;
        y = y + y_dot*dt;
        y_dot = y_dot + y_ddot*dt;
        yg = yg + yg_dot*dt;

    end
    
%     s = can_sys.s;
%     s_dot = can_sys.s_dot;
%     s_ddot = can_sys.getPhaseDDot();
%     
%     fprintf('s = %.4f, s_dot = %.5f, s_ddot = %.5f\n', s, s_dot, s_ddot);
% 
%     y_ref = model.getRefPos(s)
%     dy_ref = model.getRefVel(s, s_dot)
%     ddy_ref = model.getRefAccel(s, s_dot, s_ddot)
%     
%     y
%     y_dot
%     y_ddot
%     
%     yg
    
    fprintf('Error: pos=%e , vel=%e, accel=%e \n', norm(y - yg_new), norm(y_dot), norm(y_ddot));
    
%     figure;
%     plot(Time, Yg_filt_data);

%     figure; hold on;
%     plot(Time, s_data);
%     plot([Time(1) Time(end)], [1 1], 'LineStyle',':', 'color',0.4*[1 1 1]);
%     axis tight;

end

function args = parse_args(T, varargin)

    parser = inputParser;
    parser.KeepUnmatched = false;
    parser.PartialMatching = false;
    parser.CaseSensitive = false;

    parser.addParameter('get_target_fun', []);
%     parser.addParameter('get_time_duration_fun', @(t) T);
%     parser.addParameter('get_viapoints_fun', @(t) []);
%     parser.addParameter('ellipsoids', []);
%     parser.addParameter('via_points', []);
%     parser.addParameter('online_plot', false);
    parser.addParameter('goal_filt_coeff', inf);

    parser.parse(varargin{:});
    args = parser.Results;
    
    if (isempty(args.get_target_fun))
        error('get_target_fun must be provided as input!');
    end

end
