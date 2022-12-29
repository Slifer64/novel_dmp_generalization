function [Time, Y_data, dY_data, ddY_data, Yg_data] = simulateDMP(gmp, y0, get_target_fun, T, dt)

    %% set initial values
    Dim = gmp.numOfDoFs(); % number of DoFs

    y = y0; % position
    y_dot = zeros(Dim,1); % velocity
    y_ddot = zeros(Dim,1); % acceleration
    O_ndof = zeros(Dim, 1);

    t = 0.0;

    % dz = zeros(Dim,1);
    % z = zeros(Dim,1);

    t_end = T;
    tau = t_end;

    % phase variable, from 0 to 1
    s = 0.0;
    s_dot = 1/tau;
    s_ddot = 0; % since x_dot is constant here

    iters = 0;

    % data to log
    Time = [];
    Y_data = [];
    dY_data = [];
    ddY_data = [];
    x_data = [];
    Yg_data = [];
    
    yg = get_target_fun(t);

    gmp.setY0(y0);
    gmp.setGoal(yg);

    count = 0;
    
    ag = 20; % goal filter coeff, settling time ~= 4/ag sec

    %% simulate
    while (true)

        %% data logging
        Time = [Time t];
        Y_data = [Y_data y];
        dY_data = [dY_data y_dot];
        ddY_data = [ddY_data y_ddot];

        Yg_data = [Yg_data yg];

        %% DMP simulation
        yg_new = get_target_fun(t);
        yg_dot = ag*(yg_new - yg);
        gmp.setGoal(yg);
        
        y_s = gmp.getYd(s);
        dy_s = gmp.getYdDot(s, s_dot);
        ddy_s = gmp.getYdDDot(s, s_dot, s_ddot);
        
        if (s > 1)
            y_s = gmp.getYd(1.0);
            dy_s = 0*dy_s;
            ddy_s = 0*ddy_s;
        end

        K = 300; % set the DMP stiffness
        D = 60; % set the DMP damping

        external_signal = 0; % optionally add some external signal

        % Track it using a 2nd order dynamical system. This is actually the DMP. 
        y_ddot = ddy_s + D*(dy_s - y_dot) + K*(y_s - y) + external_signal;

        % if external_signal == 0 you can obviously set directly:
    %     y = y_s;
    %     dy = dy_s;
    %     ddy = ddy_s;

        %% Stopping criteria
        if (t>=1.0*t_end && norm(y-yg)<1e-3 && norm(y_dot)<5e-3)
            break;
        end
        
        if (t >= 1.4*t_end)
            warning('Time limit exceeded!');
            break; 
        end

        count = count + 1;

        %% Numerical integration
        iters = iters + 1;
        t = t + dt;
        s = s + s_dot*dt;
        s_dot = s_dot + s_ddot*dt;
        y = y + y_dot*dt;
        y_dot = y_dot + y_ddot*dt;
        yg = yg + yg_dot*dt;

    end
    
    fprintf('Error: pos=%e , vel=%e, accel=%e \n', norm(y - yg), norm(y_dot), norm(y_ddot));

end