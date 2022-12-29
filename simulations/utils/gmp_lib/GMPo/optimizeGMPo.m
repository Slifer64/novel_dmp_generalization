%% Optimize GMPo weights according to desired trajectory and constraints
%  @param[in] gmp_o: pointer to GMPo object.
%  @param[in] Q0: Initial orientation.
%  @param[in] Qg: target orientation.
%  @param[in] T: motion time duration.
%  @param[in] theta_bar: upper bound of angle of rotation from initial orientation. Use '[]' to ignore.
%  @param[in] qdot_bar: upper bound on the first time derivative of log(Q*inv(Q0)). Use '[]' to ignore.
%  @param[in] qddot_bar: upper bound on the second time derivative of log(Q*inv(Q0)). Use '[]' to ignore.
%  @param[in] lambda_pos: weights the position error in the optimization cost. Set to 0 to ignore position optimization.
%  @param[in] lambda_vel: weights the velocity error in the optimization cost. Set to 0 to ignore velocity optimization.
%  @param[out] exit_flag: 0 on success, 1 on premature stop, -1 on failure.
%  @param[out] exit_msg: message describing the exit result of the optimization.
function [exit_flag, exit_msg] = optimizeGMPo(gmp_o, Q0, Qg, T, x_data, theta_bar, qdot_bar, qddot_bar, lambda_pos, lambda_vel)
    
    pos_constr = 0;
    vel_constr = 0;
    accel_constr = 0;
    
    if (~isempty(theta_bar)), pos_constr = 1; end
    if (~isempty(qdot_bar)), vel_constr = 1; end
    if (~isempty(qddot_bar)), accel_constr = 1; end

    opt_pos = lambda_pos ~= 0;
    opt_vel = lambda_vel ~= 0;
    
    opt.w_p = lambda_pos;
    opt.w_v = lambda_vel;

    n_ker = gmp_o.numOfKernels();
    n_var = 3*n_ker;
    x_dot = 1/T;
    x_ddot = 0;

    % calculate cost function: J = 0.5w'Hw + f'w
    N = length(x_data);

    Q = 1e-8*eye(n_var, n_var); % for numerical stability
    f = zeros(n_var,1);
    c = 0;
    
    gmp_o.setQ0(Q0);   % set initial orientation
    gmp_o.setQg(Qg);   % set target orientation
    sc = norm( gmp_o.getScaling() );
    inv_T_sc = gmp_o.getInvScaling();

    if (opt_pos)
        H1 = zeros(n_var, n_var);
        f1 = zeros(n_var, 1);
        for i=1:N
            phi = gmp_o.regressVec(x_data(i));
            Phi = blkdiag(phi',phi',phi');
            H1 = H1 + Phi'*Phi;
            zd = inv_T_sc*gmp_o.getYd(x_data(i));
            f1 = f1 - Phi'*zd;
        end
        Q = Q + opt.w_p*H1;
        f = f + opt.w_p*f1;
    end

    if (opt_vel)
        H2 = zeros(n_var, n_var);
        f2 = zeros(n_var, 1);
        for i=1:N
            phi1 = gmp_o.regressVecDot(x_data(i), x_dot);
            Phi1 = blkdiag(phi1',phi1',phi1');
            H2 = H2 + Phi1'*Phi1;
            zd_dot = inv_T_sc*gmp_o.getYdDot(x_data(i), x_dot);
            f2 = f2 - Phi1'*zd_dot;
        end
        Q = Q + opt.w_v*H2;
        f = f + opt.w_v*f2;
    end

    % --------- inequality constraints ----------
    x_constr = 0:0.02:1;    
    M = length(x_constr);
    n_const = (pos_constr + vel_constr + accel_constr)*M;
    H_c = cell(n_const,1);
    k_c = cell(n_const,1);
    d_c = cell(n_const,1);
    
    i = 1;
    for j=1:M
        
        x = x_constr(j);
        
        if (pos_constr)
            phi = gmp_o.regressVec(x)';
            Phi = blkdiag(phi,phi,phi);
            H_c{i} = 2*(Phi'*Phi);
            k_c{i} = zeros(n_var,1);
            d_c{i} = -(theta_bar/sc)^2;
            i = i + 1;
        end
        
        if (vel_constr)
            phi = gmp_o.regressVecDot(x,x_dot)';
            Phi = blkdiag(phi,phi,phi);
            H_c{i} = 2*(Phi'*Phi);
            k_c{i} = zeros(n_var,1);
            d_c{i} = -(qdot_bar/sc)^2;
            i = i + 1;
        end
        
        if (accel_constr)
            phi = gmp_o.regressVecDDot(x,x_dot, x_ddot)';
            Phi = blkdiag(phi,phi,phi);
            H_c{i} = 2*(Phi'*Phi);
            k_c{i} = zeros(n_var,1);
            d_c{i} = -(qddot_bar/sc)^2;
            i = i + 1;
        end
    end
    
    % --------- equality constraints ----------
    q0 = zeros(3,1);
    qg = inv_T_sc*gmp_.quatLogDiff(Qg, Q0);
    beq = [q0; qg];
    
    phi_0 = gmp_o.regressVec(0)';
    phi_1 = gmp_o.regressVec(1)';
    
    Aeq = [blkdiag(phi_0,phi_0,phi_0); blkdiag(phi_1,phi_1,phi_1)];

    % tic
    
    options = optimoptions(@fmincon,'Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
        'HessianFcn',@(x,lambda)quadhess(x,lambda,Q,H_c), ...
        'Display','off');

    fun = @(x)quadobj(x,Q,f,c);
    nonlconstr = @(x)quadconstr(x,H_c,k_c,d_c);
    x0 = gmp_o.W';
    x0 = x0(:);
    [x,fval,eflag,output,lambda] = fmincon(fun,x0,...
        [],[],Aeq,beq,[],[],nonlconstr,options);

    % elaps_t = toc;
    % fprintf('===> Optimization finished! Elaps time: %f ms\n',elaps_t*1000);
    
    exit_flag = -1;
    
    if (eflag >= 1)
        exit_msg = 'Found solution! (within default tolerances)';
        exit_flag = 0;
    elseif (eflag == 0)
        exit_msg = 'Maximum iterations reached.';
        exit_flag = 1;
    elseif (eflag == -2), exit_msg = 'The problem appears to be Infeasible.';
    elseif (eflag == -3), exit_msg = 'The problem appears to be Unbounded.';
    end

    % fprintf('===> Optimization finished! Elaps time: %f ms\n',elaps_t*1000);

    % x = Q \ (-f); % unconstrained solution

    % on success, assign solution to GMP
    if (~exit_flag), gmp_o.W = reshape(x, n_ker, 3)'; end

end

%% ===========  Utility functions  ============

function [y,grady] = quadobj(x,Q,f,c)
    y = 1/2*x'*Q*x + f'*x + c;
    if nargout > 1
        grady = Q*x + f;
    end
end

function [y,yeq,grady,gradyeq] = quadconstr(x,H,k,d)
    jj = length(H); % jj is the number of inequality constraints
    y = zeros(1,jj);
    for i = 1:jj
        y(i) = 1/2*x'*H{i}*x + k{i}'*x + d{i};
    end
    yeq = [];

    if nargout > 2
        grady = zeros(length(x),jj);
        for i = 1:jj
            grady(:,i) = H{i}*x + k{i};
        end
    end
    gradyeq = [];
end

function hess = quadhess(x,lambda,Q,H)
    hess = Q;
    jj = length(H); % jj is the number of inequality constraints
    for i = 1:jj
        hess = hess + lambda.ineqnonlin(i)*H{i};
    end
end