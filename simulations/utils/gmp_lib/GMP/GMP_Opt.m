%% N-DoF GMP Optimization class
%

classdef GMP_Opt < matlab.mixin.Copyable
    
    properties (Access=public, Constant)
        
        % enum QPSolverT{
        MATLAB_QUADPROG = 0;
        OSQP = 1;
        GOLDFARB_IDNANI = 2;
        %};
        
    end
    
    methods (Access = public)
        
        %% GMP constructor.
        %  @param[in] gmp: n_DoF dmp.
        function this = GMP_Opt(gmp)
                
%             this.ex_flag_map = containers.Map('KeyType','double','ValueType','char');
%             this.ex_flag_map(-10) = 'Empty... Call ''constrOpt'' first.';
%             this.ex_flag_map(1) = 'Function converged to the solution x.';
%             this.ex_flag_map(0) = 'Number of iterations exceeded options.MaxIterations.';
% 
%             this.ex_flag_map(-2) = 'Problem is infeasible. Or, for interior-point-convex, the step size was smaller than options.StepTolerance, but constraints were not satisfied.';
%             this.ex_flag_map(-3) = 'Problem is unbounded';
%             this.ex_flag_map(2) = 'Step size was smaller than options.StepTolerance, constraints were satisfied.';
%             this.ex_flag_map(-6) = 'Nonconvex problem detected.';
%             this.ex_flag_map(-8) = 'Unable to compute a step direction.';
            
            this.exit_msg = '';
                
            this.gmp = gmp;
            
            this.tau = 8;
            
            this.setOptions(1,0,0);
            
        end
        
        function setOptions(this, pos_obj_w, vel_obj_w, accel_obj_w)
            
            this.w_p = pos_obj_w;
            this.w_v = vel_obj_w;
            this.w_a = accel_obj_w;
            
            this.setQPsolver(GMP_Opt.MATLAB_QUADPROG);

        end       
        
        %% Sets the QP-solver to use for solving the optimization problem.
        %  @param[in] solver_type: enum of type @GMP_Opt::QPSolverT.
        function setQPsolver(this, solver_type)
            
            if (solver_type == GMP_Opt.MATLAB_QUADPROG), this.qp_solver_fun = @GMP_Opt.quadprogSolver;
            elseif (solver_type == GMP_Opt.OSQP), this.qp_solver_fun = @GMP_Opt.osqpSolver;
            elseif (solver_type == GMP_Opt.GOLDFARB_IDNANI), this.qp_solver_fun = @GMP_Opt.goldfarbIdnaniSolver;
            else, error('Unsupported qp solver type');
            end
            
        end
        
        %% Trajectory optimization.
        %  @param[in] num_points: Number of discreet points to use in the objective function.
        %  @param[in] tau: Duration of motion.
        %  @param[in] pos_constr: Vector of @GMPConstr position constraints. For no constraints pass '[]'.
        %  @param[in] vel_constr: Vector of @GMPConstr velocity constraints. For no constraints pass '[]'.
        %  @param[in] accel_constr: Vector of @GMPConstr acceleration constraints. For no constraints pass '[]'.
        %  @param[out] success: true if minimum found, false otherwise.
        
        function success = optimize(this, num_points)
            
            if (nargin < 2), num_points = 150; end
            
            dx = 1 / num_points;
            x_data = 0:dx:1;
            success = this.optimize2(x_data);

        end
            
        function success = optimize2(this, x_data)

            opt_pos = this.w_p;
            opt_vel = this.w_v;
            opt_accel = this.w_a;
            
            n_ker = this.gmp.numOfKernels();
            n_dof = this.gmp.numOfDoFs();
            
            x_dot = 1/this.tau;
            x_ddot = 0;
            
            success = true;
            this.exit_msg = '';
            
            % calculate cost function: J = 0.5w'Hw + f'w
            N = length(x_data);
            
            H = 1e-8*eye(n_ker, n_ker); % for numerical stability
            f = zeros(n_ker, n_dof);
            
            if (opt_pos)
                y_offset = - this.gmp.getY0() + this.gmp.getScaling()*this.gmp.getY0d();
                H1 = zeros(n_ker, n_ker);
                f1 = zeros(n_ker, n_dof);
                for i=1:N
                    phi = this.gmp.regressVec(x_data(i));
                    H1 = H1 + phi*phi';
                    yd = this.gmp.getYd(x_data(i)) + y_offset;
                    f1 = f1 - phi*yd';
                end
                H = H + this.w_p*H1;
                f = f + this.w_p*f1;
            end
  
            if (opt_vel)
                H2 = zeros(n_ker, n_ker);
                f2 = zeros(n_ker, n_dof);
                for i=1:N
                    phi = this.gmp.regressVecDot(x_data(i), x_dot);
                    H2 = H2 + phi*phi';
                    f2 = f2 - phi*this.gmp.getYdDot(x_data(i), x_dot)';
                end
                H = H + this.w_v*H2;
                f = f + this.w_v*f2;
            end
            
            if (opt_accel)
                H3 = zeros(n_ker, n_ker);
                f3 = zeros(n_ker, n_dof);
                for i=1:N
                    phi = this.gmp.regressVecDDot(x_data(i), x_dot, x_ddot);
                    H3 = H3 + phi*phi';
                    f3 = f3 - phi*this.gmp.getYdDDot(x_data(i), x_dot, x_ddot)';
                end
                H = H + this.w_a*H3;
                f = f + this.w_a*f3;
            end
  
            % inequality constraints
            Aineq = [this.A_p; this.A_v; this.A_a];
            lb_ineq = [this.pos_lb; this.vel_lb; this.accel_lb];
            ub_ineq = [this.pos_ub; this.vel_ub; this.accel_ub];
         
            % equality constraints
            Aeq = [this.Aeq_p; this.Aeq_v; this.Aeq_a];
            beq = [this.pos_eq; this.vel_eq; this.accel_eq];
 

            % solve optimization problem
            % tic
            
            solve_coupled = 0;
            
%             disp([ 'H-sp: ' num2str(sum(H(:)== 0)) ' / ' num2str(numel(H)) ] );
%             disp([ 'A-sp: ' num2str(sum(Aineq(:)== 0)) ' / ' num2str(numel(Aineq)) ] );
            
            %% solve coupled QP (all dofs together)
            if (solve_coupled)

                W0 = this.gmp.W';
                
                [W, success, this.exit_msg] = this.qp_solver_fun(kron(eye(n_dof),H),f(:),...
                                                                        kron(eye(n_dof),Aineq),lb_ineq(:),ub_ineq(:), ...
                                                                        kron(eye(n_dof),Aeq),beq(:), W0(:));
                W = reshape(W, n_ker,n_dof)';
            
            else
            %% solve QP for each DoF separately
                
                W = zeros(n_dof, n_ker);
                W0 = this.gmp.W';
                solve_fun_pt = @this.qp_solver_fun;
                exit_msg_array = cell(n_dof,1);
                success_array = zeros(n_dof,1);
                % parfor i=1:n_dof
                for i=1:n_dof

                    beq_i = [];
                    lb_ineq_i = [];
                    ub_ineq_i = [];
                    if (~isempty(beq)), beq_i = beq(:,i); end
                    if (~isempty(lb_ineq)), lb_ineq_i = lb_ineq(:,i); end
                    if (~isempty(ub_ineq)), ub_ineq_i = ub_ineq(:,i); end

                    % solve optimization problem
                    [W(i,:), success_i, exit_status_msg] = solve_fun_pt(H,f(:,i), Aineq,lb_ineq_i,ub_ineq_i, Aeq,beq_i, W0(:,i));
                    
                    success_array(i) = success_i;

                    exit_msg_array{i} = ['DoF-' num2str(i) ': ' exit_status_msg '\n'];

                end
                this.exit_msg = [exit_msg_array{:}];
                this.exit_msg = this.exit_msg(1:end-2); % remove last '\n'
                success = isempty( find(success_array== 0) );
            end
            
            if (success)
                this.gmp.W = this.gmp.getInvScaling()*W;
            end
            
            % toc
            
        end
        
        function setMotionDuration(this, tau)
           
            this.tau = tau;
            
        end
        
        function setPosBounds(this, lower_bound, upper_bound, num_points)
            
            if (nargin < 4), num_points = 50; end
            
            n_dofs = this.gmp.numOfDoFs();
            
            if (isscalar(lower_bound)), lower_bound = lower_bound * ones(n_dofs,1); end
            if (isscalar(upper_bound)), upper_bound = upper_bound * ones(n_dofs,1); end
            
            x = linspace(0,1, num_points); 
            n = length(x);
            
            lb = repmat(lower_bound, 1, n);
            ub = repmat(upper_bound, 1, n);
            
            this.setPosConstr(x, lb, ub, [], []);
            
        end
        
        function setVelBounds(this, lower_bound, upper_bound, num_points)
            
            if (nargin < 4), num_points = 50; end
            
            n_dofs = this.gmp.numOfDoFs();
            
            if (isscalar(lower_bound)), lower_bound = lower_bound * ones(n_dofs,1); end
            if (isscalar(upper_bound)), upper_bound = upper_bound * ones(n_dofs,1); end
            
            x = linspace(0,1, num_points); 
            n = length(x);
            
            lb = repmat(lower_bound, 1, n);
            ub = repmat(upper_bound, 1, n);
            
            this.setVelConstr(x, lb, ub, [], []);
            
        end
        
        function setAccelBounds(this, lower_bound, upper_bound, num_points)
            
            if (nargin < 4), num_points = 50; end
            
            n_dofs = this.gmp.numOfDoFs();
            
            if (isscalar(lower_bound)), lower_bound = lower_bound * ones(n_dofs,1); end
            if (isscalar(upper_bound)), upper_bound = upper_bound * ones(n_dofs,1); end
            
            x = linspace(0,1, num_points); 
            n = length(x);
            
            lb = repmat(lower_bound, 1, n);
            ub = repmat(upper_bound, 1, n);
            
            this.setAccelConstr(x, lb, ub, [], []);
            
        end
        
        function setPosConstr(this, x, lb, ub, x_eq, p_eq)

            n_ker = this.gmp.numOfKernels();
            n_dof = this.gmp.numOfDoFs();
            
            y_offset = - this.gmp.getY0() + this.gmp.getScaling()*this.gmp.getY0d();
            
            % Process inequality constraints
            if (~isempty(x))
 
                m = length(x);
                
                if (size(lb,2) ~= m), error(['[GMP_Opt::setPosConstr]: Lower bounds must have ' num2str(m) ' columns (constraints).']); end
                if (size(ub,2) ~= m), error(['[GMP_Opt::setPosConstr]: Upper bounds must have ' num2str(m) ' columns (constraints).']); end
                if (size(lb,1) ~= size(ub,1)), error('[GMP_Opt::setPosConstr]: Lower and Upper bounds must have the same number of rows (DoFs)'); end
                
                this.A_p = zeros(m, n_ker);
                this.pos_lb = zeros(m, n_dof);
                this.pos_ub = zeros(m, n_dof);

                for i=1:m
                    this.pos_lb(i,:) = ( lb(:,i) + y_offset )';
                    this.pos_ub(i,:) = ( ub(:,i) + y_offset )';
                    this.A_p(i,:) = this.gmp.regressVec(x(i))';
                end
            
            end
            
            % Process equality constraints
            if (~isempty(x_eq))
                
                m = length(x_eq);
                if (size(p_eq,2) ~= m), error(['[GMP_Opt::setPosConstr]: Equality values matrix must have ' num2str(m) ' columns (constraints).']); end
                
                this.Aeq_p = zeros(m, n_ker);
                this.pos_eq = zeros(m, n_dof);

                for i=1:m
                    this.pos_eq(i,:) = ( p_eq(:,i) + y_offset )';
                    this.Aeq_p(i,:) = this.gmp.regressVec(x_eq(i))';
                end
            
            end
            
        end
        
        function setVelConstr(this, x, lb, ub, x_eq, v_eq)
            
            n_ker = this.gmp.numOfKernels();
            x_dot = 1 / this.tau;

            % Process inequality constraints
            if (~isempty(x))
 
                m = length(x);
                
                if (size(lb,2) ~= m), error(['[GMP_Opt::setVelConstr]: Lower bounds must have ' num2str(m) ' columns (constraints).']); end
                if (size(ub,2) ~= m), error(['[GMP_Opt::setVelConstr]: Upper bounds must have ' num2str(m) ' columns (constraints).']); end
                if (size(lb,1) ~= size(ub,1)), error('[GMP_Opt::setVelConstr]: Lower and Upper bounds must have the same number of rows (DoFs)'); end
                
                this.vel_lb = lb';
                this.vel_ub = ub';
                this.A_v = zeros(m, n_ker);
                for i=1:m, this.A_v(i,:) = this.gmp.regressVecDot(x(i), x_dot)'; end
            
            end
            
            % Process equality constraints
            if (~isempty(x_eq))
                
                m = length(x_eq);
                if (size(v_eq,2) ~= m), error(['[GMP_Opt::setVelConstr]: Equality values matrix must have ' num2str(m) ' columns (constraints).']); end

                this.vel_eq = v_eq';
                this.Aeq_v = zeros(m, n_ker);
                for i=1:m, this.Aeq_v(i,:) = this.gmp.regressVecDot(x_eq(i),x_dot)'; end
            
            end
            
        end
        
        function setAccelConstr(this, x, lb, ub, x_eq, a_eq)
            
            n_ker = this.gmp.numOfKernels();
            x_dot = 1 / this.tau;
            x_ddot = 0;

            % Process inequality constraints
            if (~isempty(x))
 
                m = length(x);
                
                if (size(lb,2) ~= m), error(['[GMP_Opt::setAccelConstr]: Lower bounds must have ' num2str(m) ' columns (constraints).']); end
                if (size(ub,2) ~= m), error(['[GMP_Opt::setAccelConstr]: Upper bounds must have ' num2str(m) ' columns (constraints).']); end
                if (size(lb,1) ~= size(ub,1)), error('[GMP_Opt::setAccelConstr]: Lower and Upper bounds must have the same number of rows (DoFs)'); end
                
                this.accel_lb = lb';
                this.accel_ub = ub';
                this.A_a = zeros(m, n_ker);
                for i=1:m, this.A_a(i,:) = this.gmp.regressVecDDot(x(i), x_dot, x_ddot)'; end
            
            end
            
            % Process equality constraints
            if (~isempty(x_eq))
                
                m = length(x_eq);
                if (size(a_eq,2) ~= m), error(['[GMP_Opt::setAccelConstr]: Equality values matrix must have ' num2str(m) ' columns (constraints).']); end

                this.accel_eq = a_eq';
                this.Aeq_a = zeros(m, n_ker);
                for i=1:m, this.Aeq_a(i,:) = this.gmp.regressVecDDot(x_eq(i), x_dot, x_ddot)'; end
            
            end
            
        end
        
        function clearPosConstr(this)
                
            this.A_p = [];
            this.pos_lb = [];
            this.pos_ub = [];
            
            this.Aeq_p = [];
            this.pos_eq = [];

        end
        
        function clearVelConstr(this)
                
            this.A_v = [];
            this.vel_lb = [];
            this.vel_ub = [];
            
            this.Aeq_v = [];
            this.vel_eq = [];
            
        end
        
        function clearAccelConstr(this)
            
            this.A_a = [];
            this.accel_lb = [];
            this.accel_ub = [];
            
            this.Aeq_a = [];
            this.accel_eq = [];
            
        end
        
        function clearConstr(this)
    
            this.clearPosConstr();
            this.clearVelConstr();
            this.clearAccelConstr();
            
        end
        
        function msg = getExitMsg(this)
           
            msg = this.exit_msg;
            
        end

    end
    
    
    %% ============= QP solvers ==============
    methods (Access = private, Static)
        
        %% Matlab-quatprog solver
        function [x, success, exit_status] = quadprogSolver(H, f, Aineq, lb_ineq, ub_ineq, Aeq, beq, x0)
            
            A = [Aineq;   -Aineq];
            b = [ub_ineq; -lb_ineq];
            
            opt = optimoptions('quadprog', 'Algorithm','interior-point-convex', 'StepTolerance',0, 'Display','off');
            
            %[x, ~, ex_flag, opt_output] = quadprog(H,f, A,b, Aeq,beq, [],[], x0, opt);
            [x, ~, ex_flag, opt_output] = quadprog(sparse(H),f, sparse(A),b, sparse(Aeq),beq, [],[], x0, opt);
            exit_status = GMP_Opt.exit_flag_map(ex_flag); %opt_output.message;
            success = ex_flag == 1 || ex_flag == 2; 
                    
        end
        
        %% OSQP solver
        function [x, success, exit_status] = osqpSolver(H, f, Aineq, lb_ineq, ub_ineq, Aeq, beq, x0)
            
            osqp_solver = osqp;
            osqp_solver.setup(sparse(H), f, sparse([Aineq; Aeq]), [lb_ineq; beq], [ub_ineq; beq], 'warm_start',true, 'verbose',false, 'eps_abs',1e-4, 'eps_rel',1e-5);%, 'max_iter',20000);
            res = osqp_solver.solve();
            exit_status = res.info.status;
            ex_flag = res.info.status_val;
            success = ex_flag > 0 || ex_flag==-2; % -2:max_iters
            x = res.x;
                    
        end
        
        %% Goldfarb-Idnani solver
        function [x, success, exit_status] = goldfarbIdnaniSolver(H, f, Aineq, lb_ineq, ub_ineq, Aeq, beq, x0)

            n = length(f);
            [x, ~, ex_flag, opt_output] = qpGoldfarbIdnani(H,f, [Aineq; -Aineq],[ub_ineq; -lb_ineq], Aeq,beq, -inf(n,1),inf(n,1));
            exit_status = opt_output.status{1};
            success = ex_flag == 1;
                
        end
        
    end
    
    %% =================================================
    %% ============== Private properties ===============
    %% =================================================
    
    properties (Access = private, Constant)
        
        % maps the exit flag value to the msg describing the exit flag
        exit_flag_map = containers.Map({-10,1,0,-2,-3,2,-6,-8}, ...
            {'Empty... Call ''constrOpt'' first.',...
             'Function converged to the solution x.', ...
             'Number of iterations exceeded options.MaxIterations.', ...
             'Problem is infeasible. Or, for interior-point-convex, the step size was smaller than options.StepTolerance, but constraints were not satisfied.', ...
             'Problem is unbounded', ...
             'Step size was smaller than options.StepTolerance, constraints were satisfied.', ...
             'Nonconvex problem detected.', ...
             'Unable to compute a step direction.'});
    end
    
    properties (Access = private)
        
        qp_solver_fun % pointer to the qp solver function
        
        exit_msg
        
        gmp % n_DoF GMP
   
        w_p % position objective weight
        w_v % velocity objective weight
        w_a % acceleration objective weight
        
        tau % motion duration
        
        % Constraints
        A_p
        pos_lb
        pos_ub
        
        Aeq_p
        pos_eq
        
        A_v
        vel_lb
        vel_ub
        
        Aeq_v
        vel_eq
        
        A_a
        accel_lb
        accel_ub
        
        Aeq_a
        accel_eq

    end
    
end
