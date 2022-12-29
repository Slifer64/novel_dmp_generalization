%% Quadratic Optimization Motion Primitive class
%  y_ddot(t) = Phi(x) * w
%  y_dot(t) = y_dot(t0) + int_{t0}_{t}(Phi(x) * w)
%  y(t) = y(t0) + y_dot(t0)(t-t0) + int_{t0}_{t}(int_{t0}_{t}(Phi(x) * w))
%
%  y_ddot = az*bz*(yd - y) + az*(yd_dot - y_dot) + yd_ddot
%

classdef QPMP < matlab.mixin.Copyable

    methods (Access = public)
        %% Weighted Sum of Gaussians constructor.
        %  @param[in] N_kernels: The number of kernels.
        %  @param[in] yd_ddot: Training data acceleration.
        %  @param[in] taud: Training data time duration.
        %  @param[in] yd0: Training data initial position.
        %  @param[in] gd: Training data goal position.
        %  @param[in] kernel_std_scaling: Scaling of the kernel's std. (optional, default=1.0)
        function this = QPMP(N_kernels, a_z, b_z, kernel_std_scaling)

            if (nargin < 2), kernel_std_scaling = 1.0; end
            
            this.a_z = a_z;
            this.b_z = b_z;
            this.N_kernels = N_kernels;
            
            this.zero_tol = 1e-30; %realmin;

            this.w = zeros(this.N_kernels,1);
            this.c = ((1:this.N_kernels)-1)'/(this.N_kernels-1);
            % this.c = linspace(-0.04,1.04, N_kernels)';
            
            this.h = 1./(kernel_std_scaling*(this.c(2:end)-this.c(1:end-1))).^2;
            this.h = [this.h; this.h(end)];
            
            d = this.c(2) - this.c(1);
            c_end = this.c(end);
            extra = [d; 2*d; 3*d];
            this.c = [-extra+this.c(1); this.c; extra+c_end];
            this.h = [repmat(this.h(end),length(extra),1); this.h; repmat(this.h(end),length(extra),1)];
            
            this.N_kernels = this.N_kernels + 2*length(extra);
            this.w = zeros(this.N_kernels,1);
            
            this.setIntStep(0.002);
            
        end

        
        %% Returns the number of kernels.
        %  @return The number of kernels.
        function n_ker = numOfKernels(this)
            
            n_ker = length(this.w);
            
        end
       
        
        %% Sets the integration step that is used to numerically calculate 
        %% position/velocity contraints during training.
        %  @param[in] dt: integration step
        function setIntStep(this, dt)
            
            this.dt = dt;
            
        end
        
        
        %% Sets the demo data.
        function setDemo(this, demo_accel, taud, yd0, gd)
            
            this.demo_accel = demo_accel;
            this.taud = taud;
            this.yd0 = yd0;
            this.gd = gd;
            
        end
        
        
        %% Trains the model.
        %  One can optionally pass position/velocity constraints.
        %  @param[in] tau: Duration of motion.
        %  @param[in] y0: Initial position.
        %  @param[in] ydot_0: Initial velocity.
        %  @param[in] g: Goal/Final position.
        %  @param[in] pos_constr: Vector of position constraints. Each element is a 
        %                         struct('t',t_c,'value',p_c, 'type','=/</>'). For no constraints pass '[]'.
        %  @param[in] vel_constr: Vector of velocity constraints. Each element is a 
        %                         struct('t',t_c,'value',v_c, 'type','=/</>'). For no constraints pass '[]'.
        %  @param[in] accel_constr: Vector of acceleration constraints. Each element is a 
        %                         struct('t',t_c,'value',a_c, 'type','=/</>'). For no constraints pass '[]'.
        function train(this, tau, y0, ydot_0, g, pos_constr, vel_constr, accel_constr)
 
            this.tau = tau;
            
            this.y0 = y0;
            this.ydot_0 = ydot_0;
            
            % calculate spatial/temporal scalings
            spat_s = (g - y0) / (this.gd - this.yd0);
            temp_s = this.taud / tau;
            
            % calculate scaled trajectory
            z = spat_s * temp_s^2 * this.demo_accel;

            % calculate cost function: J = 0.5w'Hw + f'w
            N = length(z);
            x = (0:(N-1)) / (N-1);
            H = zeros(this.N_kernels, this.N_kernels);
            f = zeros(this.N_kernels, 1);
            for i=1:N
                phi = this.regressVec(x(i));
                H = H + phi*phi';
                f = f - phi*z(i);
            end

            % calculate constraint matrices in canonical form
            [A,b, Aeq,beq] = this.getConstrMat(pos_constr, vel_constr, accel_constr);
            
            % solve optimization problem
            tic
            this.w = quadprog(H,f, A,b, Aeq,beq);
            toc

        end
        

        %% Initialization for trajectory generation.
        function init(this)
           
            this.t = 0;
            this.yd = this.y0;
            this.yd_dot = this.ydot_0;

        end
        
        %% =============================================================
        
        %% Returns the acceleration and the current time instant.
        function [t, y_ddot] = getAccel(this, y, y_dot, dt)

            this.t = this.t + dt;
            x = this.t / this.tau;

            yd_ddot = dot(regressVec(this, x), this.w);
            
            y_ddot = this.a_z*this.b_z*(this.yd - y) + this.a_z*(this.yd_dot - y_dot) + yd_ddot;
            t = this.t;
            
            this.yd = this.yd + this.yd_dot*dt;
            this.yd_dot = this.yd_dot + yd_ddot*dt;
            
        end
        
        function [Time, y, y_dot, y_ddot] = trajGen(this, y0, y_dot0, dt)
            
            N = length(this.demo_accel);
            
            this.init(y0, y_dot0);
            
            y = zeros(1,N);
            y_dot = zeros(1,N);
            y_ddot = zeros(1,N);
            
            for i=1:N-1
                y_ddot(i) = this.getAccel(y(i), y_dot(i), dt);
                y(i+1) = y(i) + y_dot(i)*dt;
                y_dot(i+1) = y_dot(i) + y_ddot(i)*dt;
            end
            
            Time = dt*(0:(N-1))/(N-1);
            
        end
        

    end
    
    methods (Access = private)

        %% Extracts the constraint matrices in canonical form.
        %  That is: A*w <= b, Aeq*w = beq
        %  Each input argument constraint is of the form struct('t',t_c,'value',val, 'type','=/</>').
        %  @param[in] pos_constr: Vector of position constraints.
        %  @param[in] vel_constr: Vector of velocity constraints.
        %  @param[in] accel_constr: Vector of acceleration constraints.
        %  @param[out] A: Inequality constraint matrix.
        %  @param[out] b: Inequality constraint vector.
        %  @param[out] Aeq: Equality constraint matrix.
        %  @param[out] beq: Equality constraint vector.
        function [A,b, Aeq,beq] = getConstrMat(this, pos_constr, vel_constr, accel_constr)
           
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            
            t0 = 0;

            % acceleration contraints
            if (~isempty(accel_constr))
               
                for i=1:length(accel_constr)
                    t_c = accel_constr(i).t;
                    a_c = accel_constr(i).value;
                    type = accel_constr(i).type;
                    
                    phi = this.regressVec(t_c/this.tau);
                                        
                    b_i = a_c;
                    if (type == '=')
                        Aeq = [Aeq; phi'];
                        beq = [beq; b_i];
                    elseif (type == '<')
                        A = [A; phi'];
                        b = [b; b_i];
                    elseif (type == '>')
                        A = [A; -phi'];
                        b = [b; -b_i];
                    end

                end  
                
            end
            
            % velocity contraints
            vel_constr = this.sortConstr(vel_constr);
            phi2 = zeros(this.N_kernels, 1);
            t = t0;
            if (~isempty(vel_constr))
               
                for i=1:length(vel_constr)
                    t_c = vel_constr(i).t;
                    v_c = vel_constr(i).value;
                    type = vel_constr(i).type;

                    n_steps = round((t_c-t)/this.dt);
                    for j=1:n_steps
                       phi2 = phi2 + this.regressVec(t/this.tau)*this.dt;
                       t = t + this.dt;
                    end
                    
                    b_i = v_c-this.ydot_0;
                    if (type == '=')
                        Aeq = [Aeq; phi2'];
                        beq = [beq; b_i];
                    elseif (type == '<')
                        A = [A; phi2'];
                        b = [b; b_i];
                    elseif (type == '>')
                        A = [A; -phi2'];
                        b = [b; -b_i];
                    end

                end  
                
            end
            
            % position contraints
            pos_constr = this.sortConstr(pos_constr);
            phi2 = zeros(this.N_kernels, 1);
            phi3 = zeros(this.N_kernels, 1);
            t = t0;
            if (~isempty(pos_constr))
               
                for i=1:length(pos_constr)
                    t_c = pos_constr(i).t;
                    y_c = pos_constr(i).value;
                    type = pos_constr(i).type;

                    n_steps = round((t_c-t)/this.dt);
                    for j=1:n_steps
                       phi2 = phi2 + this.regressVec(t/this.tau)*this.dt;
                       phi3 = phi3 + phi2*this.dt;
                       t = t + this.dt;
                    end
                    
                    b_i = y_c - this.y0 - (this.tau-t0)*this.ydot_0;
                    if (type == '=')
                        Aeq = [Aeq; phi3'];
                        beq = [beq; b_i];
                    elseif (type == '<')
                        A = [A; phi3'];
                        b = [b; b_i];
                    elseif (type == '>')
                        A = [A; -phi3'];
                        b = [b; -b_i];
                    end

                end  
                
            end
            
        end
        
        
        function phi = regressVec(this, x)
            
            psi = this.kernelFun(x);
            phi = psi / (sum(psi) + this.zero_tol);

        end
        
        
        %% Returns a column vector with the values of the kernel functions.
        %  @param[in] x: The phase variable.
        %  @return: Column vector with the values of the kernel functions.
        function psi = kernelFun(this, x)

            n = length(x);
            psi = zeros(this.N_kernels, n);
            
            for j=1:n
                psi(:,j) = exp(-this.h.*((x(j)-this.c).^2));
            end 

        end
        
        function constr = sortConstr(this, constr)
          
            n = length(constr);
            t = zeros(n, 1);
            for i=1:n, t(i) = constr(i).t; end
            
            [~, ind] = sort(t, 'ascend');
            
            constr = constr(ind);
            
        end
        
    end
    
    properties (Access = private)
        
        N_kernels % number of kernels (basis functions)
        a_z
        b_z
        w % N_kernels x 1 vector with the kernels' weights
        c % N_kernels x 1 vector with the kernels' centers
        h % N_kernels x 1 vector with the kernels' inverse width
        
        dt % euler integration steps for calculating position/velocity constraints
        
        zero_tol % small value used to avoid divisions with very small numbers
        
        demo_accel; % training data acceleration
        taud % training data time duration
        yd0
        gd
        
        yd
        yd_dot
        y0
        ydot_0
        
        tau
        t
        
    end
    
end
