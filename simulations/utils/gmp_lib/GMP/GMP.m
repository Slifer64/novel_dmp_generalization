%% N-DoF GMP class
%  Generalized movement primitive.
%  Encodes an N-DoF trajectory using the novel DMP formulation from 
%  the paper 'A Reversible Dynamic Movement Primitive formulation', DOI:10.1109/ICRA48506.2021.9562059
%
%  

classdef GMP < GMP_regressor
    
    methods (Access = public)
        
        %% GMP constructor.
        %  @param[in] n_dofs: number of degrees of freedom.
        %  @param[in] N_kernels: the number of kernels
        %  @param[in] kern_std_scale: Scaling for std of kernels (optional, default=1).
        %  @param[in] x_min: Minimum value of the phase variable. (optional, default=0.0)
        %  @param[in] x_max: Maximum value of the phase variable. (optional, default=1.0)
        function this = GMP(n_dofs, N_kernels, kern_std_scale, x_min, x_max)
                
            if (nargin < 1), n_dofs = 1; end
            if (nargin < 2), N_kernels = 2; end
            if (nargin < 3), kern_std_scale = 1.5; end
            if (nargin < 4), x_min = 0.0; end
            if (nargin < 5), x_max = 1.0; end
            
            this@GMP_regressor(N_kernels, kern_std_scale, x_min, x_max);

            this.K = 100 * ones(n_dofs,1);
            this.D = 30 * ones(n_dofs,1);
            
            this.W = zeros(n_dofs, N_kernels);
            this.W0 = this.W;
            
            this.Y0d = zeros(n_dofs,1);
            this.Ygd = ones(n_dofs,1);
            this.Y0 = this.Y0d;
            this.Yg = this.Ygd;
            
            this.setScaleMethod(TrajScale_Prop(n_dofs));

            this.y_dot = zeros(n_dofs,1);
            this.z_dot = zeros(n_dofs,1);

        end

        %% Returns the number of DoFs.
        %  return: number of DoFs.
        function n = numOfDoFs(this)
           
            n = size(this.W,1);
            
        end
        
        %% Returns the number of kernels.
        %  @return: number of kernels.
        function n_ker = numOfKernels(this)
            
            n_ker = length(this.c);
            
        end

        %% ==============================
        %% =========  Training  =========
        %% ==============================
        
        %% Trains the GMP.
        %  @param[in] train_method: the training method to use, as a string ('LWR', 'LS').
        %  @param[in] x: Row vector with the canonical timestamps (in [0 1]) of the training data points.
        %  @param[in] yd_data: Matrix with the desired potition for each DoF in each row.
        %  @param[out] train_error: The training error expressed as the mse error.
        %  \note The timestamps-data need not be sequencial temporarily.
        function [train_mse, Sw] = train(this, train_method, x, yd_data)
            
            if (~isempty(find(x>1 | x<0)))
               warning('[GMP::train]: The training timestamps are not normalized... Normalizing them in [0 1].');
               x = (x - min(x)) / (max(x) - min(x));
            end
            
            n_data = length(x);
            num_ker = this.numOfKernels();
            n_dofs = this.numOfDoFs();
            
            H = zeros(num_ker, n_data);
            % for j=1:n_data, Psi(:,j) = this.kernelFun(x(j)); end
            for j=1:n_data, H(:,j) = this.regressVec(x(j)); end

            if (strcmpi(train_method,'LWR') == 1)
                this.W = yd_data*H' ./ repmat(sum(H,2)',n_dofs,1);
            elseif (strcmpi(train_method,'LS') == 1)
                this.W = yd_data / H;
            else
                error('[WSoG::train]: Unsupported training method...');
            end
            
            this.W0 = this.W;
            
            this.Y0d = this.W*H(:,1);
            this.Ygd = this.W*H(:,end);
            
            this.setY0(this.Y0d);
            this.setGoal(this.Ygd);
            
            this.traj_sc.setNominalStartFinalPos(this.Y0d, this.Ygd);
            this.traj_sc.setNewStartFinalPos(this.getY0(), this.getGoal());

            if (nargout > 0)
                train_mse = mean((this.W*H - yd_data).^2, 2);
            end
            
            if (nargout > 1), Sw = inv(H*H'); end
            
        end

        %% Auto-retrains the GMP modifying the kernels according to the arguments.
        %  @param[in] N_kernels: the new number of kernels.
        %  @param[in] kern_std_scale: the scaling of the width of the kernels.
        %  @param[in] n_points: number of uniformly autogenerated points in
        %                      [0 1] for the auto-retraining. (optional, default = 200)
        %  @param[in] train_method: the training method to use, as a string
        %                           ('LWR', 'LS'). optinal (default = 'LS')
        function [train_err, Sw] = autoRetrain(this, N_kernels, kern_std_scale, n_points, train_method)
            
            if (nargin < 4), n_points = 500; end
            if (nargin < 5), train_method = 'LS'; end
            
            x = linspace(0,1, n_points);
            yd_data = zeros(this.numOfDoFs(), length(x));
            for j=1:length(x), yd_data(:,j) = this.W*this.regressVec(x(j)); end
            
            this.setKernels(N_kernels, kern_std_scale);
            [train_err, Sw] = this.train(train_method, x, yd_data);
        
        end

        
        %% =================================================
        %% =========  Set/Get initial/final state  =========
        %% =================================================
        
        %% Sets the initial position.
        %  @param[in] y0: initial position.
        function setY0(this, y0)
            
            this.traj_sc.setNewStartFinalPos(y0, this.Yg);
            this.Y0 = this.traj_sc.getY0();
            
        end
        
        %% Returns the initial position.
        %  @return the initial position.
        function Y0 = getY0(this)
            
            Y0 = this.Y0;
            
        end
        
        %% Returns the initial position learned from the training.
        %  @return the trained initial position.
        function y0d = getY0d(this)
            
            y0d = this.Y0d;
            
        end
        
        
        %% Set goal/final position.
        %  @param[in] g: goal position.
        function setGoal(this, g)

            this.traj_sc.setNewStartFinalPos(this.Y0, g);
            this.Yg = g;
        
        end

        %% Returns the goal/final position.
        %  @return goal position.
        function Yg = getGoal(this), Yg = this.Yg; end
        
        
        %% ===========================================
        %% =========  GMP output trajectory  =========
        %% ===========================================
        
        %% Returns the scaled desired position.
        %  @param[in] x: phase variable.
        %  @return: scaled desired position.
        function yd = getYd(this, x)

            yd = this.getScaling()*(this.W*this.regressVec(x) - this.Y0d) + this.Y0;
            
        end
        
        
        %% Returns the scaled desired velocity.
        %  @param[in] x: phase variable.
        %  @param[in] x_dot: 1st time derivative of the phase variable.
        %  @return: scaled desired velocity.
        function yd_dot = getYdDot(this, x, x_dot)
            
            yd_dot = this.getScaling()*this.W*this.regressVecDot(x,x_dot);

        end
        
        
        %% Returns the scaled desired acceleration.
        %  @param[in] x: phase variable.
        %  @param[in] x_dot: 1st time derivative of the phase variable.
        %  @param[in] x_ddot: 2nd time derivative of the phase variable.
        %  @return: scaled desired acceleration.
        function yd_ddot = getYdDDot(this, x, x_dot, x_ddot)
            
            if (nargin < 4), x_ddot = 0; end
            yd_ddot = this.getScaling()*this.W*this.regressVecDDot(x,x_dot,x_ddot);
            
        end

        %% ===========================================
        %% ========  Trajectory scaling ==============
        %% ===========================================
        
        function sc = getScaling(this)

            sc = this.traj_sc.getScaling();
            
        end
        
        function inv_sc = getInvScaling(this)
            
            inv_sc = this.traj_sc.getInvScaling();
            
        end

        %% Set the trajectory spatial scaling method.
        %  @param[in] traj_scale_obj: pointer to an object of type @TrajScale.
        function setScaleMethod(this, traj_scale_obj)
           
            if (this.numOfDoFs() ~= traj_scale_obj.numOfDoFs())
                error('[GMP::setScaleMethod]: Incompatible number of DoFs...');
            end
            
            this.traj_sc = traj_scale_obj;
            this.traj_sc.setNominalStartFinalPos(this.Y0d, this.Ygd);
            this.traj_sc.setNewStartFinalPos(this.Y0, this.Yg);

        end
        
        %% =================================
        %% ============  Misc ==============
        %% =================================
        
        function plotTrajectory(this, Tf, s_data)
            
            if (nargin < 3), s_data = linspace(0, 1, 100); end
            
            s_dot = 1/Tf;
            s_ddot = 0;
            
            n_data = length(s_data);
            n_dof = this.numOfDoFs();
            
            Time = s_data * Tf;
            P_data = zeros(n_dof, n_data);
            dP_data = zeros(n_dof, n_data);
            ddP_data = zeros(n_dof, n_data);
            for j=1:n_data
                s = s_data(j);
                P_data(:, j) = this.getYd(s);
                dP_data(:, j) = this.getYdDot(s, s_dot);
                ddP_data(:, j) = this.getYdDDot(s, s_dot, s_ddot);
            end
            
           if (n_dof == 3)
                plot_ = @(varargin) plot3(varargin{1}(1,:), varargin{1}(2,:), varargin{1}(3,:), varargin{2:end});
           elseif (n_dof == 2)
                plot_ = @(varargin) plot(varargin{1}(1,:), varargin{1}(2,:), varargin{2:end});
           else
               plot_ = @(varargin) 0; % do nothing
           end
             
           if (n_dof == 3 || n_dof == 2)
               figure;
               hold on;
               plot_(P_data, 'LineWidth',2, 'Color','blue', 'LineStyle','-');
               plot_(P_data(:,1), 'Marker','o', 'MarkerSize',10, 'LineWidth',4, 'Color','green', 'LineStyle','None');
               plot_(P_data(:,end), 'Marker','x', 'MarkerSize',10, 'LineWidth',4, 'Color','red', 'LineStyle','None');
           end

           figure;
           k = [1, 1+n_dof, 1+2*n_dof];
           for i=1:n_dof
               subplot(3, n_dof, k(1));
               plot(Time, P_data(i,:), 'LineWidth',2, 'Color','blue');
               if (i==1), ylabel('pos [$m$]', 'interpreter','latex', 'fontsize',14); end
               title(['DoF ' num2str(i)], 'interpreter','latex', 'fontsize',16);
               
               subplot(3, n_dof, k(2));
               plot(Time, dP_data(i,:), 'LineWidth',2, 'Color',[0 0.8 0]);
               if (i==1), ylabel('vel [$m/s$]', 'interpreter','latex', 'fontsize',14); end
               
               subplot(3, n_dof, k(3));
               plot(Time, ddP_data(i,:), 'LineWidth',2, 'Color','magenta');
               if (i==1), ylabel('accel [$m/s^2$]', 'interpreter','latex', 'fontsize',14); end
               xlabel('time [$s$]', 'interpreter','latex', 'fontsize',14);
               
               k = k + 1;
           end
            
        end
        
        %% Returns a deep copy of this object.
        %  @return: deep copy of this object.
        function cp_obj = deepCopy(this)
            
            % Make a shallow copy of all properties
            cp_obj = this.copy();
            
            % Make a deep copy of the pointers
            cp_obj.traj_sc = this.traj_sc.deepCopy();

        end
        
        function resetWeights(this)
           
            this.W = this.W0;
            
        end
        
        %% ===============================================
        %% ==========  original DMP functions ============
        %  Deprecated. Better use @getYd, @getYdDot, @getYdDDot instead.
        
        %% Calculates the time derivatives of the GMP's states.
        %  @param[in] s: Object of type @GMP_phase.
        %  @param[in] y: 'y' state of the GMP.
        %  @param[in] z: 'z' state of the GMP.
        %  @param[in] y_c: coupling term for the dynamical equation of the 'y' state (optional, default=0).
        %  @param[in] z_c: coupling term for the dynamical equation of the 'z' state (optional, default=0).
        function update(this, s, y, z, y_c, z_c)

            n_dofs = this.numOfDoFs();
            
            if (nargin < 5), y_c = zeros(n_dofs,1); end
            if (nargin < 6), z_c = zeros(n_dofs,1); end
            
            if (length(y_c) == 1), y_c = ones(n_dofs,1)*y_c(1); end
            if (length(z_c) == 1), z_c = ones(n_dofs,1)*z_c(1); end
            
            yd = this.getYd(s.x);
            yd_dot = this.getYdDot(s.x, s.x_dot);
            yd_ddot = this.getYdDDot(s.x, s.x_dot, s.x_ddot);
  
            this.y_dot = z + y_c;
            this.z_dot = this.K.*(yd - y) + this.D.*(yd_dot - z) + yd_ddot + z_c;

        end

        
        %% Returns the 'y' state time derivative.
        %  Call after @update.
        %  @return: time derivative of 'y' state.
        function y_dot = getYdot(this)
            
            y_dot = this.y_dot;
            
        end
        
        
        %% Returns the 'z' state time derivative.
        %  Call after @update.
        %  @return: time derivative of 'z' state.
        function z_dot = getZdot(this)
            
            z_dot = this.z_dot;
            
        end
        
        
        %% Returns the GMP's acceleration.
        %  Call after @update.
        %  @param[in] yc_dot: time derivative of 'y' state coupling (optional, default=0).
        %  @return: acceleration.
        function y_ddot = getYddot(this, yc_dot)
            
            n_dofs = this.numOfDoFs();
            if (nargin < 2), yc_dot = zeros(n_dofs,1); end
            if (length(yc_dot)==1), yc_dot = ones(n_dofs,1)*yc_dot(1); end
            
            y_ddot = this.getZdot() + yc_dot;

        end
        
        
        %% Calculates the GMP's acceleration.
        %  @param[in] s: Object of type @GMP_phase.
        %  @param[in] y: 'y' state of the GMP.
        %  @param[in] y_dot: time derivative of 'y' state.
        %  @param[in] y_c: coupling term for the dynamical equation of the 'y' state (optional, default=0).
        %  @param[in] z_c: coupling term for the dynamical equation of the 'z' state (optional, default=0).
        %  @param[in] yc_dot: time derivative of the 'y' state coupling (optional, default=0).
        %  @return: acceleration.
        function y_ddot = calcYddot(this, s, y, y_dot, y_c, z_c, yc_dot)

            n_dofs = this.numOfDoFs();
            if (nargin < 5), y_c = zeros(n_dofs,1); end
            if (nargin < 6), z_c = zeros(n_dofs,1); end
            if (nargin < 7), yc_dot = zeros(n_dofs,1); end
            
            if (length(y_c)==1), y_c = ones(n_dofs,1)*y_c(1); end
            if (length(z_c)==1), z_c = ones(n_dofs,1)*z_c(1); end
            if (length(yc_dot)==1), yc_dot = ones(n_dofs,1)*yc_dot(1); end

            yd = this.getYd(s.x);
            yd_dot = this.getYdDot(s.x, s.x_dot);
            yd_ddot = this.getYdDDot(s.x, s.x_dot, s.x_ddot);
            
            z = y_dot - y_c;
            z_dot = this.K.*(yd - y) + this.D.*(yd_dot - z) + yd_ddot + z_c;
            y_ddot = (z_dot + yc_dot);

        end
        
    end
    
    properties (Access = public)
        
        traj_sc % object of type @TrajScale
        
        %% weights
        W % num_DoF x num_Kernels matrix where each row contrains the weights for each DoF

        %% impedance params
        K % num_DoF x num_DoF stiffness matrix
        D % num_DoF x num_DoF damping matrix 

    end
    
    
    properties (Access = {?GMP_Update, ?gmp_})
        
        Y0 % initial position
        Yg % target position
        
        Y0d % trained initial position
        Ygd % trained target position
        
        W0 % initial weights
        
        %% output state
        %  (deprecated)
        %  Used with @update, @getY, @getZ.
        y_dot % position derivative
        z_dot % scaled velocity derivative
        
    end
    
end
