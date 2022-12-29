%% GMPo class
%  GMP for encoding Cartesian Orientation.
%  Encodes the orientation by formulating a GMP in the quaternion log space.
%  The class provides the necessary tranformation betweeen the quaternion
%  log (and its 1st and 2nd time derivatives) and the quaternion,
%  rotational velocity and acceleration.
%


classdef GMPo < GMP
    
    methods  (Access = public)
        %% Constructor. Constructs a GMP defined in the quat-log space.
        %  @param[in] N_kernels: 3x1 vector with the number of kernels for each dim of the quat log.
        %  @param[in] kernels_std_scaling: Scaling for std of kernels (optional, default=2).
        %  \note: Each of the arguments 'N_kernels', 'D', 'K' can be scalar or a 3x1 vector.
        function this = GMPo(N_kernels, kernels_std_scaling)

            if (nargin < 1), N_kernels=2; end
            if (nargin < 2), kernels_std_scaling=1; end
            
            this = this@GMP(3, N_kernels, kernels_std_scaling);
            
            this.setScaleMethod(TrajScale_Rot_min());
            
            this.setQ0([1 0 0 0]');
            this.Qd0 = [1 0 0 0]';
            
            this.K = 10 * ones(3,1);
            this.D = 2 * ones(3,1);

        end
        
        %% ==============================
        %% =========  Training  =========
        %% ==============================
        
        %% Trains the GMPo.
        %  @param[in] train_method: the training method to use, as a string ('LWR', 'LS').
        %  @param[in] Time: Row vector with the timestamps of the training data points.
        %  @param[in] Quat_data: Matrix with the desired unit quaternion (4x1 vector) in each column.
        %  @param[out] train_error: The training error expressed as the mse error.
        function train_error = train(this, train_method, Time, Quat_data)

            n_data = length(Time);
            this.setQ0(Quat_data(:,1));
            
            this.Qd0 = Quat_data(:,1);
            
            qd_data = zeros(3, n_data);
            for j=1:n_data, qd_data(:,j) = GMPo.quat2q(Quat_data(:,j), this.Q0); end
             
            if (nargout > 0), train_error = train@GMP(this, train_method, Time, qd_data);
            else, train@GMP(this, train_method, Time, qd_data); end
            
        end
        
        
        %% Sets the initial orientation.
        %  @param[in] Q0: Initial orientation (as unit quaternion).
        function setQ0(this, Q0)
            
            this.Q0 = Q0;
            
        end
        
        function Q0 = getQ0(this)
            
            Q0 = this.Q0;
            
        end

        
        %% Sets the goal/target orientation.
        %  @param[in] Qg: Goal/target orientation (as unit quaternion).
        function setQg(this, Qg)
            
            qg = GMPo.quat2q(Qg, this.Q0);
            this.setGoal(qg);
            
        end

        %% ===========================================
        %% =========  GMP output trajectory  =========
        %% ===========================================
        
        function Qd = getQd(this, x)

            Qd = GMPo.q2quat(this.getYd(x), this.Q0);
            
        end
        
        function vd = getVd(this, x, x_dot)

            %Ts = 0.002;
            %Qd = this.getQd(x);
            %Qd2 = this.getQd(x+x_dot*Ts);
            %vd = gmp_.quatLog(quatDiff(Qd2,Qd)) / Ts;
            
            Q1 = gmp_.quatExp(this.getYd(x));
            vd = gmp_.qLogDot_to_rotVel(this.getYdDot(x,x_dot), Q1);
            
        end
        
        function vd_dot = getVdDot(this, x, x_dot, x_ddot)

            %Ts = 0.002;
            %vd = this.getVd(x, x_dot);
            %vd2 = this.getVd(x + x_dot*Ts, x_dot + x_ddot*Ts);
            %vd_dot = (vd2 - vd) / Ts;
            
            Q1 = gmp_.quatExp(this.getYd(x));
            vd = gmp_.qLogDot_to_rotVel(this.getYdDot(x,x_dot), Q1);
            vd_dot = gmp_.qLogDDot_to_rotAccel(this.getYdDDot(x,x_dot,x_ddot), vd, Q1);
            
        end
        
        function [Qd, vd, vd_dot] = getRefTraj(this, x, x_dot, x_ddot)
            
            yd = this.getYd(x);
            Q1 = gmp_.quatExp(this.getYd(x));
            
            Qd = GMPo.q2quat(yd, this.Q0);
            vd = gmp_.qLogDot_to_rotVel(this.getYdDot(x,x_dot), Q1);
            vd_dot = gmp_.qLogDDot_to_rotAccel(this.getYdDDot(x,x_dot,x_ddot), vd, Q1);
            
        end
        
        
        %% ===============================================
        %% ==========  original DMP functions ============
        %  Deprecated. Use @getYd, @getYdDot, @getYdDDot instead.
        
        %% Returns the rotational velocity.
        %  Call @update first!
        %  @param[in] Q: the current orientation.
        %  @return: the rotational velocity.
        function rotVel = getRotVel(this, Q)

            Q1 = GMPo.getQ1(Q, this.Q0);
            rotVel = gmp_.qLogDot_to_rotVel(this.getYdot(), Q1);

        end
        
        
        %% Returns the rotational acceleration.
        %  Call @update first!
        %  @param[in] Q: the current orientation.
        %  @param[in] yc_dot: time derivative of 'yc' coupling term (optional, default=0).
        %  @return: the rotational acceleration.
        function rotAccel = getRotAccel(this, Q, yc_dot)

            if (nargin < 3), yc_dot=0; end

            Q1 = GMPo.getQ1(Q, this.Q0);
            qddot = this.getYddot(yc_dot);
            rotVel = gmp_.qLogDot_to_rotVel(this.getYdot(), Q1);
            rotAccel = gmp_.qLogDDot_to_rotAccel(qddot, rotVel, Q1);

        end
        
        
        %% Calculates the rotational acceleration based on the current input variables.
        %  @param[in] s: Object of type @GMP_phase.
        %  @param[in] Q: the current orientation.
        %  @param[in] rotVel: the rotational velocity.
        %  @param[in] yc: Coupling term fo 'y' state diff-equation (optional, default=0).
        %  @param[in] zc: Coupling term fo 'z' state diff-equation (optional, default=0).
        %  @param[in] yc_dot: time derivative of 'yc' coupling term (optional, default=0).
        %  @return: the rotational acceleration.
        function rotAccel = calcRotAccel(this, s, Q, rotVel, y_c, z_c, yc_dot)
            
            if (nargin < 6), y_c = zeros(3,1); end
            if (nargin < 7), z_c = zeros(3,1); end
            if (nargin < 8), yc_dot = zeros(3,1); end

            Q1 = GMPo.getQ1(Q, this.Q0);
            q = gmp_.quatLog(Q1); % q = GMPo.quat2q(Q, this.Q0);
            invQ1 = gmp_.quatInv(Q1);

            JQq = gmp_.jacob_Q_qLog(Q1);
            JQq_dot = gmp_.jacobDot_Q_qLog(Q1, rotVel);

            qdot = gmp_.rotVel_to_qLogDot(rotVel, Q1);
            qddot = this.calcYddot(s, q, qdot, y_c, z_c, yc_dot);

            rotAccel = 2*gmp_.quatProd(JQq_dot*qdot + JQq*qddot, invQ1);
            rotAccel = rotAccel(2:4);

        end

        
        %% Returns the 'y' state of the DMP based on the current orientation.
        %  @param[in] Q: Current orientation (as unit quaternion).
        function y = getY(this, Q)
            
            y = GMPo.quat2q(Q, this.Q0); 
        
        end

        
        %% Returns the 'z' state of the DMP based on the current rotational velocity and orientation.
        %  @param[in] rotVel: Current rotational velocity.
        %  @param[in] Q: Current orientation (as unit quaternion).
        function z = getZ(this, rotVel, Q)

            Q1 = GMPo.getQ1(Q, this.Q0);
            qdot = gmp_.rotVel_to_qLogDot(rotVel, Q1);
            z = qdot; 

        end
 
    end
    
    methods (Static, Access = public)
          
        %% Expresses a given quaternion w.r.t. the initial orientation. 
        %  @param[in] Q: Orientation as unit quaternion.
        %  @param[in] Q0: Initial quaternion.
        %  @return: The orientation w.r.t. Q0, i.e. Q1 = Q*Q0^{-1}.
        function Q1 = getQ1(Q, Q0), Q1 = gmp_.quatDiff(Q,Q0); end

        
        %% Returns the log of a given orientation w.r.t. the initial orientation.
        %  @param[in] Q: Orientation as unit quaternion.
        %  @param[in] Q0: Initial quaternion.
        %  @return: The logarithm of the Q w.r.t. Q0, i.e. q = log(Q*Q0^{-1}).
        function q = quat2q(Q, Q0), q = gmp_.quatLog(GMPo.getQ1(Q,Q0)); end
        
        
        %% Returns the quaternion Q given the initial orientation Q0 and the log of Q w.r.t. Q0.
        %  @param[in] q: Logarithm of orientation w.r.t. the initial orientation.
        %  @param[in] Q0: Initial orientation.
        %  @return: The orientation corresponding to log, i.e. Q = exp(q)*Q0
        function Q = q2quat(q, Q0), Q = gmp_.quatProd( gmp_.quatExp(q), Q0); end

    end
    
    
    properties (Constant)  
        
        zero_tol = 1e-15;
          
    end
    
       
    properties  (Access = {?gmp_})
        
        Q0 % initial orientation (as unit quaternion)
        
        Qd0
        
    end
    
end
