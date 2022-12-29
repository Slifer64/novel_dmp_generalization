%% DMPo class
%  For encoding Cartesian Orientation.
%


classdef DMPo < matlab.mixin.Copyable
    
    methods  (Access = public)
        %% DMPo constructor.
        %  @param[in] dmp_type: Type on DMP (see @DMP_TYPE).
        %  @param[in] a_z: Parameter 'a_z' relating to the spring-damper system.
        %  @param[in] b_z: Parameter 'b_z' relating to the spring-damper system.
        %  @param[in] can_clock_ptr: Pointer to a DMP canonical system object.
        %  @param[in] shape_attr_gating_ptr: Pointer to a DMP gating function object.
        function this = DMPo(dmp_type, N_kernels, a_z, b_z, can_clock_ptr, shape_attr_gating_ptr)

            if (nargin < 5), can_clock_ptr = CanonicalClock(); end
            if (nargin < 6), shape_attr_gating_ptr = SigmoidGatingFunction(1.0, 0.5); end
            
            if (isscalar(N_kernels)), N_kernels = ones(3,1)*N_kernels; end
            if (isscalar(a_z)), a_z = ones(3,1)*a_z; end
            if (isscalar(b_z)), b_z = ones(3,1)*b_z; end

            this.can_clock_ptr = can_clock_ptr;
            this.shape_attr_gating_ptr = shape_attr_gating_ptr;
            
            this.Q0 = [1 0 0 0]';
            
            if (dmp_type == DMP_TYPE.STD)
                for i=1:3
                    this.dmp{i} = DMP(N_kernels(i), a_z(i), b_z(i), can_clock_ptr, shape_attr_gating_ptr);
                end
            elseif (dmp_type == DMP_TYPE.BIO)
                for i=1:3
                    this.dmp{i} = DMP_bio(N_kernels(i), a_z(i), b_z(i), can_clock_ptr, shape_attr_gating_ptr);
                end
            else
                error('[DMPo::DMPo]: Unsupported DMP type!');
            end
            
            for i=1:length(this.dmp), this.dmp{i}.setY0(0); end

        end
        
        
        %% Trains the DMPo.
        %  @param[in] train_method: The training method (see dmp_::TRAIN_METHOD enum).
        %  @param[in] Time: 1xN row vector with the timestamps of the training data points.
        %  @param[in] Quat_data: 3xN matrix with desired position at each column for each timestep.
        %  @param[in] rotVel_data: 3xN matrix with desired velocity at each column for each timestep.
        %  @param[in] rotAccel_data: 3xN matrix with desired acceleration at each column for each timestep.
        %  @param[in] train_error: Optinal pointer to return the training error as norm(F-Fd)/n_data.
        function train_err = train(this, train_method, Time, Quat_data, rotVel_data, rotAccel_data)

            tau = Time(end);
            this.setTau(tau);

            n_data = length(Time);
            q_data = zeros(3, n_data);
            qdot_data = zeros(3, n_data);
            qddot_data = zeros(3, n_data);

            Q0 = Quat_data(:,1);
            
            this.setQ0(Q0);

            for j=1:n_data
               Q1 = DMPo.quatTf(Quat_data(:,j), Q0);
               q_data(:,j) = quatLog(Q1);
               qdot_data(:,j) = DMPo.rotVel2qdot(rotVel_data(:,j), Q1);
               % qddot_data(:,j) = DMPo.rotAccel2qddot(rotAccel_data(:,j), rotVel_data(:,j), Q1);
            end
            
            Ts = Time(2) - Time(1);
            for i=1:3
                qddot_data(i,:) = [diff(qdot_data(i,:)) 0]/Ts;
            end
            
            

            if (nargout > 0)
                train_err = zeros(3,1);
                for i=1:3
                    train_err(i) = this.dmp{i}.train(train_method, Time, q_data(i,:), qdot_data(i,:), qddot_data(i,:));
                end
            else
                for i=1:3
                    this.dmp{i}.train(train_method, Time, q_data(i,:), qdot_data(i,:), qddot_data(i,:));
                end
            end

        end
      

        %% Sets the initial orientation.
        %  @param[in] Q0: Initial orientation (as unit quaternion).
        function setQ0(this, Q0)
            
            this.Q0 = Q0; 
            % for i=1:length(this.dmp) this.dmp{i}.setY0(0);
            
        end

        
        %% Calculates the derivatives of the DMP states. 
        %  The derivatives can then be retrieved with 'getXdot', 'getYdot' and 'getZdot'.
        %  @param[in] x: phase variable.
        %  @param[in] Y: 'y' state of the DMP: y=log(Q*Q0^{-1}).
        %  @param[in] Z: 'z' state of the DMP (the scaled ydot state).
        %  @param[in] G: 'g' goal/target of the DMP: y=log(Qg*Q0^{-1}).
        %  @param[in] Yc: Coupling term for the deonamical equation of the 'y' state.
        %  @param[in] Zc: Coupling term for the deonamical equation of the 'z' state.
        function update(this, x, Y, Z, G, Yc, Zc)

            if (nargin < 6), Yc = zeros(3,1); end
            if (nargin < 7), Zc = zeros(3,1); end

            if (isscalar(Yc)), Yc = ones(3,1)*Yc; end
            if (isscalar(Zc)), Zc = ones(3,1)*Zc; end

            for i=1:3, this.dmp{i}.update(x, Y(i), Z(i), G(i), Yc(i), Zc(i)); end

            this.dY = zeros(3,1);
            this.dZ = zeros(3,1);
            for i=1:3
                this.dY(i) = this.dmp{i}.getYdot();
                this.dZ(i) = this.dmp{i}.getZdot();
            end
            this.dx = this.phaseDot(x);

        end
        
        
        %% Returns the time derivative of the DMP's phase variable.
        %  Call @update first!
        %  @return: time derivative of the phase variable.
        function x_dot = getXdot(this), x_dot=this.dx; end
        
        
        %% Returns the time derivative of the DMP's 'y' state.
        %  Call @update first!
        %  @return: time derivative of 'y' state.
        function dy = getYdot(this), dy=this.dY; end

        
        %% Returns the time derivative of the DMP's 'z' state.
        %  Call @update first!
        %  @return: time derivative of 'z' state.
        function z_dot = getZdot(this), z_dot=this.dZ; end

        
        %% Returns the second derivative of the DMP's 'y' state.
        %  Call @update first!
        %  @param[in] tau_dot: time derivative of time scaling (optional, default=0).
        %  @param[in] yc_dot: time derivative of 'yc' coupling term (optional, default=0).
        %  @return: second time derivative of 'y' state.
        function y_ddot = getYddot(this, tau_dot, yc_dot)

            if (nargin < 2), tau_dot=0; end
            if (nargin < 3), yc_dot=0; end
            y_ddot = (this.getZdot() - tau_dot*this.getYdot() + yc_dot)/this.getTau(); 

        end
        
        
        function y_ddot = calcYddot(this, x, y, y_dot, yg, tau_dot, yc, zc, yc_dot)

            if (nargin < 6), tau_dot = 0; end
            if (nargin < 7), yc = zeros(3,1); end
            if (nargin < 8), zc = zeros(3,1); end
            if (nargin < 9), yc_dot = zeros(3,1); end
            
            n_dim = length(this.dmp);
            y_ddot = zeros(n_dim,1);
            for i=1:n_dim
                y_ddot(i) = this.dmp{i}.calcYddot(x, y(i), y_dot(i), yg(i), tau_dot, yc(i), zc(i), yc_dot(i));
            end

        end

        %% Returns the rotational velocity.
        %  Call @update first!
        %  @param[in] Q: the current orientation.
        %  @return: the rotational velocity.
        function rotVel = getRotVel(this, Q)

            Q1 = DMPo.quatTf(Q, this.Q0);
            qdot = this.getYdot(); 
            rotVel = DMPo.qdot2rotVel(qdot, Q1);

        end
        
        
        %% Returns the rotational acceleration.
        %  Call @update first!
        %  @param[in] Q: the current orientation.
        %  @param[in] tau_dot: time derivative of time scaling (optional, default=0).
        %  @param[in] yc_dot: time derivative of 'yc' coupling term (optional, default=0).
        %  @return: the rotational acceleration.
        function rotAccel = getRotAccel(this, Q, tau_dot, yc_dot)
            
            if (nargin < 3), tau_dot=0; end
            if (nargin < 4), yc_dot=0; end

            Q1 = DMPo.quatTf(Q, this.Q0);
            qddot = this.getYddot(tau_dot, yc_dot);
            rotVel = this.getRotVel(Q);
            rotAccel = DMPo.qddot2rotAccel(qddot, rotVel, Q1);

        end
        
        
        %% Calculates the rotational acceleration based on the current input variables.
        %  @param[in] x: phase variable.
        %  @param[in] Q: the current orientation.
        %  @param[in] rotVel: the rotational velocity.
        %  @param[in] Qg: the target orientation.
        %  @param[in] tau_dot: time derivative of time scaling (optional, default=0).
        %  @param[in] yc: Coupling term fo 'y' state diff-equation (optional, default=0).
        %  @param[in] zc: Coupling term fo 'z' state diff-equation (optional, default=0).
        %  @param[in] yc_dot: time derivative of 'yc' coupling term (optional, default=0).
        %  @return: the rotational acceleration.
        function rotAccel = calcRotAccel(this, x, Q, rotVel, Qg, tau_dot, yc, zc, yc_dot)
            
            if (nargin < 6), tau_dot = 0; end
            if (nargin < 7), yc = zeros(3,1); end
            if (nargin < 8), zc = zeros(3,1); end
            if (nargin < 9), yc_dot = zeros(3,1); end

            a_z = [this.dmp{1}.a_z; this.dmp{2}.a_z; this.dmp{3}.a_z];
            b_z = [this.dmp{1}.b_z; this.dmp{2}.b_z; this.dmp{3}.b_z];
            tau = this.getTau();

            qg = DMPo.quat2q(Qg, this.Q0);

            Q1 = DMPo.quatTf(Q, this.Q0);
            q = DMPo.quat2q(Q, this.Q0);
            invQ1 = quatInv(Q1);

            JQq = DMPo.jacobQq(Q1);
            JQq_dot = DMPo.jacobDotQq(Q1, rotVel);

            fo = zeros(3,1);
            for i=1:3, fo(i) = this.dmp{i}.shapeAttractor(x, qg(i)); end

            qdot = DMPo.rotVel2qdot(rotVel, Q1);
            qddot = (a_z.*b_z.*(qg - q) - tau*a_z.*qdot + a_z.*yc + fo + tau*yc_dot - tau*tau_dot*qdot + zc)/tau^2;
            
            rotAccel = 2*quatProd(JQq_dot*qdot + JQq*qddot, invQ1);
            rotAccel = rotAccel(2:4);

        end

        
        %% Returns the 'y' state of the DMP based on the current orientation.
        %  @param[in] Q: Current orientation (as unit quaternion).
        function y = getY(this, Q), y = DMPo.quat2q(Q, this.Q0); end

        
        %% Returns the 'z' state of the DMP based on the current rotational velocity and orientation.
        %  @param[in] rotVel: Current rotational velocity.
        %  @param[in] Q: Current orientation (as unit quaternion).
        function z = getZ(this, rotVel, Q)

            Q1 = DMPo.quatTf(Q, this.Q0);
            qdot = DMPo.rotVel2qdot(rotVel, Q1);
            dy = qdot;

            z = this.getTau()*dy; 

        end

        
        %% Returns the time scaling factor.
        %  @return: The time scaling factor.
        function tau = getTau(this), tau = this.can_clock_ptr.getTau(); end
        
        
        %% Sets the time scaling factor.
        %  @param[in] tau: The time scaling factor.
        function setTau(this, tau), this.can_clock_ptr.setTau(tau); end
        
        
        %% Returns the phase variable corresponding to the given time instant.
        %  @param[in] t: The time instant.
        %  @return: The phase variable for time 't'.
        function x = phase(this, t), x = this.can_clock_ptr.getPhase(t); end
        
        
        %% Returns the derivative of the phase variable.
        %  @param[in] x: The phase variable.
        %  @return: The derivative of the phase variable.
        function dx = phaseDot(this, x), dx = this.can_clock_ptr.getPhaseDot(x); end

        
        function J = getAcellPartDev_qg_tau(this, t, Y, dY, Y0, x, Yg, tau)
            
            n_dim = length(this.dmp);
            J = zeros(n_dim, n_dim+1);
            
            for i=1:n_dim
                C = this.dmp{i}.getAcellPartDev_g_tau(t, Y(i), dY(i), Y0(i), x, Yg(i), tau);
                J(i,i) = C(1);
                J(i,end) = C(2);
            end
            
        end
        
    end
    
    methods (Static)
        
        %% Expresses a given quaternion w.r.t. the initial orientation. 
        %  @param[in] Q: Orientation as unit quaternion.
        %  @param[in] Q0: Initial quaternion.
        %  @return: Q1 = Orientation w.r.t. Q0, i.e. Q1 = Q*Q0^{-1}.
        function Q1 = quatTf(Q, Q0), Q1 = quatProd(Q, quatInv(Q0)); end

        
        %% Returns the log of a given orientation w.r.t. the initial orientation.
        %  @param[in] Q: Orientation as unit quaternion.
        %  @param[in] Q0: Initial quaternion.
        %  @return: The logarithm of the Q w.r.t. Q0, i.e. q = log(Q*Q0^{-1}).
        function q = quat2q(Q, Q0), q = quatLog(DMPo.quatTf(Q,Q0)); end
        
        
        %% Returns the quaternion Q given the initial orientation Q0 and the log of Q w.r.t. to Q0.
        %  @param[in] q: Logarithm of orientation w.r.t. the initial orientation.
        %  @param[in] Q0: Initial orientation.
        %  @return: The orientation corresponding to log, i.e. Q = exp(q)*Q0
        function Q = q2quat(q, Q0), Q = quatProd( quatExp(q), Q0); end
        

        %% Returns derivative of log given the rotational velocity and orientation (expressed w.r.t. the initial orientation)
        %  @param[in] rotVel: Rotational velocity.
        %  @param[in] Q1: Orientation expressed w.r.t. the initial orientation.
        %  @return: Derivative of log.
        function qdot = rotVel2qdot(rotVel, Q1)
            
            JqQ = DMPo.jacobqQ(Q1);
            qdot = 0.5*JqQ * quatProd([0; rotVel], Q1);

        end
        
        
        function rotVel = qdot2rotVel(qdot, Q1)
            
            JQq = DMPo.jacobQq(Q1);
            rotVel = 2 * quatProd( JQq*qdot, quatInv(Q1) );
            rotVel = rotVel(2:4);

        end
        

        function qddot = rotAccel2qddot(rotAccel, rotVel, Q1)
            
            rotVelQ = [0; rotVel];
            rotAccelQ = [0; rotAccel];

            J = DMPo.jacobqQ(Q1);
            Jdot = DMPo.jacobDotqQ(Q1, rotVel);

            qddot = 0.5 * (Jdot * quatProd(rotVelQ, Q1) + J * quatProd( rotAccelQ+0.5*quatProd(rotVelQ,rotVelQ), Q1 ) );

        end
        

        function rotAccel = qddot2rotAccel(qddot, rotVel, Q1)

            qdot = DMPo.rotVel2qdot(rotVel, Q1);
            invQ1 = quatInv(Q1);
            J = DMPo.jacobQq(Q1);
            Jdot = DMPo.jacobDotQq(Q1, rotVel);

            rotAccel = 2 * ( quatProd( Jdot*qdot+J*qddot, invQ1 ) );
            rotAccel = rotAccel(2:4);

        end
        
        
        %% Returns the Jacobian from the derivative of log to the derivative of Q.
        %  @param[in] Q1: The orientation w.r.t. the initial orientation.
        %  @return: Jacobian.
        function JQq = jacobQq(Q1)

            if ( (1-abs(Q1(1))) <= DMPo.zero_tol)
                JQq = [zeros(1, 3); eye(3,3)];
                return;
            end

            w = Q1(1);
            v = Q1(2:4);
            norm_v = norm(v);
            eta = v / norm_v;
            s_th = norm_v;
            c_th = w;
            th = atan2(s_th, c_th);
            Eta = eta*eta';

            JQq = zeros(4,3);
            JQq(1,:) = -0.5 * s_th * eta';
            JQq(2:4,:) = 0.5 * ( (eye(3,3) - Eta)*s_th/th + c_th*Eta );

        end
        
        
        %% Returns the Jacobian from the derivative of Q to the derivative of log.
        %  @param[in] Q1: The orientation w.r.t. the initial orientation.
        %  @return: Jacobian.
        function JqQ = jacobqQ(Q1)
            
            if ( (1-abs(Q1(1))) <= DMPo.zero_tol)
                JqQ = [zeros(3,1) eye(3,3)];
                return;
            end

            w = Q1(1);
            v = Q1(2:4);
            norm_v = norm(v);
            eta = v / norm_v;
            s_th = norm_v;
            c_th = w;
            th = atan2(s_th, c_th);

            JqQ = zeros(3,4);
            JqQ(:,1) = 2*eta*(th*c_th - s_th)/s_th^2;
            JqQ(:,2:4) = 2*eye(3,3)*th/s_th;

        end
       
        
        %% Returns the time derivative of the Jacobian from the derivative of log to the derivative of Q.
        %  @param[in] Q1: The orientation w.r.t. the initial orientation.
        %  @return: Jacobian time derivative.
        function JqQ_dot = jacobDotqQ(Q1, rotVel)

            qdot = DMPo.rotVel2qdot(rotVel, Q1);

            if ( (1-abs(Q1(1))) <= DMPo.zero_tol)
                JqQ_dot = [-qdot/3 zeros(3,3)];
                return;
            end

            w = Q1(1);
            v = Q1(2:4);
            norm_v = norm(v);
            eta = v / norm_v;
            s_th = norm_v;
            c_th = w;
            th = atan2(s_th, c_th);
            Eta = eta*eta';
            temp = (th*c_th-s_th)/s_th^2;

            JqQ_dot = zeros(3,4);
            JqQ_dot(:,1) = ((-th/s_th - 2*c_th*temp/s_th)*Eta + temp*(eye(3,3)-Eta)/th)*qdot;
            JqQ_dot(:,2:4) = -temp*dot(eta,qdot)*eye(3,3);

        end
        
        
        %% Returns the time derivative of the Jacobian from the derivative of Q to the derivative of log.
        %  @param[in] Q1: The orientation w.r.t. the initial orientation.
        %  @return: Jacobian time derivative.
        function JQq_dot = jacobDotQq(Q1, rotVel)
            
            qdot = DMPo.rotVel2qdot(rotVel, Q1);

            if ( (1-abs(Q1(1))) <= DMPo.zero_tol)
                JQq_dot = [-qdot'/4; zeros(3,3)];
                return;
            end

            w = Q1(1);
            v = Q1(2:4);
            norm_v = norm(v);
            eta = v / norm_v;
            s_th = norm_v;
            c_th = w;
            th = atan2(s_th, c_th);
            Eta = eta*eta';
            I_eta = eye(3,3) - Eta;
            temp = ((th*c_th-s_th)/th^2);

            JQq_dot = zeros(4,3);
            JQq_dot(1,:) = -0.25 * qdot' * (c_th*Eta + (s_th/th)*I_eta);
            JQq_dot(2:4,:) = 0.25*dot(eta,qdot)*( temp*I_eta - s_th*Eta ) + 0.25*temp*( eta*(qdot'*I_eta) + (I_eta*qdot)*eta' );

        end
          
    end
    
    
    properties (Constant)  
        
          zero_tol = 1e-15;
          
    end
    
    properties  (Access = public)
        
        dmp % 3x1 vector of DMP producing the desired trajectory for each x, y and z axes of the orientation error.

        can_clock_ptr % handle (pointer) to the canonical clock
        shape_attr_gating_ptr % pointer to gating function for the shape attractor
        
    end
       
    properties  (Access = protected)
        
        Q0 % initial orientation
        
        dZ % second derivative of the orientation error
        dY % first derivative of the orientation error
        dx % pahse variable derivative
        
    end
    
end
