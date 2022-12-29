%% DMP_orient class
%  For encoding Cartesian Orientation.
%


classdef DMP_orient < handle

    properties
        
        Q0
        Qgd % Trained target orientation as unit quaternion.
        Q0d % Trained initial orientation as unit quaternion.
        log_Qgd_invQ0d % log(Qgd * inv(Q0d)), precalculated
        tau_d % Trained time-scaling.
        
        eq_f % 3x1 vector of WSoG producing the desired orientation.
        vRot_f % 3x1 vector of WSoG producing the desired rotational velocity.
        dvRot_f % 3x1 vector of WSoG producing the desired rotational acceleration.
        
        a_z % parameter 'a_z' relating to the spring-damper system
        b_z % parameter 'b_z' relating to the spring-damper system

        can_clock_ptr % handle (pointer) to the canonical clock
        shape_attr_gating_ptr % pointer to gating function for the shape attractor

        zero_tol % tolerance value used to avoid divisions with very small numbers
        
        dphi % derivative of 'phi' state
        omega % rotational velocity
        dx % phase variable derivative

    end

    methods
        %% DMP_orient constructor.
        %  @param[in] N_kernels: 3x3 matrix where each row contains the kernels for position,
        %                       velocity and acceleration and each row for x, y and z coordinate.
        %                       If 'N_kernels' is a column vector then every coordinate will have the
        %                       same number of kernels.
        %  @param[in] a_z: Parameter 'a_z' relating to the spring-damper system.
        %  @param[in] b_z: Parameter 'b_z' relating to the spring-damper system.
        %  @param[in] can_clock_ptr: Pointer to a DMP canonical system object.
        %  @param[in] shape_attr_gating_ptr: Pointer to a DMP gating function object.
        %
        function this = DMP_orient(N_kernels, a_z, b_z, can_clock_ptr, shape_attr_gating_ptr)
            
            if (nargin < 5), can_clock_ptr = CanonicalClock(); end
            if (nargin < 6), shape_attr_gating_ptr = SigmoidGatingFunction(1.0, 0.5); end
            if (isscalar(N_kernels)), N_kernels = ones(3,1)*N_kernels; end
            if (size(N_kernels,2) == 1), N_kernels = repmat(N_kernels,1,3); end
            if (isscalar(a_z)), a_z = ones(3,1)*a_z; end
            if (isscalar(b_z)), b_z = ones(3,1)*b_z; end
            
            this.zero_tol = 1e-30; %realmin;
            this.a_z = a_z;
            this.b_z = b_z;
            this.can_clock_ptr = can_clock_ptr;
            this.shape_attr_gating_ptr = shape_attr_gating_ptr;

            shapeAttrGatFun_ptr = @(x)getOutput(this.shape_attr_gating_ptr,x);
            for i=1:3
                this.eq_f{i} = WSoG(N_kernels(1,i), shapeAttrGatFun_ptr);
                this.vRot_f{i} = WSoG(N_kernels(2,i), shapeAttrGatFun_ptr);
                this.dvRot_f{i} = WSoG(N_kernels(3,i), shapeAttrGatFun_ptr);
            end

        end
        
        
        %% Trains the DMP_orient.
        %  @param[in] train_method: The training method (see dmp_::TRAIN_METHOD enum).
        %  @param[in] Time: 1xN row vector with the timestamps of the training data points.
        %  @param[in] Qd_data: 4xN matrix with desired orientation as unit quaternion at each column for each timestep.
        %  @param[in] vRotd_data: 3xN matrix with desired angular velocity at each column for each timestep.
        %  @param[in] dvRotd_data: 3xN matrix with desired angular acceleration at each column for each timestep.
        %  @param[in] train_error: Optinal pointer to return the training error as norm(F-Fd)/n_data.
        %
        function train_err = train(this, train_method, Time, Qd_data, vRotd_data, dvRotd_data)
            
            n_data = length(Time);
            
            this.Q0d = Qd_data(:,1);
            this.Qgd = Qd_data(:,end);
            this.log_Qgd_invQ0d = quatLog( quatProd(this.Qgd, quatInv(this.Q0d) ) );
            
            eqd = zeros(3, n_data);           
            for j=1:n_data, eqd(:,j) = quatLog(quatProd(Qd_data(:,j),quatInv(this.Qgd))); end

            tau = Time(end);
            this.tau_d = tau;
            this.setTau(tau);
    
            x = this.phase(Time);
                
            if (nargout > 0)
                train_err = zeros(3,3);
                for i=1:3
                    train_err(1,i) = this.eq_f{i}.train(train_method, x, eqd(i,:));
                    train_err(2,i) = this.vRot_f{i}.train(train_method, x, vRotd_data(i,:));
                    train_err(3,i) = this.dvRot_f{i}.train(train_method, x, dvRotd_data(i,:));
                end
            else
                for i=1:3
                    this.eq_f{i}.train(train_method, x, eqd(i,:));
                    this.vRot_f{i}.train(train_method, x, vRotd_data(i,:));
                    this.dvRot_f{i}.train(train_method, x, dvRotd_data(i,:));
                end
            end

        end
        
        function setQ0(this, Q0), this.Q0 = Q0; end
        
        %% Returns the angular acceleration for the given input state defined by the timestamp,
        %  the orientation, the angular velocity and acceleration, the initial and target orientation
        %  and an optinal coupling term.
        %  @param[in] x: phase variable.
        %  @param[in] Q: Current orientation as unit quaternion.
        %  @param[in] vRot: Current rotational velocity.
        %  @param[in] Q0: initial orientation as unit quaternion.
        %  @param[in] Qg: Goal orientation as unit quaternion.
        %  @param[in] Z_c: Coupling term. (optional, default=arma::vec().zeros(3))
        %  @return dvRot: Rotational acceleration.
        %
        function dvRot = calcRotAccel(this, x, Q, vRot, Qg, Z_c)

            if (nargin < 7), Z_c=zeros(3,1); end
            
            tau = this.getTau();
            phi = vRot*tau;
            update(this, x, Q, phi, Qg, zeros(3,1), Z_c);
            dvRot = this.getDphi()/tau;
            
        end
        
        function update(this, x, Q, phi, Qg, Y_c, Z_c)

            if (nargin < 7), Y_c=zeros(3,1); end
            if (nargin < 8), Z_c=zeros(3,1); end

            eqd = zeros(3,1);
            vRotd = zeros(3,1);
            dvRotd = zeros(3,1);
            
            for i=1:3
                eqd(i) = this.eq_f{i}.output(x);
                vRotd(i) = this.vRot_f{i}.output(x);
                dvRotd(i) = this.dvRot_f{i}.output(x);
            end 

            Qd = quatProd( quatExp(eqd), this.Qgd );
            
            Q0 = this.Q0;
            
            tau = this.getTau();
            tau_d = this.tau_d;
            Qgd = this.Qgd;
            Q0d = this.Q0d;
            % kt = this.tau_d / tau;
            
            ks = quatLog( quatProd(Qg, quatInv(Q0) ) ) ./ this.log_Qgd_invQ0d;
            % ks = ones(3,1);
            
            Qd_s = quatExp( ks.*quatLog(Qd) );
            Qgd_s = quatExp( ks.*quatLog(Qgd) );
            
%             Qrot = quatProd( quatProd(Qg,quatInv(Q0)) , quatInv( quatProd(Qgd,quatInv(Q0d)) ));
%             Qrot = quatProd( Qg, quatInv(Qgd));
%             
%             Qd = quatProd(Qrot, Qd);
%             vQ = quatProd( quatProd(Qrot, [0;vRotd]), quatInv(Qrot) );
%             vRotd = vQ(2:4);
%             dvQ = quatProd( quatProd(Qrot, [0;dvRotd]), quatInv(Qrot) );
%             dvRotd = dvQ(2:4);
% 
%             eo = quatLog( quatProd(Q, quatInv(Qd)) );
            
            %% y-g+ks*gd-ks*yd
%             eo = quatLog( quatProd( quatProd(Q, quatInv(Qg)), quatProd(Qgd_s, quatInv(Qd_s)) ) );
            
            %% ks * log( Q * Qg^{-1} * Qgd * Qd^{-1} )
%             eo = ks .* quatLog( quatProd( quatProd(Q,quatInv(Qg)), quatProd(Qgd,quatInv(Qd)) ) );

            %% (y - y0) - ks*(yd - y0d)
%             eo = quatLog( quatProd( quatProd(Q,quatInv(Q0)), quatInv(quatExp(ks.*quatLog(quatProd(Qd,quatInv(Q0d))))) ) );
            DQ = quatProd(Qd,quatInv(Q0d));
            DQs = quatExp(ks.*quatLog(DQ));
            Qds = quatProd(DQs,Q0);
            eo = quatLog( quatProd(Q, quatInv(Qds) ) );
%             eo = quatLog(quatProd(Q,quatInv(Qd)));
            
            %% (y-ks*yd) - (g-ks*gd)
%             QQd = quatProd( Q, quatInv(Qd_s) );     
%             QgQgd = quatProd( Qg, quatInv(Qgd_s) );
%             eo = quatLog( quatProd( QQd, quatInv(QgQgd) ) );
            
            %% (y-g) - ks*(yd-gd)
%             QQg = quatProd(Q,quatInv(Qg));
%             inv_exp_QdQgd = quatInv( quatExp( ks.* quatLog( quatProd(Qd,quatInv(this.Qgd)) ) ) );
%             eo = quatLog( quatProd(QQg, inv_exp_QdQgd) );

            %% -log(Qg*Q^{-1}) + log(Qgd*Qd^{-1}) 
%             eo = -quatLog(quatProd(Qg,quatInv(Q))) + ks.*quatLog(quatProd(Qgd,quatInv(Qd)));

%             eqd2 = zeros(3,1);    
%             for i=1:3, eqd(i) = this.eq_f{i}.output(x+0.01); end 
%             Qd2 = quatProd( quatExp(eqd2), this.Qgd );
%             
%             v = quatLog( quatProd(Qd2, quatInv(Qd)) );
%             v2 = quatLog( quatProd( quatExp(ks.*quatLog(quatProd(Qd2,quatInv(Q0d)))), quatExp(ks.*quatLog(quatProd(Qd,quatInv(Q0d)))) ) );
%             ks = v2./v;
   
            this.dphi = ( -(this.a_z.*this.b_z).*eo - this.a_z.*phi + tau_d^2*ks.*dvRotd + tau_d*this.a_z.*ks.*vRotd + Z_c ) / tau;
            this.omega = (phi + Y_c) / tau;
            this.dx = this.phaseDot(x);

        end
        
        function dx = getDx(this), dx=this.dx; end
        function dphi = getDphi(this), dphi=this.dphi; end
        function omega = getOmega(this), omega=this.omega; end

        %% Returns the time scaling factor.
        %  @return: The time scaling factor.
        %
        function tau = getTau(this)

            tau = this.can_clock_ptr.getTau();

        end
        
        
        %% Sets the time scaling factor.
        %  @param[in] tau: The time scaling factor.
        %
        function setTau(this, tau)

            this.can_clock_ptr.setTau(tau);

        end
        
        
        %% Returns the phase variable corresponding to the given time instant.
        %  @param[in] t: The time instant.
        %  @return: The phase variable for time 't'.
        %
        function x = phase(this, t)
            
            x = this.can_clock_ptr.getPhase(t);

        end
        
        
        %% Returns the derivative of the phase variable.
        %  @param[in] x: The phase variable.
        %  @return: The derivative of the phase variable.
        %
        function dx = phaseDot(this, x)
            
            dx = this.can_clock_ptr.getPhaseDot(x);

        end

    end
    
end
