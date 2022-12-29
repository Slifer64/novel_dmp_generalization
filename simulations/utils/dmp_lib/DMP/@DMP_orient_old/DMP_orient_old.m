%% DMP_orient_old class
%  For encoding Cartesian Orientation.
%


classdef DMP_orient_old < handle % : public DMP_
    
    properties
        can_clock_ptr % handle (pointer) to the canonical clock
        shape_attr_gating_ptr % pointer to gating function for the shape attractor
        dmp
    end

    methods
        %% DMP constructor.
        %  @param[in] N_kernels: the number of kernels
        %  @param[in] a_z: Parameter 'a_z' relating to the spring-damper system.
        %  @param[in] b_z: Parameter 'b_z' relating to the spring-damper system.
        %  @param[in] can_clock_ptr: Pointer to a DMP canonical system object.
        function this = DMP_orient_old(dmp_type, N_kernels, a_z, b_z, can_clock_ptr, shape_attr_gating_ptr)
            
            if (nargin < 5), can_clock_ptr = CanonicalClock(); end
            if (nargin < 6), shape_attr_gating_ptr = SigmoidGatingFunction(1.0, 0.5); end
            
            if (isscalar(N_kernels)), N_kernels = ones(3,1)*N_kernels; end
            if (isscalar(a_z)), a_z = ones(3,1)*a_z; end
            if (isscalar(b_z)), b_z = ones(3,1)*b_z; end
            
            this.can_clock_ptr = can_clock_ptr;
            this.shape_attr_gating_ptr = shape_attr_gating_ptr;
            
            if (dmp_type == DMP_TYPE.STD)
                for i=1:3
                    this.dmp{i} = DMP(N_kernels(i), a_z(i), b_z(i), can_clock_ptr, shape_attr_gating_ptr);
                end
            elseif (dmp_type == DMP_TYPE.BIO)
                for i=1:3
                    this.dmp{i} = DMP_bio(N_kernels(i), a_z(i), b_z(i), can_clock_ptr, shape_attr_gating_ptr);
                end
            else
                error('[DMP_orient_old::DMP_orient_old]: Unsupported DMP type!');
            end
  
        end

        %% Trains the DMP_orient_old.
        %  @param[in] Time: Row vector with the timestamps of the training data points.
        %  @param[in] yd_data: Row vector with the desired potition.
        %  @param[in] dyd_data: Row vector with the desired velocity.
        %  @param[in] ddyd_data: Row vector with the desired accelaration.
        function [train_error] = train(this, train_method, Time, Qd_data, vRot_data, dvRot_data)

            n_data = length(Time);

            tau = Time(end);
            this.setTau(tau);
            
            Qgd = Qd_data(:,end);
            
            Y_data = zeros(3, n_data);
            for j=1:n_data, Y_data(:,j) = -this.quatLog( this.quatProd( Qgd, this.quatInv(Qd_data(:,j)) ) ); end
            
            train_error = zeros(3,1);
            for i=1:3
                train_error(i) = this.dmp{i}.train(train_method, Time, Y_data(i,:), vRot_data(i,:), dvRot_data(i,:));  
            end

        end
        
        
        %% Returns the derivatives of the DMP states.
        %  @param[in] x: phase variable.
        %  @param[in] y: \a y state of the this.
        %  @param[in] z: \a z state of the this.
        %  @param[in] y0: initial position.
        %  @param[in] g: Goal position.
        %  @param[in] y_c: coupling term for the dynamical equation of the \a y state.
        %  @param[in] z_c: coupling term for the dynamical equation of the \a z state.
        %  @param[out] dy: derivative of the \a y state of the this.
        %  @param[out] dz: derivative of the \a z state of the this.
        %  @param[out] dx: derivative of the phase variable of the this.
        function dvRot = calcRotAccel(this, x, Q, vRot, Qg, Q0, y_c, z_c)
            
            if (nargin < 7), y_c = zeros(3,1); end
            if (nargin < 8), z_c = zeros(3,1); end
            
            tau = this.getTau();
            
            Y = -this.quatLog( this.quatProd(Qg, this.quatInv(Q)) );
            Y0 = -this.quatLog( this.quatProd(Qg, this.quatInv(Q0)) );
            Yg = zeros(3,1);
            
            deta = zeros(3,1);
            eta = vRot * tau;
            
            dy = zeros(3,1);
            
            for i=1:3
               this.dmp{i}.setY0(Y0(i));
               this.dmp{i}.update(x, Y(i), eta(i), Yg(i), y_c(i), z_c(i));
               deta(i) = this.dmp{i}.getZdot();
            end
            
            dvRot = deta/tau;
            
            dx = this.phaseDot(x);
            
            % if (nargout>1), dQ = 0.5 * quatProd([0; eta], Q) / tau; end
            % if (nargout>2), dx = this.phaseDot(x); end

        end
        
        
        %% Returns the time scale of the DMP_orient_old.
        %  @param[out] tau: The time scale of the this.
        function tau = getTau(this)

            tau = this.can_clock_ptr.getTau();

        end
        
        
        %% Sets the time scale of the DMP_orient_old.
        function tau = setTau(this, tau)

            this.can_clock_ptr.setTau(tau);

        end
        
        
        %% Returns the phase variable.
        %  @param[in] t: The time instant.
        %  @param[out] x: The phase variable for time 't'.
        function x = phase(this, t)
            
            x = this.can_clock_ptr.getPhase(t);

        end
        
        
        %% Returns the derivative of the phase variable.
        %  @param[in] x: The phase variable.
        %  @param[out] dx: The derivative of the phase variable.
        function dx = phaseDot(this, x)
            
            dx = this.can_clock_ptr.getPhaseDot(x);

        end


        %% Returns the goal attractor of the this.
        %  @param[in] x: The phase variable.
        %  @param[in] y: \a y state of the this.
        %  @param[in] z: \a z state of the this.
        %  @param[in] g: Goal position.
        %  @param[out] goal_attr: The goal attractor of the this.
        function goal_attr = goalAttractor(this, x, Q, eta, Qg)

            goal_attr = zeros(3,1);
            
            Y = -quatLog( quatProd(Qg, quatInv(Q)) );
            Yg = zeros(3,1);
            
            for i=1:3
                goal_attr(i) = this.dmp{i}.goalAttractor(x, Y(i), eta(i), Yg(i));
            end

        end
        
        function shape_attr = shapeAttractor(this, x, Qg, Q0)

            shape_attr = zeros(3,1);
            
            Y0 = -quatLog( quatProd(Qg, quatInv(Q0)) );
            Yg = zeros(3,1);
            
            for i=1:3
                shape_attr(i) = this.dmp{i}.shapeAttractor(x, Y0(i), Yg(i));
            end

        end
     
        
        %% Returns the shape attractor gating factor.
        %  @param[in] x: The phase variable.
        function sAttrGat = shapeAttrGating(this, x)

            sAttrGat = this.shape_attr_gating_ptr.getOutput(x);
            sAttrGat(sAttrGat<0) = 0.0;

        end
        
        
        %% Returns the goal attractor gating factor.
        %  @param[in] x: The phase variable.
        function gAttrGat = goalAttrGating(this, x)
            
            gAttrGat = 1.0;

        end


    end
    
    methods (Static)
        
        function v_rot = quatLog(Q)

            n = Q(1);
            e = Q(2:4);
            norm_e = norm(e);

            if (norm_e > 1e-16)
                theta = 2*real(atan2(norm_e,n));
                v_rot = theta*e/norm_e;
            else
                v_rot = zeros(size(e));
            end

        end
        
        function Q = quatExp(v_rot)

            norm_v_rot = norm(v_rot);
            theta = norm_v_rot;

            if (norm_v_rot > 1e-16)
                Q(1) = cos(theta/2);
                Q(2:4) = sin(theta/2)*v_rot/norm_v_rot;
            else
                Q = [1 0 0 0]';
            end

        end
        
        function Q12 = quatProd(Q1, Q2)

            Q1 = Q1(:);
            Q2 = Q2(:);

            n1 = Q1(1);
            e1 = Q1(2:4);

            n2 = Q2(1);
            e2 = Q2(2:4);

            Q12 = zeros(4,1);
            Q12(1) = n1*n2 - e1'*e2;
            Q12(2:4) = n1*e2 + n2*e1 + cross(e1,e2);

        end
        
        function invQ = quatInv(Q)

            invQ = zeros(size(Q));

            invQ(1) = Q(1);
            invQ(2:4) = -Q(2:4);

        end

    end
    
end
