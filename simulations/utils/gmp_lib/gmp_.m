%% GMP namespace
%   Contains general utility functions defined in the "gmp_" namespace
%

classdef gmp_
    
    %% ========================================
    %% ==============   MATH   ================
    %% ========================================    
    methods (Access = public, Static)

        %% quatProd
        %  Returns the quaternion product between two quaternions.
        %  @param[in] Q1: First quaterion as a 4x1 vector.
        %  @param[in] Q2: Second quaterion as a 4x1 vector.
        %  @return: the quaternion product Q1*Q2.
        %  \note All quaternions must have the form Q = [u; v] where u is
        %  the scalar and v the 3x1 vector part. The quaternions need not
        %  be unit.
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

            % if Q12(1)<0, Q12=-Q12; end

        end
        
        %% quatInv
        %  Returns the inverse of a quaternion.
        %  @param[in] Q: Quaterion as a 4x1 vector [u; v], with u the scalar and v the vector part.
        %  @return: the inverse quaternion Q^{-1}.
        %  \note The quaternion need not be unit.
        function invQ = quatInv(Q)
            
            invQ = zeros(size(Q));

            n = (Q'*Q);
            if (n < 1e-16)
                invQ = zeros(4,1);
                return;
            end

            invQ(1) = Q(1);
            invQ(2:4) = -Q(2:4);

            invQ = invQ / n;

        end
        
        %% quatLog
        %  Returns the logarithm map of a unit quaternion.
        %  @param[in] Q: Unit quaterion as a 4x1 vector [u; v], with u the scalar and v the vector part.
        %  @return: the quaternion logarithm log(Q).
        function log_Q = quatLog(Q)

            n = Q(1);
            e = Q(2:4);
            norm_e = norm(e);

            if (norm_e > 1e-16)
                theta = 2*real(atan2(norm_e,n));
                log_Q = theta*e/norm_e;
            else
                log_Q = zeros(size(e));
            end

        end
        
        %% quatExp
        %  Returns the quaternion exponential map.
        %  @param[in] v: 3x1 vector.
        %  @return: the quaternion exponential exp(log_Q).
        function Q = quatExp(v)

            norm_v_rot = norm(v);
            theta = norm_v_rot;

            Q = zeros(4,1);

            if (norm_v_rot > 1e-16)
            Q(1) = cos(theta/2);
            Q(2:4) = sin(theta/2)*v/norm_v_rot;
            else
              Q = [1 0 0 0]';
            end

        end

        %% quatDiff
        %  Returns the quaternion difference between two quaternions 
        %  as Q1*(Q2^{-1})
        %  @param[in] Q1: First quaterion as a 4x1 vector.
        %  @param[in] Q2: Second quaterion as a 4x1 vector.
        %  @return: the quaternion difference Q1*(Q2^{-1}).
        %  \note All quaternions must have the form Q = [u; v] where u is
        %  the scalar and v the 3x1 vector part.
        function Qdiff = quatDiff(Q1, Q2)
            
            if ( dot(Q1,Q2) < 0), Q2 = -Q2; end
            Qdiff = gmp_.quatProd(Q1, gmp_.quatInv(Q2));

        end
        
        %% quatLogDiff
        % Returns the quaternion logarithm of the difference between two 
        % quaternions as log( Q1*(Q2^{-1}) )
        function v = quatLogDiff(Q1, Q2)

            v = gmp_.quatLog(gmp_.quatDiff(Q1,Q2));

        end
        
        %% rotVel_to_qLogDot
        %  Returns the time derivative of the quaternion logarithm of a given 
        %  orientation given a rotational velocity in Cartesian space.
        %  @param[in] rotVel: rotational velocity in Cartesian space.
        %  @param[in] Q1: an orientation as unit quaternion.
        %  @return: time derivative of the quaternion logarithm of Q1.
        function qdot = rotVel_to_qLogDot(rotVel, Q1)
            
            w = Q1(1); % quaternion scalar part
            if ( (1-abs(w)) <= 1e-15)
                qdot = rotVel;
                return;
            end

            v = Q1(2:4); % quaternion vector part
            norm_v = norm(v);
            k = v / norm_v; % axis of rotation
            s_th = norm_v; % sin(theta/2)
            c_th = w;      % cos(theta/2)
            th = atan2(s_th, c_th); % theta/2
            Pk_rotVel = k*dot(k,rotVel); % project rotVel on k
            qdot = Pk_rotVel + (rotVel - Pk_rotVel)*th*c_th/s_th - th*cross(k,rotVel);
            % k_cross = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
            % inv_Jn = k*k' + (eye(3,3) - k*k')*th2*c_th2/s_th2 - th2*k_cross;
            % qdot = inv_Jn * rotVel;

        end
        
        %% qLogDot_to_rotVel
        %  Returns the Cartesian rotational velocity given a unit quaternion
        %  and its quaternion logarithm time derivative.
        %  @param[in] qdot: time derivative of quaternion logarithm of Q1.
        %  @param[in] Q1: orientation as unit quaternion.
        %  @return: Cartesian rotational velocity.
        function rotVel = qLogDot_to_rotVel(qdot, Q1)
            
            w = Q1(1); % quaternion scalar part
            if ( (1-abs(w)) <= 1e-15)
                rotVel = qdot;
                return;
            end

            v = Q1(2:4); % quaternion vector part
            norm_v = norm(v);
            k = v / norm_v; % axis of rotation
            s_th = norm_v; % sin(theta/2)
            c_th = w;      % cos(theta/2)
            th = atan2(s_th, c_th); % theta/2
            % k_cross = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
            % Jn = k*k' + (eye(3,3) - k*k')*s_th2*c_th2/th2 + (s_th2^2/th2)*k_cross;
            Pk_qdot = k*dot(k,qdot); % project rotVel on k
            rotVel = Pk_qdot + (qdot - Pk_qdot)*s_th*c_th/th + (s_th^2/th)*cross(k,qdot);
            % rotVel = Jn * qdot;

        end

        %% rotAccel_to_qLogDDot
        %  Returns the 2nd time derivative of the quaternion logarithm of a given 
        %  orientation given a rotational velocity and acceleration in Cartesian space.
        %  @param[in] rotAccel: rotational acceleration in Cartesian space.
        %  @param[in] rotVel: rotational velocity in Cartesian space.
        %  @param[in] Q1: an orientation as unit quaternion.
        %  @return: 2nd time derivative of the quaternion logarithm of Q1.
        function qddot = rotAccel_to_qLogDDot(rotAccel, rotVel, Q1)
            
            w = Q1(1); % quaternion scalar part
            if ( (1-abs(w)) <= 1e-15)
                qddot = rotAccel;
                return;
            end

            v = Q1(2:4); % quaternion vector part
            norm_v = norm(v);
            k = v / norm_v; % axis of rotation
            s_th = norm_v; % sin(theta/2)
            c_th = w;      % cos(theta/2)
            th = atan2(s_th, c_th); % theta/2
            
            th_cth_over_sth = th*c_th/s_th;
            
            Pk_rotVel = k*dot(k,rotVel); % projection of rotVel on k
            rotVel_bot = rotVel - Pk_rotVel;
            qdot = Pk_rotVel + rotVel_bot*th_cth_over_sth - th*cross(k,rotVel);
            % qdot = gmp_.rotVel_to_qLogDot(rotVel, Q1);
            
            k_dot = 0.5*(qdot - dot(k,qdot)*k)/th;
            th2_dot = 0.5*dot(k,qdot);
            
            invJnDot_rotVel = (1 - th_cth_over_sth)*(dot(k,rotVel)*k_dot + dot(k_dot,rotVel)*k) + ...
                    -th*cross(k_dot,rotVel) - th2_dot*cross(k,rotVel) ...
                    +( (s_th*c_th - th) / s_th^2 )*th2_dot*rotVel_bot;
                
            Pk_rotAccel = k*dot(k,rotAccel); % project rotAccel on k
            invJn_rotAccel = Pk_rotAccel + (rotAccel - Pk_rotAccel)*th_cth_over_sth - th*cross(k,rotAccel); 
            %invJn_rotAccel = gmp_.rotVel_to_qLogDot(rotAccel, Q1);
            
            qddot = invJn_rotAccel + invJnDot_rotVel;

        end

        %% qLogDDot_to_rotAccel
        % Returns the Cartesian rotational acceleration given a unit quaternion, its 
        % quaternion logarithm 2nd time derivative and a Cartesian rotational velocity.
        %  @param[in] qddot: 2nd time derivative of quaternion logarithm of Q1.
        %  @param[in] rotVel: rotational velocity in Cartesian space.
        %  @param[in] Q1: orientation as unit quaternion.
        %  @return: Cartesian rotational acceleration.
        function rotAccel = qLogDDot_to_rotAccel(qddot, rotVel, Q1)

            w = Q1(1); % quaternion scalar part
            if ( (1-abs(w)) <= 1e-15)
                rotAccel = qddot;
                return;
            end

            v = Q1(2:4); % quaternion vector part
            norm_v = norm(v);
            k = v / norm_v; % axis of rotation
            s_th = norm_v; % sin(theta/2)
            c_th = w;      % cos(theta/2)
            th = atan2(s_th, c_th); % theta/2
            
            Pk_rotVel = k*dot(k,rotVel); % projection of rotVel on k
            qdot = Pk_rotVel + (rotVel - Pk_rotVel)*th*c_th/s_th - th*cross(k,rotVel);
            % qdot = gmp_.rotVel_to_qLogDot(rotVel, Q1);
            
            qdot_bot = qdot - dot(k,qdot)*k; % projection of qdot on plane normal to k
            k_dot = 0.5*qdot_bot/th;
            th2_dot = 0.5*dot(k,qdot);
            
            sc_over_th = (s_th * c_th) / th;
            
            JnDot_qdot = (1 - sc_over_th)*(dot(k,qdot)*k_dot + dot(k_dot,qdot)*k) + ...
                    (s_th^2/th)*cross(k_dot,qdot) + ...
                    ( (1 - 2*s_th^2)/th - sc_over_th/th )*th2_dot*qdot_bot + ...
                    (2*sc_over_th - (s_th/th)^2)*th2_dot*cross(k,qdot);
            
            Pk_qddot = dot(k,qddot)*k; % projection of qddot on k
            Jn_qddot = Pk_qddot + (qddot - Pk_qddot)*sc_over_th + (s_th^2/th)*cross(k,qddot);
            %Jn_qddot = gmp_.qLogDot_to_rotVel(qddot, Q1);
            
            rotAccel = Jn_qddot + JnDot_qdot;

        end
        
        %% torque_to_qLogTorque
        %  Converts a Cartesian torque in the quaternion logarithm space.
        %  @param[in] torque: torque in Cartesian space.
        %  @param[in] Q1: an orientation as unit quaternion.
        %  @return: torque in the quaternion logarithm space.
        function logTorq = torque_to_qLogTorque(torque, Q1)
            
            w = Q1(1); % quaternion scalar part
            if ( (1-abs(w)) <= 1e-15)
                logTorq = torque;
                return;
            end
            
            v = Q1(2:4);
            norm_v = norm(v);
            k = v / norm_v;
            s_th = norm_v;
            c_th = w;
            th = atan2(s_th, c_th);
            Pk_torq = k*dot(k,torque); % project torque on k
            logTorq = Pk_torq + (torque - Pk_torq)*s_th*c_th/th - (s_th^2/th)*cross(k,torque);

        end
        
        %% qLogTorque_to_torque
        %  Converts a torque from the quaternion logarithm to the Cartesian space.
        %  @param[in] logTorq: torque in the quaternion logarithm space.
        %  @param[in] Q1: an orientation as unit quaternion.
        %  @return: torque in Cartesian space.
        function torque = qLogTorque_to_torque(logTorq, Q1)
           
            w = Q1(1); % quaternion scalar part
            if ( (1-abs(w)) <= 1e-15)
                torque = logTorq;
                return;
            end

            v = Q1(2:4); % quaternion vector part
            norm_v = norm(v);
            k = v / norm_v; % axis of rotation
            s_th = norm_v; % sin(theta/2)
            c_th = w;      % cos(theta/2)
            th = atan2(s_th, c_th); % theta/2
            Pk_logTorq = k*dot(k,logTorq); % project logTorq on k
            torque = Pk_logTorq + (logTorq - Pk_logTorq)*th*c_th/s_th + th*cross(k,logTorq);
            
        end
        
    end
    
    %% ======================================
    %% ==============   IO   ================
    %% ======================================
    methods (Access = public, Static)
        
        %% write
        %  Write the GMP model to a file.
        %  @param[in] gmp: Pointer to a @GMP or @GMPo object.
        %  @param[in] fid: Filename string or object of type @FileIO associated with the file.
        %  @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
        function write(gmp, fid, prefix)
            
            if (nargin < 3), prefix=''; end
            
            c = class(gmp);
            if ( strcmpi(c, 'GMP') ), gmp_.writeGMP(gmp, fid, prefix)
            elseif ( strcmpi(c, 'GMPo') ), gmp_.writeGMPo(gmp, fid, prefix)
            else, error('[gmp_::write]: Unsupported class type');
            end
            
        end
        
        %% read
        %  Reads the GMP model from a file.
        %  @param[in] gmp: Pointer to a @GMP or @GMPo object.
        %  @param[in] fid: Filename string or object of type @FileIO associated with the file.
        %  @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
        function read(gmp, fid, prefix)
            
            if (nargin < 3), prefix=''; end
            
            c = class(gmp);
            if ( strcmpi(c, 'GMP') ), gmp_.readGMP(gmp, fid, prefix)
            elseif ( strcmpi(c, 'GMPo') ), gmp_.readGMPo(gmp, fid, prefix)
            else, error('[gmp_::read]: Unsupported class type');
            end
            
        end
        
    end
    
    methods (Access = private, Static)
        
        %% writeGMP
        %  Write the GMP model to a file.
        %  @param[in] gmp: Pointer to a @GMP object.
        %  @param[in] fid: Filename string or object of type @FileIO associated with the file.
        %  @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
        function writeGMP(gmp, fid, prefix)
            
            if (nargin < 3), prefix=''; end
            
            if (ischar(fid))
               filename = fid;
               fid = FileIO(filename, bitor(FileIO.out,FileIO.trunc));
            end
            
            N_kernels = gmp.numOfKernels();
            n_dofs = gmp.numOfDoFs();
            
            fid.write([prefix 'weights'], gmp.W);
            fid.write([prefix 'damping'], gmp.D);
            fid.write([prefix 'stiffness'], gmp.K);
            fid.write([prefix 'N_kernels'], uint32(N_kernels) );
            fid.write([prefix 'N_DoFs'], uint32(n_dofs) );
            fid.write([prefix 'scale_type'], int32(gmp.traj_sc.getScaleType()) );
            fid.write([prefix 'c'], gmp.getCenters());
            fid.write([prefix 'h'], gmp.getInvWidths());
            
        end
        
        %% readGMP
        %  Reads the GMP model from a file.
        %  @param[in] gmp: Pointer to a @GMP object.
        %  @param[in] fid: Filename string or object of type @FileIO associated with the file.
        %  @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
        function readGMP(gmp, fid, prefix)
            
            if (nargin < 3), prefix=''; end
            
            if (ischar(fid))
               filename = fid;
               fid = FileIO(filename, FileIO.in);
            end
            
            gmp.W = fid.read([prefix 'weights']);
            gmp.W0 = gmp.W;
            gmp.D = fid.read([prefix 'damping']);
            gmp.K = fid.read([prefix 'stiffness']);
            %N_kernels = fid.read([prefix 'N_kernels']);
            %n_dofs = fid.read([prefix 'N_DoFs']);
            scale_type = fid.read([prefix 'scale_type']);
            c = fid.read([prefix 'c']);
            h = fid.read([prefix 'h']);
            gmp.setKernels2(c,h);
            
            gmp.Y0d = gmp.W*gmp.regressVec(0);
            gmp.Ygd = gmp.W*gmp.regressVec(1);
            
            gmp.setY0(gmp.Y0d);
            gmp.setGoal(gmp.Ygd);
            
            n_dofs = gmp.numOfDoFs();
            gmp.y_dot = zeros(n_dofs,1);
            gmp.z_dot = zeros(n_dofs,1);
            
            gmp.gmp_up = GMP_Update(gmp);
            gmp.gmp_up.initExpSigmaw(0.01);
            gmp.gmp_up.recursiveUpdate(false);
            
            if (scale_type == TrajScale.PROP_SCALE), gmp.setScaleMethod(TrajScale_Prop(n_dofs));
            elseif (scale_type == TrajScale.ROT_MIN_SCALE), gmp.setScaleMethod(TrajScale_Rot_min());
            elseif (scale_type == TrajScale.ROT_WB_SCALE), gmp.setScaleMethod(TrajScale_Rot_wb());
                elseif (scale_type == TrajScale.NONE), gmp.setScaleMethod(TrajScale_None(n_dofs));
            else, error(['[GMP_IO::read]: Unsupported scale type ''' num2str(scale_type) '''...']);
            end
            
        end
        
        %% writeGMPo
        %  Write the GMP model to a file.
        %  @param[in] gmp: Pointer to a @GMPo object.
        %  @param[in] fid: Filename string or object of type @FileIO associated with the file.
        %  @param[in] prefix: The prefix that will be used for writing the names of all GMP params (optional, default="").
        function writeGMPo(gmp, fid, prefix)
            
            if (nargin < 3), prefix=''; end
            
            if (ischar(fid))
                filename = fid;
                fid = FileIO(filename, bitor(FileIO.out,FileIO.trunc));
                gmp_.writeGMP(gmp, fid, prefix);
                fid.write([prefix 'Qd0'], gmp.Qd0);
            else
                gmp_.writeGMP(gmp, fid, prefix);
            end
            
        end
        
        %% readGMPo
        %  Reads the GMP model from a file.
        %  @param[in] gmp: Pointer to a @GMPo object.
        %  @param[in] fid: Filename string or object of type @FileIO associated with the file.
        %  @param[in] prefix: The prefix that will be used for reading the names of all GMP params (optional, default="").
        function readGMPo(gmp, fid, prefix)
            
            if (nargin < 3), prefix=''; end
            
            if (ischar(fid))
                filename = fid;
                fid = FileIO(filename, FileIO.in);
                gmp_.readGMP(gmp, fid, prefix);
                gmp.Qd0 = fid.read([prefix 'Qd0']);
            else
                gmp_.readGMP(gmp, fid, prefix);
            end
            
        end
        
    end

end
