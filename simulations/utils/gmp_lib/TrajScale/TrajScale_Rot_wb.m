%% TrajScale_Rot_wb class
%
% For generalizing a 3-DoF trajectory to new target/initial positions using 
% a rotation and scaling base on
% scale type 3 from 
% 'A novel DMP formulation for global and frame independent spatial scaling in the task space'
% DOI: 10.1109/RO-MAN47096.2020.9223500
% 
% scaling matrix: 
% Ks = ( ||g - y0|| / ||gd - yd0|| ) * {rotation matrix that aligns gd - yd0 with g - y0 and the rotation axis is the closest to the workbench normal}
% 
% where g, y0 are the new target and initial positions
%      gd, yd0 are the demonstrated target and initial positions
%

classdef TrajScale_Rot_wb < TrajScale
   
    methods (Access = public)
        
        %% Constructor.
        %  @param[in] n_dof: degrees of freedom.
        function this = TrajScale_Rot_wb()
    
            this@TrajScale(3);
            
            this.setWorkBenchNormal([0; 0; 1]);
            
        end
        
        %% Sets the normal to the work bench vector. Applies only to 'ROT_WB_SCALE'.
        %  @param[in] n_wb: 3x1 vector of the normal to the workbench.
        function setWorkBenchNormal(this, n_wb)
            
            this.n_wb = n_wb;
            
            this.calcScaling();
            this.calcInvScaling();
            
        end
        

    end
    
    methods (Access = public) % Abstract implementations
        
        function scale_type = getScaleType(this)
            
            scale_type = TrajScale.ROT_WB_SCALE;
            
        end
        
    end
    
    methods (Access = protected)
       
        % ------------------------------------------
        
        function sc = calcScaling(this)
        
            nd = this.Ygd - this.Y0d;  nd = nd/norm(nd);
            n = this.Yg - this.Y0;  n = n/norm(n);
            dot_n_nd = dot(n,nd);
            if (abs(abs(dot_n_nd) - 1) < 1e-12)
                R = eye(3,3);
            else
                p = (n - nd) / norm(n - nd);
                k = this.n_wb - dot(p,this.n_wb)*p;
                k = k / norm(k);
                
                % Ik = eye(3,3) - k*k';
                nd2 = nd - dot(k,nd)*k; % Ik*nd
                n2 = n - dot(k,n)*k; % Ik*n
                cos_n = dot(n2,nd2) / (norm(n2)*norm(nd2));
                if ( abs(abs(cos_n)-1) < 1e-15), theta = pi;
                else, theta = acos( cos_n );
                end

                R = axang2rotm([k' theta]);
                if (norm(n - R*nd) > 1e-8), R = axang2rotm([k' -theta]); end

%                 k = this.calcAxis(nd, n);
%                 [sth, cth] = this.calcAngle(nd, n, k);
%                 theta = atan2(sth, cth);
%                 R = axang2rotm([k' theta]);
                
                if (norm(n - R*nd) > 1e-6)
                   error('[TrajScale_Rot_wb::calcScaling]: R produces significant error...') ;
                end
            end
            
            sc = R*norm(this.Yg - this.Y0)/norm(this.Ygd - this.Y0d);
            
        end
        
        function sc = calcInvScaling(this)
        
            sc = inv(this.calcScaling());
            
        end

        % ------------------------------------------
        
        function axis = calcAxis(this, v1, v2)

            n1 = v1 / norm(v1);
            n2 = v2 / norm(v2);

            nAvg = (n1 + n2) / 2;
            nCross = cross(n1, n2);

            nPlane = cross(nAvg, nCross);
            nPlane = nPlane / norm(nPlane);

            if nPlane(3) < 0
                nPlane = -nPlane;
            end

            if nPlane(1)^2 > 0 || nPlane(2)^2 > 0

                nPlane2 = cross([0; 0; 1], nPlane);

                if nPlane(1)^2 > 0

                    num = (nPlane2(1) * nPlane(3) - nPlane2(3) * nPlane(1));
                    denom = (nPlane2(2) * nPlane(1) - nPlane2(1) * nPlane(2));
                    q = num / denom;

                    axis = [-(nPlane(3) + nPlane(2) * q) / nPlane(1); q; 1];

                else

                    axis = [(nPlane2(2) * nPlane(3) - nPlane2(3) * nPlane(2)) / (nPlane2(1) * nPlane(2)); -nPlane(3) / nPlane(2)];

                end

                axis = axis / norm(axis);

            else

                axis = [0; 0; 1];

            end

        end
        
        function [sth, cth] = calcAngle(this, v1, v2, axis )

            n1 = v1 / norm(v1);
            n2 = v2 / norm(v2);

            K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
            K1 = K * n1;
            K2 = K * K * n1;
            d = n2 - n1;

            if K1(1)^2 > 0

                if (-K1(2) * K2(1) + K1(1) * K2(2))^2 > 0

                    vth = (K1(1) * d(2) - d(1) * K1(2))  / (-K1(2) * K2(1) + K1(1) * K2(2));

                else

                    vth = (K1(1) * d(3) - d(1) * K1(3))  / (-K1(3) * K2(1) + K1(1) * K2(3));        

                end
                sth = (-K2(1) * vth + d(1)) / K1(1);

            elseif K1(2) ^2 > 0

                 if (-K1(1) * K2(2) + K1(2) * K2(1))^2 > 0

                    vth = (K1(2) * d(1) - d(2) * K1(1))  / (-K1(1) * K2(2) + K1(2) * K2(1));

                else

                    vth = (K1(2) * d(3) - d(2) * K1(3))  / (-K1(3) * K2(2) + K1(2) * K2(3));        

                 end
                 sth = (-K2(2) * vth + d(2)) / K1(2);

            else

                if (K1(1) * K2(3) + K1(3) * K2(1))^2 > 0

                    vth = (K1(3) * d(1) - d(3) * K1(1))  / (-K1(1) * K2(3) + K1(3) * K2(1));

                else

                    vth = (K1(3) * d(2) - d(3) * K1(2))  / (-K1(2) * K2(3) + K1(3) * K2(2));        

                end
                sth = (-K2(3) * vth + d(3)) / K1(3);

            end

            cth = 1 - vth;

        end

    end
    
    properties (Access = private)

        n_wb % workbench normal

    end
    
end
