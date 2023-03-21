%% TrajScale_Rot_min class
%
% For generalizing a 3-DoF trajectory to new target/initial positions using 
% a rotation and scaling base on
% scale type 1 from 
% 'A novel DMP formulation for global and frame independent spatial scaling in the task space'
% DOI: 10.1109/RO-MAN47096.2020.9223500
% 
% scaling matrix: 
% Ks = ( ||g - y0|| / ||gd - yd0|| ) * {rotation matrix that aligns gd - yd0 with g - y0 using the minimum angle of rotation}
% 
% where g, y0 are the new target and initial positions
%      gd, yd0 are the demonstrated target and initial positions
%

classdef TrajScale_roto_dial < TrajScale
    
    methods (Access = public)
        
        %% Constructor.
        %  @param[in] n_dof: degrees of freedom.
        function this = TrajScale_roto_dial(n_dofs)
    
            this@TrajScale(n_dofs);
            
        end

    end
    
    methods (Access = public) % Abstract implementations
        
        function scale_type = getScaleType(this)
            
            scale_type = TrajScale.ROT_DIAL_SCALE;
            
        end
        
    end
    
    methods (Access = protected)
       
        % ------------------------------------------
        
        function sc = calcScaling(this, Y0, Yg)
            
            this.Y0 = Y0;
            this.Yg = Yg;
            
            % Write the inputs as 1D array
            x0 = this.Ygd - this.Y0d;
            x1 = this.Yg - this.Y0;

            % Extract the norms
            norm_0 = norm(x0);
            norm_1 = norm(x1);

            % Normalize the two vectors
            x0_norm = x0 / norm_0;
            x1_norm = x1 / norm_1;
            Mx0 = this.fnAR(x0_norm);
            Mx1 = this.fnAR(x1_norm);
            M = Mx1' * Mx0;
            sc = M * norm_1 / norm_0;
            
        end
        
        function sc = calcInvScaling(this)
        
            sc = inv(this.calcScaling(this.Y0, this.Yg));
            
        end
        
        function R = fnAR(this, x)
        % Accelerated rotation of vector x to the direction of axis x_0.
            x = x(:);
            n = length(x);
            R = eye(n, n);
            step = 1;
            while (step < n)
                A = eye(n, n);
                it = 1;
                while(it < n - step + 1)
                    r2 = x(it) * x(it) + x(it + step) * x(it + step);
                    if (r2 > 0)
                        r = sqrt(r2);
                        pcos = x(it) / r;
                        psin = -x(it + step) / r;

                        % Base 2-dimensional rotation
                        A(it, it) = pcos;
                        A(it, it + step) = - psin;
                        A(it + step, it) = psin;
                        A(it + step, it + step) = pcos;
                    end
                    it = it + 2 * step;
                end
                step = step * 2;
                x = A * x;
                R = A * R;
            end

        end

        % ------------------------------------------
        
    end
    
end
