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

classdef TrajScale_Rot_min < TrajScale
    
    methods (Access = public)
        
        %% Constructor.
        %  @param[in] n_dof: degrees of freedom.
        function this = TrajScale_Rot_min()
    
            this@TrajScale(3);
            
        end

    end
    
    methods (Access = public) % Abstract implementations
        
        function scale_type = getScaleType(this)
            
            scale_type = TrajScale.ROT_MIN_SCALE;
            
        end
        
    end
    
    methods (Access = protected)
       
        % ------------------------------------------
        
        function sc = calcScaling(this)
        
            nd = this.Ygd - this.Y0d;  nd = nd/norm(nd);
            n = this.Yg - this.Y0;  n = n/norm(n);
            dot_n_nd = dot(n,nd);
            if (abs(abs(dot_n_nd) - 1) < 1e-14)
                R = eye(3,3);
            else
                k = cross(nd,n);
                theta = acos(dot_n_nd);
                R = axang2rotm([k' theta]);
            end
            
            sc = R*norm(this.Yg - this.Y0)/norm(this.Ygd - this.Y0d);
            
        end
        
        function sc = calcInvScaling(this)
        
            nd = this.Ygd - this.Y0d;  nd = nd/norm(nd);
            n = this.Yg - this.Y0;  n = n/norm(n);
            dot_n_nd = dot(n,nd);
            if (abs(abs(dot_n_nd) - 1) < 1e-14)
                R = eye(3,3);
            else
                k = cross(nd,n);
                theta = acos(dot_n_nd);
                R = axang2rotm([k' theta]);
            end
            sc = R'*norm(this.Ygd - this.Y0d)/norm(this.Yg - this.Y0);
            
        end

        % ------------------------------------------
        
    end
    
end
