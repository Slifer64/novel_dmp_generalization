%% TrajScale_Prop class
%
% For generalizing an n-DoF trajectory to new target/initial positions using 
% proportional scaling (as in the original DMP):
% 
% scaling matrix: Ks = diag( (g - y0) ./ (gd - yd0) )
% 
% where g, y0 are the new target and initial positions
%      gd, yd0 are the demonstrated target and initial positions
%

classdef TrajScale_Prop < TrajScale
    
    methods (Access = public)
        
        %% Constructor.
        %  @param[in] n_dof: degrees of freedom.
        function this = TrajScale_Prop(n_dof)
    
            this@TrajScale(n_dof);
            
        end

    end
    
    methods (Access = public) % Abstract implementations
        
        function scale_type = getScaleType(this)
            
            scale_type = TrajScale.PROP_SCALE;
            
        end
        
    end
    
    methods (Access = protected) % Abstract implementations

        function sc = calcScaling(this)

            sc = diag( (this.Yg - this.Y0) ./ (this.Ygd - this.Y0d) );
            
        end
        
        % ------------------------------------------
        
        function inv_sc = calcInvScaling(this)
            
            inv_sc = diag( (this.Ygd - this.Y0d) ./ (this.Yg - this.Y0) );
            
        end
        
    end
    
end
