%% TrajScale_None class
%
% Produces a unit scaling (so essentially no scaling)
%

classdef TrajScale_None < TrajScale
    
    methods (Access = public)
        
        %% Constructor.
        %  @param[in] n_dof: degrees of freedom.
        function this = TrajScale_None(n_dof)
    
            this@TrajScale(n_dof);
            
        end

    end
    
    methods (Access = public) % Abstract implementations
        
        function scale_type = getScaleType(this)
            
            scale_type = TrajScale.NONE;
            
        end
        
    end
    
    methods (Access = protected) % Abstract implementations

        function sc = calcScaling(this, Y0, Yg)

            sc = eye(this.n_dof, this.n_dof);
            
        end
        
        % ------------------------------------------
        
        function inv_sc = calcInvScaling(this)
            
            inv_sc = eye(this.n_dof, this.n_dof);
            
        end
        
    end
    
end
