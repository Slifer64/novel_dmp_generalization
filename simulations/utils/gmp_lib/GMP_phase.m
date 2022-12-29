%% GMP phase class
%  Generalized movement primitive phase variable.
%

classdef GMP_phase < matlab.mixin.Copyable
    
    methods (Access = public)
        
        %% GMP constructor.
        %  @param[in] x: phase variable.
        %  @param[in] x_dot: phase variable 1st time derivative.
        %  @param[in] x_ddot: phase variable 2nd time derivative.
        function this = GMP_phase(x, x_dot, x_ddot)
                
            this.x = x;
            this.x_dot = x_dot;
            this.x_ddot = x_ddot;
            
        end

    end
    
    properties (Access = public)

        x % phase variable
        x_dot % phase variable 1st time derivative
        x_ddot % phase variable 2nd time derivative

    end

    
end
