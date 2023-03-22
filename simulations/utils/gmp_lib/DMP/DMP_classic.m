% Implements the DMP from:
% "Dynamical Movement Primitives: Learning Attractor Models for Motor Behaviors"
% DOI: 10.1162/NECO_a_00393

classdef DMP_classic < DMP_
    
    methods (Access = public)
        
        function this = DMP_classic(n_dofs, N_kernels)
            
            this@DMP_(n_dofs, N_kernels, 'diag', 'gaussian', 'linear');
            
        end

        function f_g = goal_attractor(this, y, y_dot, tau)
           
            f_g = (this.K*(this.g - y) - tau*this.D*y_dot) / tau^2;
            
        end
        
        function f_s = shape_attractor(this, s, s_dot, s_ddot, tau)
           
            f_s = this.forcing_term(s) / tau^2;

        end
        
    end
    
    methods (Access = protected)
        
        function Fd = calcFd(this, s, y, dy, ddy, g, y0, tau)
            
            f_gat = this.gating_fun(s);
            Fd = (ddy*tau^2 + this.D*tau*dy - this.K*(g - y)) ./ f_gat;
            
        end 

    end
    
end